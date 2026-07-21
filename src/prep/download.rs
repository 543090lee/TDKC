use std::path::{Path, PathBuf};
use std::pin::Pin;
use std::sync::{Arc, Mutex};
use std::task::{Context, Poll};
use std::time::Duration;

use anyhow::{anyhow, Context as _, Result};
use async_compression::tokio::bufread::GzipDecoder;
use async_trait::async_trait;
use rustc_hash::FxHashMap;
use sha2::{Digest, Sha256};
use tokio::io::{AsyncBufReadExt, AsyncRead, AsyncReadExt, AsyncWriteExt, BufReader, BufWriter, ReadBuf};
use tokio_util::io::StreamReader;
use futures::StreamExt;


#[derive(Clone, Copy, Debug)]
pub enum BackendKind {
    Http,
}

#[async_trait]
pub trait Backend: Send + Sync {
    async fn fetch_range(&self, path: &str, start: u64) -> Result<reqwest::Response>;
}

pub struct HttpBackend {
    client: reqwest::Client,
}

impl HttpBackend {
    const BASE: &'static str = "https://ftp.ncbi.nlm.nih.gov";

    fn new() -> Result<Self> {
        let client = reqwest::Client::builder()
            .pool_max_idle_per_host(32)
            .connect_timeout(Duration::from_secs(30))
            .build()
            .context("building reqwest client")?;
        Ok(Self { client })
    }
}

#[async_trait]
impl Backend for HttpBackend {
    async fn fetch_range(&self, path: &str, start: u64) -> Result<reqwest::Response> {
        let url = format!("{}{}", Self::BASE, path);
        let mut req = self.client.get(&url);
        if start > 0 {
            req = req.header(reqwest::header::RANGE, format!("bytes={}-", start));
        }
        let resp = req.send().await.with_context(|| format!("GET {}", url))?;
        let status = resp.status();
        if start > 0 {
            if status != reqwest::StatusCode::PARTIAL_CONTENT {
                return Err(anyhow!("GET {} range from {} expected 206, got {}", url, start, status));
            }
        } else if !status.is_success() {
            return Err(anyhow!("GET {} returned status {}", url, status));
        }
        Ok(resp)
    }
}

pub fn make_backend(kind: BackendKind) -> Result<Arc<dyn Backend>> {
    Ok(match kind {
        BackendKind::Http => Arc::new(HttpBackend::new()?) as Arc<dyn Backend>,
    })
}

pub struct Sha256Reader<R> {
    inner: R,
    hasher: Arc<Mutex<Sha256>>,
    bytes: Arc<Mutex<u64>>,
}

impl<R: AsyncRead + Unpin> Sha256Reader<R> {
    pub fn new(inner: R) -> (Self, ShaHandle) {
        let hasher = Arc::new(Mutex::new(Sha256::new()));
        let bytes = Arc::new(Mutex::new(0u64));
        let reader = Self {
            inner,
            hasher: Arc::clone(&hasher),
            bytes: Arc::clone(&bytes),
        };
        let handle = ShaHandle { hasher, bytes };
        (reader, handle)
    }
}

impl<R: AsyncRead + Unpin> AsyncRead for Sha256Reader<R> {
    fn poll_read(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &mut ReadBuf<'_>,
    ) -> Poll<std::io::Result<()>> {
        let before = buf.filled().len();
        let res = Pin::new(&mut self.inner).poll_read(cx, buf);
        if let Poll::Ready(Ok(())) = &res {
            let after = buf.filled().len();
            if after > before {
                let new_bytes = &buf.filled()[before..after];
                if let Ok(mut hasher) = self.hasher.lock() {
                    hasher.update(new_bytes);
                }
                if let Ok(mut b) = self.bytes.lock() {
                    *b += (after - before) as u64;
                }
            }
        }
        res
    }
}

#[derive(Clone)]
pub struct ShaHandle {
    hasher: Arc<Mutex<Sha256>>,
    bytes: Arc<Mutex<u64>>,
}

impl ShaHandle {
    pub fn finalize(&self) -> Option<(String, u64)> {
        let hasher_clone = self.hasher.lock().ok()?.clone();
        let bytes = *self.bytes.lock().ok()?;
        Some((hex::encode(hasher_clone.finalize()), bytes))
    }
}

pub type BoxedReader = Pin<Box<dyn AsyncRead + Send + Unpin>>;

async fn reconnect(backend: &Arc<dyn Backend>, path: &str, offset: u64) -> Result<reqwest::Response> {
    with_retry(&format!("resume[{}]", path), 10, || backend.fetch_range(path, offset)).await
}

/// Adapt a possibly-interrupted HTTP body into a continuous byte stream: when
/// the connection drops mid-transfer, reconnect with a `Range` request from the
/// last received byte and keep going. If a total length is known, a short read
/// at the end is reported as an error so truncation never passes silently.
fn resume_stream(
    backend: Arc<dyn Backend>,
    path: String,
    initial: reqwest::Response,
    total: Option<u64>,
) -> impl futures::Stream<Item = std::io::Result<bytes::Bytes>> {
    // Give up if repeated reconnects make no forward progress, so a server that
    // keeps dropping the body at the same offset can't spin forever.
    const MAX_STALLED_RESUMES: u32 = 10;
    async_stream::stream! {
        let mut offset: u64 = 0;
        let mut resp = initial;
        let mut stalled: u32 = 0;
        loop {
            let start_offset = offset;
            let mut body = resp.bytes_stream();
            let mut interrupted = false;
            while let Some(item) = body.next().await {
                match item {
                    Ok(chunk) => {
                        offset += chunk.len() as u64;
                        yield Ok(chunk);
                    }
                    Err(_) => {
                        interrupted = true;
                        break;
                    }
                }
            }
            if !interrupted {
                if let Some(expected) = total {
                    if offset != expected {
                        yield Err(std::io::Error::other(format!(
                            "{}: truncated download, got {} of {} bytes", path, offset, expected)));
                    }
                }
                return;
            }
            if offset == start_offset {
                stalled += 1;
                if stalled >= MAX_STALLED_RESUMES {
                    yield Err(std::io::Error::other(format!(
                        "{}: no progress after {} resume attempts at byte {}", path, stalled, offset)));
                    return;
                }
            } else {
                stalled = 0;
            }
            match reconnect(&backend, &path, offset).await {
                Ok(next) => resp = next,
                Err(e) => {
                    yield Err(std::io::Error::other(format!("{}: resume failed: {:#}", path, e)));
                    return;
                }
            }
        }
    }
}

async fn open_resilient(backend: Arc<dyn Backend>, path: &str) -> Result<(BoxedReader, ShaHandle)> {
    let resp = backend.fetch_range(path, 0).await?;
    let total = resp.content_length();
    let stream = Box::pin(resume_stream(Arc::clone(&backend), path.to_string(), resp, total));
    let (rdr, h) = Sha256Reader::new(StreamReader::new(stream));
    Ok((Box::pin(rdr), h))
}

pub async fn stream_raw(backend: &Arc<dyn Backend>, path: &str) -> Result<(BoxedReader, ShaHandle)> {
    open_resilient(Arc::clone(backend), path).await
}

pub async fn stream_gz(backend: &Arc<dyn Backend>, path: &str) -> Result<(BoxedReader, ShaHandle)> {
    let (raw, h) = open_resilient(Arc::clone(backend), path).await?;
    let decoder = GzipDecoder::new(BufReader::new(raw));
    Ok((Box::pin(decoder), h))
}

pub async fn with_retry<F, Fut, T>(label: &str, max_attempts: u32, mut f: F) -> Result<T>
where
    F: FnMut() -> Fut,
    Fut: std::future::Future<Output = Result<T>>,
{
    let mut delay = Duration::from_secs(2);
    let mut last_err: Option<anyhow::Error> = None;
    for attempt in 1..=max_attempts {
        match f().await {
            Ok(v) => return Ok(v),
            Err(e) => {
                if attempt == max_attempts {
                    last_err = Some(e);
                    break;
                }
                eprintln!(
                    "[retry] {} attempt {}/{} failed: {:#}; sleeping {:?}",
                    label, attempt, max_attempts, e, delay
                );
                tokio::time::sleep(delay).await;
                delay = delay.mul_f32(1.5);
                last_err = Some(e);
            }
        }
    }
    Err(last_err.unwrap_or_else(|| anyhow!("retry: no error captured")))
}

pub struct TaxdumpPaths {
    pub nodes_dmp: PathBuf,
    pub names_dmp: PathBuf,
}

pub async fn download_taxdump(
    backend: &Arc<dyn Backend>,
    output_dir: &Path,
) -> Result<TaxdumpPaths> {
    let tax_dir = output_dir.join("taxonomy");
    tokio::fs::create_dir_all(&tax_dir).await?;

    let nodes = tax_dir.join("nodes.dmp");
    let names = tax_dir.join("names.dmp");
    if nodes.exists() && names.exists() {
        eprintln!("Taxonomy already present, skipping download.");
        return Ok(TaxdumpPaths { nodes_dmp: nodes, names_dmp: names });
    }

    eprintln!("Downloading taxonomy information ...");

    let buf = with_retry("taxdump", 5, || async {
        let (mut rdr, sha) = stream_raw(backend, "/pub/taxonomy/taxdump.tar.gz").await?;
        let mut buf = Vec::with_capacity(80 * 1024 * 1024);
        rdr.read_to_end(&mut buf).await
            .context("reading taxdump body")?;
        let _ = sha.finalize();
        Ok(buf)
    })
    .await?;

    let tax_dir_owned = tax_dir.clone();
    let (nodes, names) = tokio::task::spawn_blocking(move || -> Result<(PathBuf, PathBuf)> {
        let gz = flate2::read::GzDecoder::new(std::io::Cursor::new(buf));
        let mut archive = tar::Archive::new(gz);
        let mut nodes_path = None;
        let mut names_path = None;

        for entry in archive.entries()? {
            let mut entry = entry?;
            let path_in_tar = entry.path()?.to_path_buf();
            let fname = match path_in_tar.file_name().and_then(|s| s.to_str()) {
                Some(n) => n.to_string(),
                None => continue,
            };

            match fname.as_str() {
                "nodes.dmp" => {
                    let out = tax_dir_owned.join("nodes.dmp");
                    let mut f = std::fs::File::create(&out)?;
                    std::io::copy(&mut entry, &mut f)?;
                    nodes_path = Some(out);
                }
                "names.dmp" => {
                    let out = tax_dir_owned.join("names.dmp");
                    let mut f = std::fs::File::create(&out)?;
                    std::io::copy(&mut entry, &mut f)?;
                    names_path = Some(out);
                }
                _ => {} 
            }
        }
        Ok((
            nodes_path.ok_or_else(|| anyhow!("nodes.dmp not in taxdump.tar.gz"))?,
            names_path.ok_or_else(|| anyhow!("names.dmp not in taxdump.tar.gz"))?,
        ))
    })
    .await
    .context("taxdump extraction task")??;

    Ok(TaxdumpPaths { nodes_dmp: nodes, names_dmp: names })
}

#[derive(Debug, Clone)]
pub struct GenomeEntry {
    pub assembly_accession: String,
    pub taxid: u32,
    pub asm_name: String,
    pub ftp_path: String, // the format is like https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../GCF_xxx_yyy
}

impl GenomeEntry {
    
    pub fn fna_gz_path(&self) -> Result<String> {
        // ftp_path is something like https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2
        let path_only = strip_ncbi_host(&self.ftp_path)
            .ok_or_else(|| anyhow!("unrecognized ftp_path: {}", self.ftp_path))?;

        let trimmed_path = path_only.trim_end_matches('/');
        let basename = trimmed_path.rsplit('/').next().unwrap_or("");
        if basename.is_empty() {
            return Err(anyhow!("empty basename in ftp_path: {}", self.ftp_path));
        }
        Ok(format!("{}/{}_genomic.fna.gz", path_only, basename))
    }
}

fn strip_ncbi_host(url: &str) -> Option<&str> {
    for prefix in [
        "https://ftp.ncbi.nlm.nih.gov",
        "http://ftp.ncbi.nlm.nih.gov",
        "ftp://ftp.ncbi.nlm.nih.gov",
    ] {
        if let Some(rest) = url.strip_prefix(prefix) {
            return Some(rest);
        }
    }
    None
}

pub async fn fetch_assembly_summary(
    backend: &Arc<dyn Backend>,
    domain: &str,
) -> Result<Vec<GenomeEntry>> {
    let remote_dir = if domain == "human" {
        "vertebrate_mammalian/Homo_sapiens".to_string()
    } else {
        domain.to_string()
    };
    let path = format!("/genomes/refseq/{}/assembly_summary.txt", remote_dir);

    let body = with_retry(&format!("assembly_summary[{}]", domain), 5, || async {
        let (mut rdr, _sha) = stream_raw(backend, &path).await?;
        let mut s = String::with_capacity(8 * 1024 * 1024);
        rdr.read_to_string(&mut s).await
            .context("reading assembly_summary body")?;
        Ok(s)
    })
    .await?;

    let mut out = Vec::new();
    for line in body.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 20 {
            continue;
        }
        let acc = cols[0].trim();
        let taxid: u32 = match cols[5].trim().parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let assembly_level = cols[11].trim();
        let asm_name = cols[15].trim();
        let ftp_path = cols[19].trim();

        // Match Kraken2 : keep only Complete Genome /Chromosome assemblies
        if assembly_level != "Complete Genome" && assembly_level != "Chromosome" {
            continue;
        }
        // Match Kraken2 : the human library keeps only lines mentioning the GRC submitters
        if domain == "human" && !line.contains("Genome Reference Consortium") {
            continue;
        }
        if ftp_path.is_empty() || ftp_path == "na" {
            continue;
        }

        out.push(GenomeEntry {
            assembly_accession: acc.to_string(),
            taxid,
            asm_name: asm_name.to_string(),
            ftp_path: ftp_path.to_string(),
        });
    }
    Ok(out)
}

const ACC2TAXID_SOURCES: [&str; 2] =
    ["nucl_gb.accession2taxid.gz", "nucl_wgs.accession2taxid.gz"];

pub async fn prefetch_acc2taxid(backend: &Arc<dyn Backend>, tax_dir: &Path) -> Result<()> {
    for src in ACC2TAXID_SOURCES {
        let local = tax_dir.join(src.trim_end_matches(".gz"));
        if !local.exists() {
            let remote = format!("/pub/taxonomy/accession2taxid/{}", src);
            ensure_acc2taxid_on_disk(backend, &remote, &local).await?;
        }
    }
    Ok(())
}

pub async fn fetch_acc2taxid_filtered(
    backend: &Arc<dyn Backend>,
    wanted: &std::collections::HashSet<String>,
    tax_dir: &Path,
) -> Result<FxHashMap<String, u32>> {
    let mut out: FxHashMap<String, u32> = FxHashMap::default();
    if wanted.is_empty() {
        return Ok(out);
    }

    for src in ACC2TAXID_SOURCES {
        if out.len() >= wanted.len() {
            break;
        }
        let remote = format!("/pub/taxonomy/accession2taxid/{}", src);
        let local = tax_dir.join(src.trim_end_matches(".gz"));

        if !local.exists() {
            ensure_acc2taxid_on_disk(backend, &remote, &local).await?;
        }

        let f = tokio::fs::File::open(&local).await
            .with_context(|| format!("open acc2taxid {}", local.display()))?;
        let new_entries = stream_acc2taxid_into(f, wanted).await?;

        for (k, v) in new_entries {
            out.entry(k).or_insert(v);
        }
    }
    Ok(out)
}

async fn ensure_acc2taxid_on_disk(
    backend: &Arc<dyn Backend>,
    remote_path: &str,
    local_path: &Path,
) -> Result<()> {
    if let Some(parent) = local_path.parent() {
        tokio::fs::create_dir_all(parent).await?;
    }
    let tmp = local_path.with_extension("partial");
    let tmp_owned = tmp.clone();
    with_retry(&format!("acc2taxid-download[{}]", remote_path), 5, || {
        let tmp_owned = tmp_owned.clone();
        async move {
            let (mut rdr, _sha) = stream_gz(backend, remote_path).await?;
            let out_file = tokio::fs::File::create(&tmp_owned).await
                .with_context(|| format!("create {}", tmp_owned.display()))?;
            let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, out_file);
            tokio::io::copy(&mut rdr, &mut writer).await
                .with_context(|| format!("decompressing {}", remote_path))?;
            writer.flush().await?;
            Ok(())
        }
    })
    .await?;
    tokio::fs::rename(&tmp, local_path).await
        .with_context(|| format!("rename {} -> {}", tmp.display(), local_path.display()))?;
    Ok(())
}

async fn stream_acc2taxid_into<R>(
    rdr: R,
    wanted: &std::collections::HashSet<String>,
) -> Result<Vec<(String, u32)>>
where
    R: AsyncRead + Unpin,
{
    let mut buffered = BufReader::with_capacity(4 * 1024 * 1024, rdr);
    let mut line = String::new();
    let mut first = true;
    let mut hits = Vec::new();
    loop {
        line.clear();
        let n = tokio::io::AsyncBufReadExt::read_line(&mut buffered, &mut line).await?;
        if n == 0 {
            break;
        }
        if first {
            first = false;
            if line.starts_with("accession") {
                continue; // header
            }
        }
        let mut it = line.split('\t');
        let _acc = it.next();
        let acc_ver = match it.next() {
            Some(s) => s,
            None => continue,
        };
        let taxid_str = match it.next() {
            Some(s) => s.trim(),
            None => continue,
        };
        if !wanted.contains(acc_ver) {
            continue;
        }
        let taxid: u32 = match taxid_str.parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        hits.push((acc_ver.to_string(), taxid));
    }
    Ok(hits)
}

pub struct UnivecMetadata {
    pub fetch_path: String,
    pub filename: String,
    pub taxid: u32,
}

pub fn resolve_univec(domain: &str) -> UnivecMetadata {
    let filename = if domain == "univec" { "UniVec" } else { "UniVec_Core" };
    UnivecMetadata {
        fetch_path: format!("/pub/UniVec/{}", filename),
        filename: filename.to_string(),
        taxid: 28384, // I am following the same label as Kraken2, they get special label!
    }
}

/// Some RefSeq libraries (plasmid, plastid, mitochondrion) have no assembly_summary.txt
pub async fn list_refseq_catalog_files(
    backend: &Arc<dyn Backend>,
    dir: &str,
) -> Result<Vec<String>> {
    let path = format!("/genomes/refseq/{}/", dir);
    let body = with_retry(&format!("listing[{}]", dir), 5, || async {
        let (mut rdr, _sha) = stream_raw(backend, &path).await?;
        let mut s = String::with_capacity(64 * 1024);
        rdr.read_to_string(&mut s).await
            .context("reading directory listing")?;
        Ok(s)
    })
    .await?;

    let needle = ".genomic.fna.gz";
    let prefix = format!("{}.", dir);
    let bytes = body.as_bytes();
    let mut seen = std::collections::HashSet::new();
    let mut files = Vec::new();
    for (idx, _) in body.match_indices(needle) {
        let end = idx + needle.len();
        // Walk back to the start of the filename token
        let mut start = idx;
        while start > 0 {
            let c = bytes[start - 1];
            if c == b'"' || c == b'\'' || c == b'>' || c == b'/' {
                break;
            }
            start -= 1;
        }
        let name = &body[start..end];
        if name.starts_with(&prefix) && seen.insert(name.to_string()) {
            files.push(name.to_string());
        }
    }
    files.sort();
    if files.is_empty() {
        return Err(anyhow!("no {}*.genomic.fna.gz files found in /genomes/refseq/{}/", dir, dir));
    }
    Ok(files)
}

pub async fn fetch_and_decompress_to(
    backend: &Arc<dyn Backend>,
    remote_path: &str,
    local_path: &Path,
) -> Result<Vec<String>> {
    let local_owned = local_path.to_path_buf();
    let accs = with_retry(&format!("catalog[{}]", remote_path), 5, || {
        let local_owned = local_owned.clone();
        async move {
            let (rdr, _sha) = stream_gz(backend, remote_path).await?;
            let mut buffered = BufReader::with_capacity(4 * 1024 * 1024, rdr);
            let out_file = tokio::fs::File::create(&local_owned).await
                .with_context(|| format!("create {}", local_owned.display()))?;
            let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, out_file);

            let mut accs = Vec::new();
            let mut line = String::new();
            loop {
                line.clear();
                let n = buffered.read_line(&mut line).await
                    .with_context(|| format!("reading {}", remote_path))?;
                if n == 0 {
                    break;
                }
                if line.starts_with('>') {
                    if let Some(tok) = line[1..].split_whitespace().next() {
                        accs.push(tok.to_string());
                    }
                }
                writer.write_all(line.as_bytes()).await?;
            }
            writer.flush().await?;
            Ok(accs)
        }
    })
    .await?;
    Ok(accs)
}

pub async fn stream_local_raw(path: &str) -> Result<(BoxedReader, ShaHandle)> {
    let f = tokio::fs::File::open(path).await
        .with_context(|| format!("failed to open local file {}", path))?;
    let (rdr, h) = Sha256Reader::new(f);
    Ok((Box::pin(rdr), h))
}

pub async fn stream_local_gz(path: &str) -> Result<(BoxedReader, ShaHandle)> {
    let f = tokio::fs::File::open(path).await
        .with_context(|| format!("failed to open local file {}", path))?;
    let (rdr, h) = Sha256Reader::new(f);
    let buffered = BufReader::new(rdr);
    let decoder = GzipDecoder::new(buffered);
    Ok((Box::pin(decoder), h))
}
