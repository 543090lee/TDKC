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
use tokio::io::{AsyncRead, AsyncReadExt, BufReader, ReadBuf};
use tokio_util::io::StreamReader;


#[derive(Clone, Copy, Debug)]
pub enum BackendKind {
    Http,
}

#[async_trait]
pub trait Backend: Send + Sync {
    async fn fetch(&self, path: &str) -> Result<reqwest::Response>;
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
    async fn fetch(&self, path: &str) -> Result<reqwest::Response> {
        let url = format!("{}{}", Self::BASE, path);
        let resp = self
            .client
            .get(&url)
            .send()
            .await
            .with_context(|| format!("GET {}", url))?;
        if !resp.status().is_success() {
            return Err(anyhow!("GET {} returned status {}", url, resp.status()));
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

pub fn response_to_reader(
    resp: reqwest::Response,
) -> (Sha256Reader<StreamReader<ReqwestByteStream, bytes::Bytes>>, ShaHandle) {
    let stream = ReqwestByteStream::new(resp);
    let reader = StreamReader::new(stream);
    Sha256Reader::new(reader)
}

pub struct ReqwestByteStream {
    inner: Pin<
        Box<
            dyn futures::Stream<Item = std::result::Result<bytes::Bytes, reqwest::Error>>
                + Send
                + 'static,
        >,
    >,
}

impl ReqwestByteStream {
    pub fn new(resp: reqwest::Response) -> Self {
        Self {
            inner: Box::pin(resp.bytes_stream()),
        }
    }
}

impl futures::Stream for ReqwestByteStream {
    type Item = std::io::Result<bytes::Bytes>;
    fn poll_next(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
    ) -> Poll<Option<Self::Item>> {
        match self.inner.as_mut().poll_next(cx) {
            Poll::Ready(Some(Ok(b))) => Poll::Ready(Some(Ok(b))),
            Poll::Ready(Some(Err(e))) => Poll::Ready(Some(Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                e,
            )))),
            Poll::Ready(None) => Poll::Ready(None),
            Poll::Pending => Poll::Pending,
        }
    }
}

pub type BoxedReader = Pin<Box<dyn AsyncRead + Send + Unpin>>;

pub async fn stream_raw(
    backend: &dyn Backend,
    path: &str,
) -> Result<(BoxedReader, ShaHandle)> {
    let resp = backend.fetch(path).await?;
    let (rdr, h) = response_to_reader(resp);
    Ok((Box::pin(rdr), h))
}

/// Streaming download and on-the-fly gzip decompression

pub async fn stream_gz(
    backend: &dyn Backend,
    path: &str,
) -> Result<(BoxedReader, ShaHandle)> {
    let resp = backend.fetch(path).await?;
    let (rdr, h) = response_to_reader(resp);
    let buffered = BufReader::new(rdr);
    let decoder = GzipDecoder::new(buffered);
    Ok((Box::pin(decoder), h))
}

pub async fn with_retry<F, Fut, T>(label: &str, max_attempts: u32, mut f: F) -> Result<T>
where
    F: FnMut() -> Fut,
    Fut: std::future::Future<Output = Result<T>>,
{
    let mut delay = Duration::from_secs(1);
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
                    "  [retry] {} attempt {}/{} failed: {:#}; sleeping {:?}",
                    label, attempt, max_attempts, e, delay
                );
                tokio::time::sleep(delay).await;
                delay *= 2;
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
    backend: &dyn Backend,
    output_dir: &Path,
) -> Result<TaxdumpPaths> {
    let tax_dir = output_dir.join("taxonomy");
    tokio::fs::create_dir_all(&tax_dir).await?;
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
    backend: &dyn Backend,
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
        let version_status = cols[10].trim();
        let genome_rep = cols[13].trim();
        let asm_name = cols[15].trim();
        let ftp_path = cols[19].trim();

        if version_status != "latest" || genome_rep != "Full" {
            continue;
        }
        if domain == "human" && !asm_name.contains("GRCh") {
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
    eprintln!("  [{}] {} genomes after filter", domain, out.len());
    Ok(out)
}

pub async fn fetch_acc2taxid_filtered(
    backend: &dyn Backend,
    wanted: &std::collections::HashSet<String>,
) -> Result<FxHashMap<String, u32>> {
    let mut out: FxHashMap<String, u32> = FxHashMap::default();
    if wanted.is_empty() {
        return Ok(out);
    }

    for src in ["nucl_gb.accession2taxid.gz", "nucl_wgs.accession2taxid.gz"] {
        if out.len() >= wanted.len() {
            break;
        }
        let path = format!("/pub/taxonomy/accession2taxid/{}", src);

        let new_entries: Vec<(String, u32)> =
            with_retry(&format!("acc2taxid[{}]", src), 5, || async {
                let (rdr, _sha) = stream_gz(backend, &path).await?;
                stream_acc2taxid_into(rdr, wanted).await
            })
            .await?;

        for (k, v) in new_entries {
            out.entry(k).or_insert(v);
        }
    }
    Ok(out)
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