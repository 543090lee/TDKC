use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::process::Stdio;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context as _, Result};
use rustc_hash::FxHashMap;
use tokio::fs::File as TokioFile;
use tokio::io::{AsyncReadExt, AsyncWriteExt, BufWriter as TokioBufWriter};
use tokio::process::Command;
use tokio::sync::{mpsc, Semaphore};
use tokio::task::JoinHandle;

use super::download::{stream_gz, stream_raw, stream_local_raw, stream_local_gz, with_retry, Backend};

// Gotta stream the file and split at record boundaries when chunk size exceeds the threshold,
// otherwise, dustmasking windows will break

#[derive(Clone, Debug)]
pub struct GenomeSource {
    pub fetch_path: String,
    pub gzipped: bool,
    pub taxid: u32,
    pub assembly_accession: String,
}

struct ByteChunk {
    bytes: Vec<u8>,
    taxid: u32,
    assembly_accession: String,
}

struct MaskedRecord {
    record_id: String,
    header: Vec<u8>,
    seq: Vec<u8>,
}

struct MaskedChunk {
    records: Vec<MaskedRecord>,
    taxid: u32,
    assembly_accession: String,
}

#[derive(Debug, Clone)]
pub struct ManifestRow {
    pub assembly_accession: String,
    pub source_label: String,
    pub taxid: u32,
    pub bytes_received: u64,
    pub sha256_hex: String,
}

pub struct SharedWriters {
    pub prelim_map_tx: mpsc::UnboundedSender<(String, u32)>,
    pub manifest_tx: mpsc::UnboundedSender<ManifestRow>,
}

pub struct SharedWriterHandles {
    pub prelim_map: JoinHandle<Result<()>>,
    pub manifest: JoinHandle<Result<()>>,
}

impl SharedWriters {
    pub fn new(genome_dir: &Path, domain_name: &str) -> Result<(Arc<Self>, SharedWriterHandles)> {
        let prelim_path = genome_dir.join(format!("_prelim.{}.txt", domain_name));
        let manifest_path = genome_dir.join(format!("_manifest.{}.tsv", domain_name));

        let (prelim_tx, mut prelim_rx) =
            mpsc::unbounded_channel::<(String, u32)>();
        let (manifest_tx, mut manifest_rx) = mpsc::unbounded_channel::<ManifestRow>();

        let prelim_handle: JoinHandle<Result<()>> = tokio::spawn(async move {
            let f = TokioFile::create(&prelim_path)
                .await
                .with_context(|| format!("create {}", prelim_path.display()))?;
            let mut w = TokioBufWriter::with_capacity(4 * 1024 * 1024, f);
            while let Some((acc, taxid)) = prelim_rx.recv().await {
                let line = format!("{}\t{}\n", acc, taxid);
                w.write_all(line.as_bytes()).await?;
            }
            w.flush().await?;
            Ok(())
        });

        let manifest_handle: JoinHandle<Result<()>> = tokio::spawn(async move {
            let f = TokioFile::create(&manifest_path)
                .await
                .with_context(|| format!("create {}", manifest_path.display()))?;
            let mut w = TokioBufWriter::with_capacity(1 * 1024 * 1024, f);
            while let Some(row) = manifest_rx.recv().await {
                let line = format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    row.assembly_accession,
                    row.source_label,
                    row.taxid,
                    row.bytes_received,
                    row.sha256_hex
                );
                w.write_all(line.as_bytes()).await?;
            }
            w.flush().await?;
            Ok(())
        });

        Ok((
            Arc::new(Self {
                prelim_map_tx: prelim_tx,
                manifest_tx,
            }),
            SharedWriterHandles {
                prelim_map: prelim_handle,
                manifest: manifest_handle,
            },
        ))
    }

    fn send_prelim(&self, record_id: String, taxid: u32) -> Result<()> {
        self.prelim_map_tx
            .send((record_id, taxid))
            .map_err(|_| anyhow!("prelim_map writer task died"))
    }

    fn send_manifest(&self, row: ManifestRow) -> Result<()> {
        self.manifest_tx
            .send(row)
            .map_err(|_| anyhow!("manifest writer task died"))
    }
}

pub struct PipelineConfig {
    pub domain_name: String,
    pub source_label: String,
    pub output_genome_dir: PathBuf,
    pub backend: Arc<dyn Backend>,
    pub concurrent_downloads: usize,
    pub max_in_flight_chunks: usize,
    pub batch_threshold_bytes: usize,
    pub relevant_taxids: Arc<HashSet<u32>>,
    pub shared: Arc<SharedWriters>,
    pub custom_map: Option<Arc<FxHashMap<String, u32>>>,
    pub no_mask: bool,
    /// When a custom_map is present, drop records whose id has no taxid mapping
    /// instead of emitting them with the fallback taxid (used for plasmid).
    pub drop_unmapped: bool,
}

#[derive(Debug, Default, Clone)]
pub struct PipelineStats {
    pub genomes_processed: u64,
    pub records_processed: u64,
    pub records_in_target: u64,
    pub bytes_received: u64,
    pub chunks_masked: u64,
}

pub async fn run_pipeline(
    cfg: PipelineConfig,
    sources: Vec<GenomeSource>,
) -> Result<PipelineStats> {
    let (chunks_tx, mut chunks_rx) =
        mpsc::channel::<ByteChunk>(cfg.max_in_flight_chunks.max(1));
    let (mask_tx, mut mask_rx) =
        mpsc::channel::<JoinHandle<Result<MaskedChunk>>>(cfg.max_in_flight_chunks.max(1));

    let backend = Arc::clone(&cfg.backend);
    let dl_concurrent = cfg.concurrent_downloads.max(1);
    let dl_threshold = cfg.batch_threshold_bytes;
    let dl_shared = Arc::clone(&cfg.shared);
    let dl_source_label = cfg.source_label.clone();
    let dl_domain = cfg.domain_name.clone();
    let downloader_handle: JoinHandle<Result<u64>> = tokio::spawn(async move {
        run_downloaders(
            sources,
            backend,
            dl_concurrent,
            dl_threshold,
            chunks_tx,
            dl_shared,
            dl_source_label,
            dl_domain,
        )
        .await
    });

    let in_flight_sem = Arc::new(Semaphore::new(cfg.max_in_flight_chunks.max(1)));
    let no_mask = cfg.no_mask;
    let masker_handle: JoinHandle<Result<u64>> = {
        let sem = Arc::clone(&in_flight_sem);
        tokio::spawn(async move {
            let mut chunk_count: u64 = 0;
            while let Some(chunk) = chunks_rx.recv().await {
                let permit = sem.clone().acquire_owned().await
                    .map_err(|_| anyhow!("semaphore closed"))?;
                let worker_handle: JoinHandle<Result<MaskedChunk>> =
                    tokio::spawn(async move {
                        let res = eof_worker(chunk, no_mask).await;
                        drop(permit);
                        res
                    });
                if mask_tx.send(worker_handle).await.is_err() {
                    return Err(anyhow!("sink dropped, masker exiting"));
                }
                chunk_count += 1;
            }
            Ok(chunk_count)
        })
    };

    let domain_path = cfg.output_genome_dir.join(format!("{}.fna", cfg.domain_name));
    let target_shard_path = cfg
        .output_genome_dir
        .join(format!("_target.{}.fna", cfg.domain_name));
    let sink_relevant = Arc::clone(&cfg.relevant_taxids);
    let sink_shared = Arc::clone(&cfg.shared);
    let sink_domain = cfg.domain_name.clone();

    let sink_custom_map = cfg.custom_map.clone();
    let sink_drop_unmapped = cfg.drop_unmapped;

    let sink_handle: JoinHandle<Result<PipelineStats>> = tokio::spawn(async move {
        run_sink(
            &mut mask_rx,
            domain_path,
            target_shard_path,
            sink_relevant,
            sink_shared,
            sink_domain,
            sink_custom_map,
            sink_drop_unmapped,
        )
        .await
    });

    let bytes_received = downloader_handle.await
        .map_err(|e| anyhow!("downloader task panicked: {e}"))??;
    let _chunks_dispatched = masker_handle.await
        .map_err(|e| anyhow!("masker task panicked: {e}"))??;
    let mut stats = sink_handle.await
        .map_err(|e| anyhow!("sink task panicked: {e}"))??;
    stats.bytes_received = bytes_received;

    Ok(stats)
}

async fn run_downloaders(
    sources: Vec<GenomeSource>,
    backend: Arc<dyn Backend>,
    concurrent: usize,
    threshold: usize,
    chunks_tx: mpsc::Sender<ByteChunk>,
    shared: Arc<SharedWriters>,
    source_label: String,
    domain_name: String,
) -> Result<u64> {
    let sem = Arc::new(Semaphore::new(concurrent));
    let total_bytes = Arc::new(std::sync::atomic::AtomicU64::new(0));

    let mut handles: Vec<JoinHandle<Result<()>>> = Vec::with_capacity(sources.len());
    for source in sources {
        let permit_sem = Arc::clone(&sem);
        let backend_c = Arc::clone(&backend);
        let chunks_tx_c = chunks_tx.clone();
        let shared_c = Arc::clone(&shared);
        let source_label_c = source_label.clone();
        let domain_c = domain_name.clone();
        let total_bytes_c = Arc::clone(&total_bytes);

        let h = tokio::spawn(async move {
            let _permit = permit_sem.acquire_owned().await
                .map_err(|_| anyhow!("download semaphore closed"))?;
            download_one(
                source,
                backend_c.as_ref(),
                threshold,
                &chunks_tx_c,
                &shared_c,
                &source_label_c,
                &domain_c,
                &total_bytes_c,
            )
            .await
        });
        handles.push(h);
    }

    drop(chunks_tx);

    for h in handles {
        h.await.map_err(|e| anyhow!("download task panicked: {e}"))??;
    }
    Ok(total_bytes.load(std::sync::atomic::Ordering::Relaxed))
}

async fn download_one(
    source: GenomeSource,
    backend: &dyn Backend,
    threshold: usize,
    chunks_tx: &mpsc::Sender<ByteChunk>,
    shared: &SharedWriters,
    source_label: &str,
    domain_name: &str,
    total_bytes: &std::sync::atomic::AtomicU64,
) -> Result<()> {
    let label = format!("download[{}/{}]", domain_name, source.assembly_accession);

    let (mut reader, sha) = with_retry(&label, 5, || async {
        // lets check if this path actually exists on your computer
        let is_local = std::path::Path::new(&source.fetch_path).exists();
        
        let pair = if is_local {
            if source.gzipped {
                stream_local_gz(&source.fetch_path).await?
            } else {
                stream_local_raw(&source.fetch_path).await?
            }
        } else {
            // It's not local, so fetch it from NCBI
            if source.gzipped {
                stream_gz(backend, &source.fetch_path).await?
            } else {
                stream_raw(backend, &source.fetch_path).await?
            }
        };
        Ok::<_, anyhow::Error>(pair)
    })
    .await?;

    let mut accum: Vec<u8> = Vec::with_capacity(threshold + 4 * 1024 * 1024);
    let mut read_buf = vec![0u8; 256 * 1024];
    let mut emitted_any = false;

    loop {
        let n = reader.read(&mut read_buf).await
            .with_context(|| format!("reading {}", source.fetch_path))?;
        if n == 0 {
            break;
        }
        let read_slice = &read_buf[..n];
        accum.extend_from_slice(read_slice);
        
        while accum.len() >= threshold {
            let split_at = find_record_boundary(&accum, !emitted_any);
            match split_at {
                Some(0) => {
                    break;
                }
                Some(idx) => {
                    let chunk_bytes: Vec<u8> = accum.drain(..idx).collect();
                    send_chunk(
                        chunks_tx,
                        chunk_bytes,
                        source.taxid,
                        source.assembly_accession.clone(),
                    )
                    .await?;
                    emitted_any = true;
                }
                None => {
                    bail!(
                        "no FASTA record boundary in {} after {} bytes (malformed input?)",
                        source.fetch_path,
                        accum.len()
                    );
                }
            }
        }
    }

    if !accum.is_empty() {
        if accum[0] != b'>' {
            bail!("trailing data in {} doesn't start with '>'",source.fetch_path);
        }
        send_chunk(
            chunks_tx,
            accum,
            source.taxid,
            source.assembly_accession.clone(),
        )
        .await?;
    }

    let (sha_hex, bytes) = sha.finalize().unwrap_or_else(|| ("".to_string(), 0));
    total_bytes.fetch_add(bytes, std::sync::atomic::Ordering::Relaxed);
    shared.send_manifest(ManifestRow {
        assembly_accession: source.assembly_accession.clone(),
        source_label: source_label.to_string(),
        taxid: source.taxid,
        bytes_received: bytes,
        sha256_hex: sha_hex,
    })?;

    Ok(())
}

#[inline]
async fn send_chunk(
    chunks_tx: &mpsc::Sender<ByteChunk>,
    bytes: Vec<u8>,
    taxid: u32,
    assembly_accession: String,
) -> Result<()> {
    chunks_tx
        .send(ByteChunk {
            bytes,
            taxid,
            assembly_accession,
        })
        .await
        .map_err(|_| anyhow!("masker stage closed channel"))
}

#[inline]
fn find_record_boundary(buf: &[u8], is_first_chunk_of_genome: bool) -> Option<usize> {
    if buf.is_empty() {
        return None;
    }
    
    if buf[0] != b'>' {
        return None;
    }
    if let Some(nl_pos) = rfind_nl_gt(buf) {
        return Some(nl_pos + 1); 
    }
    if is_first_chunk_of_genome {
        Some(0)
    } else {
        Some(0)
    }
}

#[inline]
fn rfind_nl_gt(buf: &[u8]) -> Option<usize> {
    if buf.len() < 2 {
        return None;
    }
    let mut i = buf.len() - 2;
    loop {
        if buf[i] == b'\n' && buf[i + 1] == b'>' {
            return Some(i);
        }
        if i == 0 {
            return None;
        }
        i -= 1;
    }
}

async fn eof_worker(chunk: ByteChunk, no_mask: bool) -> Result<MaskedChunk> {
    if no_mask {
        let taxid = chunk.taxid;
        let assembly = chunk.assembly_accession;
        let chunk_bytes = chunk.bytes;
        let parsed = tokio::task::spawn_blocking(move || -> Result<Vec<MaskedRecord>> {
            parse_masked_fasta(&chunk_bytes)
        })
        .await
        .context("parse_masked_fasta task")??;
        return Ok(MaskedChunk {
            records: parsed,
            taxid,
            assembly_accession: assembly,
        });
    }

    let mut child = Command::new("dustmasker")
        .args(["-in", "-", "-outfmt", "fasta"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .context("spawning dustmasker (is BLAST+ installed?)")?;

    let mut stdin = child.stdin.take().ok_or_else(|| anyhow!("dustmasker stdin missing"))?;
    let mut stdout = child.stdout.take().ok_or_else(|| anyhow!("dustmasker stdout missing"))?;
    let mut stderr = child.stderr.take().ok_or_else(|| anyhow!("dustmasker stderr missing"))?;

    let chunk_bytes = chunk.bytes;
    let writer_fut = async move {
        stdin
            .write_all(&chunk_bytes)
            .await
            .context("writing to dustmasker stdin")?;
        stdin.shutdown().await.ok(); 
        drop(stdin);
        Ok::<usize, anyhow::Error>(chunk_bytes.len())
    };
    let reader_fut = async {
        let mut out = Vec::with_capacity(64 * 1024 * 1024);
        tokio::io::copy(&mut stdout, &mut out)
            .await
            .context("reading dustmasker stdout")?;
        Ok::<Vec<u8>, anyhow::Error>(out)
    };
    let stderr_fut = async {
        let mut err_out = Vec::new();
        let _ = stderr.read_to_end(&mut err_out).await;
        err_out
    };
    let (write_res, read_res, err_bytes) =
        tokio::join!(writer_fut, reader_fut, stderr_fut);
    let _written = write_res?;
    let masked_bytes = read_res?;

    let status = child.wait().await.context("waiting on dustmasker")?;
    if !status.success() {
        let err_msg = String::from_utf8_lossy(&err_bytes);
        bail!(
            "dustmasker exited with {}; stderr:\n{}",
            status,
            err_msg.trim()
        );
    }

    let taxid = chunk.taxid;
    let assembly = chunk.assembly_accession;
    let parsed = tokio::task::spawn_blocking(move || -> Result<Vec<MaskedRecord>> {
        parse_masked_fasta(&masked_bytes)
    })
    .await
    .context("parse_masked_fasta task")??;

    Ok(MaskedChunk {
        records: parsed,
        taxid,
        assembly_accession: assembly,
    })
}

fn parse_masked_fasta(bytes: &[u8]) -> Result<Vec<MaskedRecord>> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'>' {
            i += 1;
            continue;
        }
        i += 1;
        let header_start = i;
        while i < bytes.len() && bytes[i] != b'\n' {
            i += 1;
        }
        let header_end = if i > header_start && bytes[i - 1] == b'\r' {
            i - 1
        } else {
            i
        };
        let header = bytes[header_start..header_end].to_vec();
        if i < bytes.len() {
            i += 1;
        }

        let record_id = match std::str::from_utf8(&header) {
            Ok(s) => s.split_whitespace().next().unwrap_or("").to_string(),
            Err(_) => continue, 
        };
        if record_id.is_empty() {
        }

        let mut seq = Vec::with_capacity(64);
        while i < bytes.len() && bytes[i] != b'>' {
            let b = bytes[i];
            if b != b'\n' && b != b'\r' {
                seq.push(b);
            }
            i += 1;
        }

        if !record_id.is_empty() {
            out.push(MaskedRecord {
                record_id,
                header,
                seq,
            });
        }
    }
    Ok(out)
}

async fn run_sink(
    mask_rx: &mut mpsc::Receiver<JoinHandle<Result<MaskedChunk>>>,
    domain_path: PathBuf,
    target_shard_path: PathBuf,
    relevant: Arc<HashSet<u32>>,
    shared: Arc<SharedWriters>,
    domain_name: String,
    custom_map: Option<Arc<FxHashMap<String, u32>>>,
    drop_unmapped: bool,
) -> Result<PipelineStats> {
    let domain_file = TokioFile::create(&domain_path)
        .await
        .with_context(|| format!("create {}", domain_path.display()))?;
    let mut domain_w = TokioBufWriter::with_capacity(8 * 1024 * 1024, domain_file);

    let target_file = TokioFile::create(&target_shard_path)
        .await
        .with_context(|| format!("create {}", target_shard_path.display()))?;
    let mut target_w = TokioBufWriter::with_capacity(8 * 1024 * 1024, target_file);

    let mut seen_assemblies: FxHashMap<String, ()> = FxHashMap::default();

    let mut stats = PipelineStats::default();

    while let Some(handle) = mask_rx.recv().await {
        let masked = handle.await
            .map_err(|e| anyhow!("[{}] masker task panicked: {e}", domain_name))??;
        stats.chunks_masked += 1;
        seen_assemblies.entry(masked.assembly_accession.clone()).or_insert(());

        for mut rec in masked.records {
            lowercase_to_x(&mut rec.seq);

            let mut actual_taxid = masked.taxid;

            if let Some(cmap) = &custom_map {
                match cmap.get(&rec.record_id) {
                    Some(&tid) => actual_taxid = tid,
                    None if drop_unmapped => continue,
                    None => {}
                }
            }

            let is_target = relevant.contains(&actual_taxid);
            
            domain_w.write_all(b">").await?;
            domain_w.write_all(&rec.header).await?;
            domain_w.write_all(b"\n").await?;
            write_wrapped_fasta_seq(&mut domain_w, &rec.seq, 80).await?;

            if is_target {
                target_w.write_all(b">").await?;
                target_w.write_all(&rec.header).await?;
                target_w.write_all(b"\n").await?;
                write_wrapped_fasta_seq(&mut target_w, &rec.seq, 80).await?;
                stats.records_in_target += 1;
            }

            shared.send_prelim(rec.record_id, actual_taxid)?;
            stats.records_processed += 1;
        }
    }

    domain_w.flush().await?;
    target_w.flush().await?;
    stats.genomes_processed = seen_assemblies.len() as u64;
    Ok(stats)
}

#[inline]
fn lowercase_to_x(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if (b'a'..=b'z').contains(b) {
            *b = b'x';
        }
    }
}

async fn write_wrapped_fasta_seq<W: AsyncWriteExt + Unpin>(
    w: &mut W,
    seq: &[u8],
    width: usize,
) -> Result<()> {
    let mut i = 0;
    while i < seq.len() {
        let end = (i + width).min(seq.len());
        w.write_all(&seq[i..end]).await?;
        w.write_all(b"\n").await?;
        i = end;
    }
    Ok(())
}

fn collect_shards(genome_dir: &Path, prefix: &str, suffix: &str) -> Result<Vec<PathBuf>> {
    let mut shards = Vec::new();
    let rd = match std::fs::read_dir(genome_dir) {
        Ok(rd) => rd,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(shards),
        Err(e) => return Err(e).with_context(|| format!("read_dir {}", genome_dir.display())),
    };
    for entry in rd {
        let entry = entry?;
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if name.starts_with(prefix) && name.ends_with(suffix) {
            shards.push(entry.path());
        }
    }
    shards.sort();
    Ok(shards)
}

async fn concat_shards(
    shards: &[PathBuf],
    out_path: &Path,
    header: Option<&[u8]>,
) -> Result<()> {
    let out_file = TokioFile::create(out_path)
        .await
        .with_context(|| format!("create {}", out_path.display()))?;
    let mut w = TokioBufWriter::with_capacity(8 * 1024 * 1024, out_file);
    if let Some(h) = header {
        w.write_all(h).await?;
    }
    for shard in shards {
        let mut f = TokioFile::open(shard)
            .await
            .with_context(|| format!("open shard {}", shard.display()))?;
        tokio::io::copy(&mut f, &mut w)
            .await
            .with_context(|| format!("copying shard {}", shard.display()))?;
    }
    w.flush().await?;
    Ok(())
}

pub async fn finalize_outputs(out_dir: &Path, genome_dir: &Path) -> Result<()> {
    let target_shards = collect_shards(genome_dir, "_target.", ".fna")?;
    let prelim_shards = collect_shards(genome_dir, "_prelim.", ".txt")?;
    let manifest_shards = collect_shards(genome_dir, "_manifest.", ".tsv")?;

    concat_shards(&target_shards, &out_dir.join("target.fasta"), None).await?;
    concat_shards(
        &prelim_shards,
        &out_dir.join("prelim_map.txt"),
        Some(b"accession\ttaxid\n"),
    )
    .await?;
    concat_shards(
        &manifest_shards,
        &genome_dir.join("manifest.tsv"),
        Some(b"assembly_accession\tsource\ttaxid\tbytes_received\tsha256\n"),
    )
    .await?;

    Ok(())
}

pub async fn finalize_shared(
    shared: Arc<SharedWriters>,
    handles: SharedWriterHandles,
) -> Result<()> {

    drop(shared);
    handles.prelim_map.await
        .map_err(|e| anyhow!("prelim_map writer panicked: {e}"))??;
    handles.manifest.await
        .map_err(|e| anyhow!("manifest writer panicked: {e}"))??;
    Ok(())
}