
use anyhow::{anyhow, bail, Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use flate2::{Decompress, FlushDecompress, Status};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, Read};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread::JoinHandle;

pub fn read_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::FileExt;
        file.read_exact_at(buf, offset)
    }
    #[cfg(windows)]
    {
        use std::os::windows::fs::FileExt;
        let mut done = 0usize;
        while done < buf.len() {
            let n = file.seek_read(&mut buf[done..], offset + done as u64)?;
            if n == 0 {
                return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "short read"));
            }
            done += n;
        }
        Ok(())
    }
}

#[inline]
fn plausible_gzip_header(b: &[u8]) -> bool {
    if b.len() < 10 || b[0] != 0x1f || b[1] != 0x8b || b[2] != 0x08 {
        return false;
    }
    if b[3] & 0xE0 != 0 {
        return false;
    }
    let xfl = b[8];
    if xfl != 0 && xfl != 2 && xfl != 4 {
        return false;
    }
    let os = b[9];
    os <= 13 || os == 255
}

fn parse_gzip_header(b: &[u8]) -> Result<usize> {
    if b.len() < 10 || b[0] != 0x1f || b[1] != 0x8b {
        bail!("not a gzip member header");
    }
    if b[2] != 8 {
        bail!("unsupported gzip compression method {}", b[2]);
    }
    let flg = b[3];
    let mut p = 10usize;
    if flg & 0x04 != 0 {
        if p + 2 > b.len() {
            bail!("truncated gzip FEXTRA");
        }
        let xlen = u16::from_le_bytes([b[p], b[p + 1]]) as usize;
        p += 2 + xlen;
    }
    if flg & 0x08 != 0 {
        let z = memchr::memchr(0, &b[p.min(b.len())..]).ok_or_else(|| anyhow!("truncated gzip FNAME"))?;
        p += z + 1;
    }
    if flg & 0x10 != 0 {
        let z = memchr::memchr(0, &b[p.min(b.len())..]).ok_or_else(|| anyhow!("truncated gzip FCOMMENT"))?;
        p += z + 1;
    }
    if flg & 0x02 != 0 {
        p += 2;
    }
    if p > b.len() {
        bail!("truncated gzip header");
    }
    Ok(p)
}

fn scan_gzip_candidates(file: &File, file_len: usize) -> Result<Vec<usize>> {
    const BLOCK: usize = 32 << 20;
    const OVERLAP: usize = 16;
    let mut cands = Vec::new();
    let mut buf = vec![0u8; BLOCK + OVERLAP];
    let mut off = 0usize;
    while off < file_len {
        let n = (BLOCK + OVERLAP).min(file_len - off);
        read_at(file, &mut buf[..n], off as u64)?;
        let limit = if off + n >= file_len { n } else { n - OVERLAP };
        for p in memchr::memmem::find_iter(&buf[..n], &[0x1f, 0x8b, 0x08]) {
            if p < limit && plausible_gzip_header(&buf[p..n]) {
                cands.push(off + p);
            }
        }
        off += BLOCK;
    }
    Ok(cands)
}

struct InputWindow<'a> {
    file: &'a File,
    file_len: usize,
    base: usize,
    buf: Vec<u8>,
}

impl<'a> InputWindow<'a> {
    const REFILL: usize = 4 << 20;

    fn new(file: &'a File, file_len: usize, start: usize, initial: usize, mut buf: Vec<u8>) -> Result<Self> {
        let n = initial.min(file_len - start);
        buf.clear();
        buf.resize(n, 0);
        read_at(file, &mut buf, start as u64)?;
        Ok(Self { file, file_len, base: start, buf })
    }

    fn into_buf(self) -> Vec<u8> {
        self.buf
    }

    fn ensure(&mut self, pos: usize, min_len: usize) -> Result<()> {
        let want_end = (pos + min_len).min(self.file_len);
        let have_end = self.base + self.buf.len();
        if want_end <= have_end {
            return Ok(());
        }
        if pos - self.base > (8 << 20) {
            self.buf.drain(..pos - self.base);
            self.base = pos;
        }
        let have_end = self.base + self.buf.len();
        let new_end = (have_end + Self::REFILL).max(want_end).min(self.file_len);
        let old_len = self.buf.len();
        self.buf.resize(new_end - self.base, 0);
        read_at(self.file, &mut self.buf[old_len..], have_end as u64)?;
        Ok(())
    }

    #[inline]
    fn slice_from(&self, pos: usize) -> &[u8] {
        &self.buf[pos - self.base..]
    }

    #[inline]
    fn buffered_end(&self) -> usize {
        self.base + self.buf.len()
    }
}

struct Chunk {
    task_idx: u64,
    seq: u64,
    last: bool,
    actual_end: usize,
    data: Vec<u8>,
    err: Option<String>,
}

struct Input {
    file: File,
    file_len: usize,
    task_starts: Vec<usize>,
    tasks: Vec<(usize, usize)>,
}

const STREAM_CHUNK: usize = 16 << 20;

fn take_buf(pool: &Receiver<Vec<u8>>, cap: usize) -> Vec<u8> {
    match pool.try_recv() {
        Ok(mut v) => {
            v.clear();
            v
        }
        Err(_) => Vec::with_capacity(cap),
    }
}

fn decode_task(input: &Input, task_idx: usize, tx: &Sender<Chunk>, pool: &Receiver<Vec<u8>>, mut data: Vec<u8>, win_buf: &mut Vec<u8>) -> Result<()> {
    let (start, end_hint) = input.tasks[task_idx];
    let file_len = input.file_len;
    let mut win = InputWindow::new(&input.file, file_len, start, (end_hint - start) + 1024, std::mem::take(win_buf))?;
    let mut pos = start;
    let mut seq: u64 = 0;
    data.clear();
    data.reserve(6 * (end_hint - start).max(256 << 10));
    let mut inflater = Decompress::new(false);
    let mut crc = crc32fast::Hasher::new();
    let overshoot_cap = end_hint + 256 * 1024 * 1024;

    loop {
        if pos >= file_len {
            break;
        }
        win.ensure(pos, 64 * 1024)?;
        let hdr = parse_gzip_header(win.slice_from(pos)).with_context(|| format!("gzip header at offset {}", pos))?;
        let body_start = pos + hdr;
        inflater.reset(false);
        crc.reset();
        let member_out_start = data.len();
        let mut consumed = 0usize;
        loop {
            if data.capacity() - data.len() < 64 * 1024 {
                data.reserve_exact((data.len() / 4).max(1 << 20));
            }
            let in_pos = body_start + consumed;
            win.ensure(in_pos, 64 * 1024)?;
            let before_out = data.len();
            let before_in = inflater.total_in();
            let status = inflater
                .decompress_vec(win.slice_from(in_pos), &mut data, FlushDecompress::None)
                .with_context(|| format!("inflate error in member at offset {}", pos))?;
            consumed += (inflater.total_in() - before_in) as usize;
            match status {
                Status::StreamEnd => break,
                Status::Ok => {
                    if body_start + consumed >= file_len && data.len() == before_out {
                        bail!("truncated gzip member at offset {}", pos);
                    }
                }
                Status::BufError => {
                    if body_start + consumed >= file_len {
                        bail!("truncated gzip member at offset {}", pos);
                    }
                    if data.len() == before_out && body_start + consumed < win.buffered_end() {
                        data.reserve_exact(1 << 20);
                    }
                }
            }
        }
        let trailer = body_start + consumed;
        if trailer + 8 > file_len {
            bail!("truncated gzip trailer at offset {}", trailer);
        }
        win.ensure(trailer, 8)?;
        let member_out = &data[member_out_start..];
        crc.update(member_out);
        let tb = win.slice_from(trailer);
        let crc_expect = u32::from_le_bytes(tb[0..4].try_into().unwrap());
        let isize_expect = u32::from_le_bytes(tb[4..8].try_into().unwrap());
        if isize_expect != (member_out.len() as u32) {
            bail!("gzip ISIZE mismatch in member at offset {}", pos);
        }
        if crc.clone().finalize() != crc_expect {
            bail!("gzip CRC mismatch in member at offset {}", pos);
        }
        pos = trailer + 8;

        let done = pos >= file_len || (pos >= end_hint && input.task_starts.binary_search(&pos).is_ok());
        if done {
            break;
        }
        if pos > overshoot_cap {
            bail!("gzip member chain starting at {} never reached a task boundary", start);
        }
        if data.len() >= STREAM_CHUNK {
            let chunk = std::mem::replace(&mut data, take_buf(pool, STREAM_CHUNK + (1 << 20)));
            if tx.send(Chunk { task_idx: task_idx as u64, seq, last: false, actual_end: pos, data: chunk, err: None }).is_err() {
                *win_buf = win.into_buf();
                return Ok(());
            }
            seq += 1;
        }
    }
    let _ = tx.send(Chunk { task_idx: task_idx as u64, seq, last: true, actual_end: pos, data, err: None });
    *win_buf = win.into_buf();
    Ok(())
}

pub struct ParallelGzReader {
    rx: Receiver<Chunk>,
    pool_tx: Sender<Vec<u8>>,
    input: Arc<Input>,
    pending: BTreeMap<(u64, u64), Chunk>,
    errors: BTreeMap<u64, String>,
    cur: Vec<u8>,
    cur_pos: usize,
    cur_task: u64,
    cur_seq: u64,
    eof: bool,
    failed: Arc<AtomicBool>,
    threads: Vec<JoinHandle<()>>,
}

impl ParallelGzReader {
    pub fn open(path: &str, threads: usize) -> Result<Option<Self>> {
        let file = File::open(path).with_context(|| format!("Cannot open reads file {}", path))?;
        let file_len = file.metadata()?.len() as usize;
        if file_len < 18 {
            return Ok(None);
        }
        let mut head = [0u8; 10];
        read_at(&file, &mut head, 0)?;
        if head[0] != 0x1f || head[1] != 0x8b || head[2] != 8 {
            return Ok(None);
        }
        let threads = threads.max(1);
        let mut cands = scan_gzip_candidates(&file, file_len)?;
        if cands.first() != Some(&0) {
            cands.insert(0, 0);
        }
        let span = (file_len / (threads * 8)).clamp(512 << 10, 2 << 20);
        let mut starts: Vec<usize> = vec![0];
        let mut last_start = 0usize;
        for &c in &cands[1..] {
            if c - last_start >= span {
                starts.push(c);
                last_start = c;
            }
        }
        let mut tasks = Vec::with_capacity(starts.len());
        for i in 0..starts.len() {
            let end = if i + 1 < starts.len() { starts[i + 1] } else { file_len };
            tasks.push((starts[i], end));
        }
        let input = Arc::new(Input { file, file_len, task_starts: starts, tasks });
        let (tx, rx) = bounded::<Chunk>(3);
        let pool_cap = threads + 4;
        let (pool_tx, pool_rx) = bounded::<Vec<u8>>(pool_cap);
        for _ in 0..pool_cap {
            let _ = pool_tx.try_send(Vec::new());
        }
        let cursor = Arc::new(AtomicUsize::new(0));
        let failed = Arc::new(AtomicBool::new(false));
        let mut handles = Vec::new();
        for _ in 0..threads {
            let inp = input.clone();
            let cursor = cursor.clone();
            let tx = tx.clone();
            let pool_rx = pool_rx.clone();
            let failed = failed.clone();
            handles.push(std::thread::spawn(move || {
                let mut win_buf: Vec<u8> = Vec::new();
                loop {
                if failed.load(Ordering::Relaxed) {
                    break;
                }
                let i = cursor.fetch_add(1, Ordering::Relaxed);
                if i >= inp.tasks.len() {
                    break;
                }
                let data = match pool_rx.recv() {
                    Ok(v) => v,
                    Err(_) => break,
                };
                if let Err(e) = decode_task(&inp, i, &tx, &pool_rx, data, &mut win_buf) {
                    let _ = tx.send(Chunk {
                        task_idx: i as u64,
                        seq: u64::MAX,
                        last: true,
                        actual_end: 0,
                        data: Vec::new(),
                        err: Some(format!("{:#}", e)),
                    });
                }
                }
            }));
        }
        drop(tx);
        drop(pool_rx);
        Ok(Some(Self {
            rx,
            pool_tx,
            input,
            pending: BTreeMap::new(),
            errors: BTreeMap::new(),
            cur: Vec::new(),
            cur_pos: 0,
            cur_task: 0,
            cur_seq: 0,
            eof: false,
            failed,
            threads: handles,
        }))
    }

    fn advance(&mut self) -> io::Result<bool> {
        loop {
            if let Some(c) = self.pending.remove(&(self.cur_task, self.cur_seq)) {
                if c.last {
                    let next_pos = c.actual_end;
                    if next_pos >= self.input.file_len {
                        self.eof = true;
                    } else {
                        self.cur_task = match self.input.task_starts.binary_search(&next_pos) {
                            Ok(i) => i as u64,
                            Err(_) => {
                                return Err(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("gzip member chain ended at {} which is not a task start", next_pos),
                                ))
                            }
                        };
                        self.cur_seq = 0;
                        let stale: Vec<(u64, u64)> = self.pending.range(..(self.cur_task, 0)).map(|(k, _)| *k).collect();
                        for k in stale {
                            if let Some(c) = self.pending.remove(&k) {
                                let _ = self.pool_tx.try_send(c.data);
                            }
                        }
                        let stale_err: Vec<u64> = self.errors.range(..self.cur_task).map(|(k, _)| *k).collect();
                        for k in stale_err {
                            self.errors.remove(&k);
                        }
                    }
                } else {
                    self.cur_seq += 1;
                }
                let done = std::mem::replace(&mut self.cur, c.data);
                if done.capacity() > 0 {
                    let _ = self.pool_tx.try_send(done);
                }
                self.cur_pos = 0;
                return Ok(true);
            }
            if self.eof {
                return Ok(false);
            }
            if let Some(e) = self.errors.get(&self.cur_task) {
                let msg = e.clone();
                self.failed.store(true, Ordering::SeqCst);
                return Err(io::Error::new(io::ErrorKind::InvalidData, msg));
            }
            match self.rx.recv() {
                Ok(c) => {
                    if c.task_idx < self.cur_task {
                        let _ = self.pool_tx.try_send(c.data);
                        continue;
                    }
                    if let Some(e) = c.err {
                        if c.task_idx == self.cur_task {
                            self.failed.store(true, Ordering::SeqCst);
                            return Err(io::Error::new(io::ErrorKind::InvalidData, e));
                        }
                        self.errors.insert(c.task_idx, e);
                        continue;
                    }
                    self.pending.insert((c.task_idx, c.seq), c);
                }
                Err(_) => {
                    return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "gzip decode pool ended before the file was fully read"));
                }
            }
        }
    }
}

impl Read for ParallelGzReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        while self.cur_pos >= self.cur.len() {
            if !self.advance()? {
                return Ok(0);
            }
        }
        let n = (self.cur.len() - self.cur_pos).min(buf.len());
        buf[..n].copy_from_slice(&self.cur[self.cur_pos..self.cur_pos + n]);
        self.cur_pos += n;
        Ok(n)
    }
}

impl Drop for ParallelGzReader {
    fn drop(&mut self) {
        self.failed.store(true, Ordering::SeqCst);
        let rx = std::mem::replace(&mut self.rx, bounded(1).1);
        drop(rx);
        let pool_tx = std::mem::replace(&mut self.pool_tx, bounded(1).0);
        drop(pool_tx);
        for t in self.threads.drain(..) {
            let _ = t.join();
        }
    }
}
