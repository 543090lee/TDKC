use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result};

struct IndexEntry {
    length: usize,
    offset: u64,
    line_bases: usize,
    line_width: usize,
}

/// Indexed FASTA reader for random access by sequence name.
pub struct FastaIndex {
    path: String,
    index: HashMap<String, IndexEntry>,
}

impl FastaIndex {
    pub fn new(fasta_path: &str) -> Result<Self> {
        let fai_path = format!("{}.fai", fasta_path);

        let index = if Path::new(&fai_path).exists() {
            Self::load_fai(&fai_path)?
        } else {
            eprintln!("Creating FASTA index...");
            let idx = Self::build_fai(fasta_path)?;
            Self::save_fai(&fai_path, &idx)?;
            idx
        };

        eprintln!("Indexed {} sequences", index.len());

        Ok(Self {
            path: fasta_path.to_string(),
            index,
        })
    }

    /// Fetch a sequence by name. Returns None if not found.
    pub fn get_sequence(&self, name: &str) -> Option<String> {
        let entry = self.index.get(name)?;

        let mut file = File::open(&self.path).ok()?;
        file.seek(SeekFrom::Start(entry.offset)).ok()?;

        let mut sequence = String::with_capacity(entry.length);
        let mut reader = BufReader::new(file);
        let mut bases_read = 0;

        loop {
            let mut line = String::new();
            match reader.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {}
                Err(_) => break,
            }

            let line = line.trim_end();
            if line.is_empty() || line.starts_with('>') {
                break;
            }

            for c in line.chars() {
                if c != '\n' && c != '\r' {
                    sequence.push(c);
                    bases_read += 1;
                    if bases_read >= entry.length {
                        return Some(sequence);
                    }
                }
            }
        }

        Some(sequence)
    }

    fn load_fai(path: &str) -> Result<HashMap<String, IndexEntry>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut index = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let name = parts[0].to_string();
                let entry = IndexEntry {
                    length: parts[1].parse()?,
                    offset: parts[2].parse()?,
                    line_bases: parts[3].parse()?,
                    line_width: parts[4].parse()?,
                };
                index.insert(name, entry);
            }
        }

        Ok(index)
    }

    fn build_fai(fasta_path: &str) -> Result<HashMap<String, IndexEntry>> {
        let file = File::open(fasta_path).context("Cannot open FASTA")?;
        let reader = BufReader::new(file);
        let mut index = HashMap::new();

        let mut current_name: Option<String> = None;
        let mut seq_length: usize = 0;
        let mut seq_start_offset: u64 = 0;
        let mut line_bases: usize = 0;
        let mut line_width: usize = 0;
        let mut offset: u64 = 0;

        for line in reader.lines() {
            let line = line?;
            let line_len = line.len() as u64 + 1; // +1 for newline

            if line.starts_with('>') {
                // Save previous entry
                if let Some(name) = current_name.take() {
                    index.insert(
                        name,
                        IndexEntry {
                            length: seq_length,
                            offset: seq_start_offset,
                            line_bases,
                            line_width,
                        },
                    );
                }

                // Parse new name (up to first space)
                let header = &line[1..];
                let name = header
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                current_name = Some(name);
                seq_length = 0;
                seq_start_offset = offset + line_len;
                line_bases = 0;
                line_width = 0;
            } else if !line.is_empty() {
                let bases = line.chars().filter(|c| *c != '\r').count();
                seq_length += bases;
                if line_bases == 0 {
                    line_bases = bases;
                    line_width = line.len() + 1; // include newline
                }
            }

            offset += line_len;
        }

        // Save last entry
        if let Some(name) = current_name {
            index.insert(
                name,
                IndexEntry {
                    length: seq_length,
                    offset: seq_start_offset,
                    line_bases,
                    line_width,
                },
            );
        }

        Ok(index)
    }

    fn save_fai(path: &str, index: &HashMap<String, IndexEntry>) -> Result<()> {
        use std::io::Write;
        let mut file = File::create(path)?;
        for (name, entry) in index {
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}",
                name, entry.length, entry.offset, entry.line_bases, entry.line_width
            )?;
        }
        Ok(())
    }
}