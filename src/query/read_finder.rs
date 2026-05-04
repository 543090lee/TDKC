use std::collections::HashMap;
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};

pub struct Sample {
    pub name: String,
    pub r1: PathBuf,
    pub r2: Option<PathBuf>,
}

// This will automatically detect reads, and run them sequentially
pub fn discover_samples(dir: &Path) -> Result<Vec<Sample>> {
    let mut sample_map: HashMap<String, Sample> = HashMap::new();

    for entry in std::fs::read_dir(dir).context("Failed to read input dir")? {
        let entry = entry?;
        let path = entry.path();

        if !path.is_file() {
            continue;
        }

        let file_name = path.file_name().unwrap_or_default().to_string_lossy().to_string();
        
        let is_fastq = file_name.ends_with(".fastq") 
            || file_name.ends_with(".fq") 
            || file_name.ends_with(".fastq.gz") 
            || file_name.ends_with(".fq.gz");

        if !is_fastq {
            continue;
        }

        let (sample_name, is_r1, is_r2) = if let Some(idx) = file_name.find("_R1") {
            (file_name[..idx].to_string(), true, false)
        } else if let Some(idx) = file_name.find("_R2") {
            (file_name[..idx].to_string(), false, true)
        } else {
            // now just assume this is single ended
            let name = file_name.split(".f").next().unwrap_or(&file_name).to_string();
            (name, true, false)
        };

        let entry = sample_map.entry(sample_name.clone()).or_insert(Sample {
            name: sample_name,
            r1: PathBuf::new(),
            r2: None,
        });

        if is_r1 {
            entry.r1 = path;
        } else if is_r2 {
            entry.r2 = Some(path);
        }
    }

    let mut valid_samples: Vec<Sample> = sample_map.into_values()
        .filter(|s| s.r1.exists())
        .collect();

    // sort so it processes from A01 through H12 cleanly
    valid_samples.sort_by(|a, b| a.name.cmp(&b.name));
    
    Ok(valid_samples)
}