use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};
use crate::database::KmerDatabase;

pub struct InspectConfig {
    pub db_dir: String,
    pub output: String,
}

pub fn run_inspect(config: InspectConfig) -> Result<()> {
    let db_prefix = format!("{}/db", config.db_dir);

    // I am kind of lazy so I am just reusing KmerDatabase's load fn here
    // So this also loads heavy mphf too. I will be modifying this so load fn
    // is customizable which specific files to load. :)

    let db = KmerDatabase::load(&db_prefix, false)?;

    let n = db.index_to_taxid.len();
    let mut counts: Vec<usize> = vec![0; n];
    for &idx in db.taxid_indices() {
        let i = idx as usize;
        if i < counts.len() {
            counts[i] += 1;
        }
    }

    let total: usize = counts.iter().sum();

    let mut rows: Vec<(u32, String, usize)> = counts
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, &c)| {
            let taxid = db.true_taxid(i as u8);
            let name = db.get_taxid_name(taxid);
            (taxid, name, c)
        })
        .collect();

    rows.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)));

    let f = File::create(&config.output)
        .with_context(|| format!("Cannot create output file: {}", config.output))?;
    let mut writer = BufWriter::new(f);

    writeln!(writer, "TaxID\tName\tMinimizer_Count\tRatio")?;
    writeln!(writer, "Total\t\t{}\t{}", total as f64, 1.00)?;

    for (taxid, name, count) in &rows {
        let ratio = if total > 0 { *count as f64 / total as f64 } else { 0.0 };
        writeln!(writer, "{}\t{}\t{}\t{:.4}", taxid, name, count, ratio)?;
    }
    writer.flush()?;
    eprintln!("Done!");
    Ok(())
}