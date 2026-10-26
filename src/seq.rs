use std::fs::File;
use std::io::{BufRead, BufReader};

pub struct PaddedSeq {
    data: Vec<u8>,
    orig_len: usize,
}

impl PaddedSeq {
    pub fn new(raw: &[u8], sentinel: u8) -> Self {
        let mut data = Vec::with_capacity(raw.len() + 16);
        data.extend_from_slice(raw);
        // Pad with 16 distinct sentinel bytes so 64-bit word loads never read OOB
        for _ in 0..16 {
            data.push(sentinel);
        }
        Self {
            orig_len: raw.len(),
            data,
        }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.orig_len
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.orig_len == 0
    }

    #[inline(always)]
    pub fn as_slice(&self) -> &[u8] {
        &self.data[..self.orig_len]
    }

    #[inline(always)]
    pub fn padded_slice(&self) -> &[u8] {
        &self.data
    }

    #[inline(always)]
    pub fn get(&self, idx: usize) -> u8 {
        self.data[idx]
    }
}

pub fn read_fasta(path: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    let file = File::open(path).map_err(|e| format!("Failed to open {}: {}", path, e))?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut cur_header = String::new();
    let mut cur_seq = Vec::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(|e| e.to_string())?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with('>') {
            if !cur_header.is_empty() {
                records.push((cur_header, cur_seq));
                cur_seq = Vec::new();
            }
            let name = line[1..].split_whitespace().next().unwrap_or("seq");
            cur_header = name.to_string();
        } else {
            for b in line.bytes() {
                if b > b' ' {
                    cur_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }

    if !cur_header.is_empty() {
        records.push((cur_header, cur_seq));
    }

    Ok(records)
}
