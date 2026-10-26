use std::env;
use std::process;
use wfa::{biwfa_align, dp_align, read_fasta, wfa_align, PaddedSeq, WfaParams};

fn print_usage() {
    eprintln!(
        "Usage: wfa [options] -q <query.fa> -t <target.fa>\n\n\
         Options:\n\
           -q, --query <file>     Query sequence FASTA (required)\n\
           -t, --target <file>    Target sequence FASTA (required)\n\
           -m, --mode <mode>      Alignment mode: wfa (default), biwfa, dp\n\
           -x <int>               Mismatch score penalty (default: 4)\n\
           -o <int>               Gap open penalty (default: 6)\n\
           -e <int>               Gap extend penalty (default: 2)\n\
           -s, --show-alignment   Print visual alignment\n"
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }

    let mut query_path = String::new();
    let mut target_path = String::new();
    let mut mode = "wfa".to_string();
    let mut params = WfaParams::default();
    let mut show_alignment = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-q" | "--query" => {
                i += 1;
                if i < args.len() { query_path = args[i].clone(); }
            }
            "-t" | "--target" => {
                i += 1;
                if i < args.len() { target_path = args[i].clone(); }
            }
            "-m" | "--mode" => {
                i += 1;
                if i < args.len() { mode = args[i].clone(); }
            }
            "-x" => {
                i += 1;
                if i < args.len() { params.x = args[i].parse().unwrap_or(4); }
            }
            "-o" => {
                i += 1;
                if i < args.len() { params.o = args[i].parse().unwrap_or(6); }
            }
            "-e" => {
                i += 1;
                if i < args.len() { params.e = args[i].parse().unwrap_or(2); }
            }
            "-s" | "--show-alignment" => {
                show_alignment = true;
            }
            _ => {}
        }
        i += 1;
    }

    if query_path.is_empty() || target_path.is_empty() {
        eprintln!("Error: -q <query.fa> and -t <target.fa> are required.");
        print_usage();
        process::exit(1);
    }

    let queries = read_fasta(&query_path).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        process::exit(1);
    });
    let targets = read_fasta(&target_path).unwrap_or_else(|e| {
        eprintln!("Error reading target: {}", e);
        process::exit(1);
    });

    for (q_name, q_bytes) in &queries {
        for (t_name, t_bytes) in &targets {
            let q_pad = PaddedSeq::new(q_bytes, b'#');
            let t_pad = PaddedSeq::new(t_bytes, b'$');

            let res = match mode.as_str() {
                "dp" => dp_align(&q_pad, &t_pad, &params),
                "biwfa" => biwfa_align(&q_pad, &t_pad, &params),
                _ => wfa_align(&q_pad, &t_pad, &params),
            };

            println!(
                "Query: {}\tTarget: {}\tScore: {}\tCIGAR: {}",
                q_name, t_name, res.score, res.cigar
            );

            if show_alignment {
                println!("  Score: {} (X={}, O={}, E={})", res.score, params.x, params.o, params.e);
                println!("  CIGAR: {}", res.cigar);
            }
        }
    }
}
