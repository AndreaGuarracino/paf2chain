extern crate clap;

use clap::{App, Arg}; //, SubCommand};

mod paf;

use crate::paf::{PafFile, paf_query_end, paf_target, paf_target_length, paf_target_begin, paf_target_end, paf_query, paf_query_length, paf_query_is_rev, paf_query_begin};

fn main() {
    let matches = App::new("paf2chain")
        .version("0.1.1")
        .author("Andrea Guarracino")
        .about("Generate a CHAIN format file from a PAF format file")
        .arg(
            Arg::with_name("INPUT")
                .required(true)
                .takes_value(true)
                .short("i")
                .long("input")
                .help("input PAF file"),
        )
        .get_matches();

    let filename = matches.value_of("INPUT").unwrap();
    let paf = PafFile::new(filename);

    let mut id = 0;
    let cigar_to_chain = |line: &str| {
        let query_rev = paf_query_is_rev(line);
        //let (x, y) = paf.global_start(line, query_rev);
        //let mut target_pos = x;
        //let mut query_pos = y;
        // find and walk the cigar string
        //println!("{}", line);

        let mut cigars = line
            .split('\t')
            .filter(|s| s.starts_with("cg:Z:"))
            .map(|s| s.strip_prefix("cg:Z:").unwrap())
            .collect::<Vec<&str>>();

        // Compute the fake CIGAR, in case the real one is missing in the line
        let query_len = paf_query_end(line) - paf_query_begin(line);
        let target_len = paf_target_end(line) - paf_target_begin(line);
        let mut fake_cigar = String::new();
        if cigars.is_empty() {
            fake_cigar.push_str(
                &format!(
                    "{}M",
                    if query_len < target_len {
                        query_len
                    } else {
                        target_len
                    }
                )
            );
            // If query and target have different lengths, we manage
            // the excess length to obtain a valid fake CIGAR string
            if query_len != target_len {
                if query_len < target_len {
                    fake_cigar.push_str(&format!("{}D", target_len - query_len));
                } else {
                    fake_cigar.push_str(&format!("{}I", query_len - target_len));
                }
            }
            cigars.push(fake_cigar.as_str());
        }

        let score = 255;
        let t_name = paf_target(line);
        let t_size = paf_target_length(line);
        let t_strand = "+";
        let t_start = paf_target_begin(line);
        let t_end = paf_target_end(line);
        let q_name = paf_query(line);
        let q_size = paf_query_length(line);
        let q_strand = if query_rev { "-" } else { "+" };
        let q_start = paf_query_begin(line);
        let q_end = paf_query_end(line);
        for mut cigar in cigars {
            //println!("{}", cigar);

            let mut first: usize = 0;

            // Trim INDELs at the beginning and at the end of the CIGAR string
            let mut first_not_indel_found = false;
            let mut trim_from: usize = 0;
            let mut query_from_delta: usize = 0;
            let mut target_from_delta: usize = 0;
            let mut trim_to: usize = 0;
            let mut query_to_delta: usize = 0;
            let mut target_to_delta: usize = 0;
            for (i, b) in cigar.bytes().enumerate() {
                let c = b as char;
                //println!("{} {}", i, b as char);
                match c {
                    'M' | '=' | 'X' => {
                        first = i + 1;

                        trim_to = i + 1;
                        query_to_delta = 0;
                        target_to_delta = 0;

                        first_not_indel_found = true;
                    }
                    'D' => {
                        let n = cigar[first..i].parse::<usize>().unwrap();
                        first = i + 1;

                        if !first_not_indel_found {
                            trim_from = i + 1;
                            target_from_delta += n;
                        }

                        target_to_delta += n;
                    }
                    'I' => {
                        let n = cigar[first..i].parse::<usize>().unwrap();
                        first = i + 1;

                        if !first_not_indel_found {
                            trim_from = i + 1;
                            query_from_delta += n;
                        }

                        query_to_delta += n;
                    }
                    _ => {}
                }
            }

            // Check if there are ungapped alignments to write
            //eprintln!("{} - {}", trim_from, trim_to);
            if trim_from < trim_to {
                cigar = &cigar[trim_from..trim_to];

                // Manage query start and end positions
                let chain_q_start;
                let chain_q_end;
                if query_rev {
                    chain_q_start = q_size - (q_end - query_from_delta);
                    chain_q_end = q_size - (q_start + query_to_delta);
                } else {
                    chain_q_start = q_start + query_from_delta;
                    chain_q_end = q_end - query_to_delta;
                }

                // Header line
                println!(
                    "chain\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    score,
                    t_name, t_size, t_strand,
                    t_start + target_from_delta,
                    t_end - target_to_delta,
                    q_name, q_size, q_strand,
                    chain_q_start, chain_q_end,
                    id
                );
                id += 1;

                let mut ungapped_alignment_len: usize = 0;
                let mut query_delta: usize = 0;
                let mut target_delta: usize = 0;

                first = 0;
                for (i, b) in cigar.bytes().enumerate() {
                    let c = b as char;
                    //println!("{} {}", i, b as char);
                    match c {
                        'M' | '=' | 'X' => {
                            let n = cigar[first..i].parse::<usize>().unwrap() as usize;
                            //query_pos += if query_rev { 0 - n } else { n };
                            //target_pos += n;
                            first = i + 1;

                            if ungapped_alignment_len > 0 && (query_delta > 0 || target_delta > 0) {
                                println!("{}\t{}\t{}", ungapped_alignment_len, target_delta, query_delta);

                                ungapped_alignment_len = 0;
                            }

                            target_delta = 0;
                            query_delta = 0;

                            ungapped_alignment_len += n;
                        }
                        'D' => {
                            let n = cigar[first..i].parse::<usize>().unwrap();
                            //target_pos += n;
                            first = i + 1;

                            target_delta += n;
                        }
                        'I' => {
                            let n = cigar[first..i].parse::<usize>().unwrap();
                            //query_pos += if query_rev { 0 - n } else { n };
                            first = i + 1;

                            query_delta += n;
                        }
                        _ => {}
                    }
                }

                println!("{}", ungapped_alignment_len);
            }
        }
        println!();
    };
    paf.for_each_line_in_file(cigar_to_chain);
}
