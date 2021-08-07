use boomphf::Mphf;
use itertools::Itertools;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
struct AlignedSeq {
    // name of the given sequence
    name: String,
    // its length
    length: usize,
    // its rank among other seqs in the query or target set
    rank: usize,
    // its start offset in the global all-to-all alignment matrix
    offset: usize,
}

impl AlignedSeq {
    fn new() -> Self {
        AlignedSeq {
            name: String::new(),
            length: 0,
            rank: 0,
            offset: 0,
        }
    }
}

pub struct PafFile {
    // our input file
    filename: String,
    // each target name in order of first appearance
    targets: Vec<AlignedSeq>,
    // each query name in order of first appearance
    queries: Vec<AlignedSeq>,
    // maps from sequence name to internal id
    target_mphf: Mphf<String>,
    // maps from sequence name to internal id
    query_mphf: Mphf<String>,
    // target axis length
    target_length: f64,
    // query axis length
    query_length: f64,
}

fn paf_query(line: &str) -> String {
    line.split('\t').next().unwrap().into()
}

fn paf_query_length(line: &str) -> usize {
    line.split('\t').nth(1).unwrap().parse::<usize>().unwrap()
}

fn paf_query_begin(line: &str) -> usize {
    line.split('\t').nth(2).unwrap().parse::<usize>().unwrap()
}

fn paf_query_end(line: &str) -> usize {
    line.split('\t').nth(3).unwrap().parse::<usize>().unwrap()
}

fn paf_query_is_rev(line: &str) -> bool {
    line.split('\t').nth(4).unwrap() == "-"
}

fn paf_target(line: &str) -> String {
    line.split('\t').nth(5).unwrap().into()
}

fn paf_target_length(line: &str) -> usize {
    line.split('\t').nth(6).unwrap().parse::<usize>().unwrap()
}

fn paf_target_begin(line: &str) -> usize {
    line.split('\t').nth(7).unwrap().parse::<usize>().unwrap()
}

fn paf_target_end(line: &str) -> usize {
    line.split('\t').nth(8).unwrap().parse::<usize>().unwrap()
}

fn for_each_line_in_file(paf_filename: &str, mut callback: impl FnMut(&str)) {
    let file = File::open(paf_filename).unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        callback(&line.unwrap());
    }
}

impl PafFile {
    pub fn new(filename: &str) -> Self {
        let mut query_names: Vec<String> = Vec::new();
        let mut target_names: Vec<String> = Vec::new();
        for_each_line_in_file(filename, |l: &str| {
            query_names.push(paf_query(l));
            target_names.push(paf_target(l));
        });
        query_names = query_names.into_iter().sorted().dedup().collect();
        target_names = target_names.into_iter().sorted().dedup().collect();

        let query_mphf = Mphf::new(1.7, &query_names);
        let target_mphf = Mphf::new(1.7, &target_names);
        let mut seen_queries = vec![false; query_names.len()];
        let mut seen_targets = vec![false; target_names.len()];
        let mut queries: Vec<AlignedSeq> = vec![AlignedSeq::new(); query_names.len()];
        let mut targets: Vec<AlignedSeq> = vec![AlignedSeq::new(); target_names.len()];
        for_each_line_in_file(filename, |l: &str| {
            let query_name: String = paf_query(l);
            let query_id = query_mphf.hash(&query_name) as usize;
            if !seen_queries[query_id] {
                seen_queries[query_id] = true;
                let mut query = &mut queries[query_id];
                query.name = query_name;
                query.length = paf_query_length(l);
            }
            let target_name: String = paf_target(l);
            let target_id = target_mphf.hash(&target_name) as usize;
            if !seen_targets[target_id] {
                seen_targets[target_id] = true;
                let mut target = &mut targets[target_id];
                target.name = target_name;
                target.length = paf_target_length(l);
            }
        });

        let mut targets_sort = targets.clone();
        targets_sort.sort_by(|a, b| b.length.partial_cmp(&a.length).unwrap());
        let mut target_idx: usize = 0;
        let mut target_offset: usize = 0;
        targets_sort.iter().for_each(|t| {
            let target_id = target_mphf.hash(&t.name) as usize;
            let mut target = &mut targets[target_id];
            target.rank = target_idx;
            target_idx += 1;
            target.offset = target_offset;
            target_offset += target.length;
        });
        let mut queries_sort = queries.clone();
        queries_sort.sort_by(|a, b| b.length.partial_cmp(&a.length).unwrap());
        let mut query_idx: usize = 0;
        let mut query_offset: usize = 0;
        queries_sort.iter().for_each(|q| {
            let query_id = query_mphf.hash(&q.name) as usize;
            let mut query = &mut queries[query_id];
            query.rank = query_idx;
            query_idx += 1;
            query.offset = query_offset;
            query_offset += query.length;
        });

        PafFile {
            filename: filename.to_string(),
            targets,
            queries,
            target_mphf,
            query_mphf,
            target_length: target_offset as f64,
            query_length: query_offset as f64,
        }
    }

    fn query_range(self: &PafFile, name: &str, start: usize, end: usize) -> (usize, usize) {
        //println!("query_range {} {}", start, end);
        let query_id = self.query_mphf.hash(&name.into()) as usize;
        //let length = self.query_length(query_id);
        let gstart = self.global_query_start(query_id);
        //println!("global query start {}", gstart);
        //println!("query length {}", length);
        let final_start = gstart + start;
        let final_end = gstart + end;
        (final_start, final_end)
    }
    fn target_range(self: &PafFile, name: &str, start: usize, end: usize) -> (usize, usize) {
        //println!("target_range {} {}", start, end);
        let target_id = self.target_mphf.hash(&name.into()) as usize;
        //let length = self.target_length(target_id);
        let gstart = self.global_target_start(target_id);
        //println!("global target start {}", gstart);
        let final_start = gstart + start;
        let final_end = gstart + end;
        (final_start, final_end)
    }
    fn global_query_start(self: &PafFile, idx: usize) -> usize {
        self.queries[idx].offset
    }
    fn query_length(self: &PafFile, idx: usize) -> usize {
        self.queries[idx].length
    }
    fn global_target_start(self: &PafFile, idx: usize) -> usize {
        self.targets[idx].offset
    }
    fn target_length(self: &PafFile, idx: usize) -> usize {
        self.targets[idx].length
    }
    fn global_start(self: &PafFile, line: &str, query_rev: bool) -> (usize, usize) {
        let query_id = self.query_mphf.hash(&paf_query(line)) as usize;
        let target_id = self.target_mphf.hash(&paf_target(line)) as usize;
        (
            self.global_target_start(target_id) + paf_target_begin(line),
            if query_rev {
                self.global_query_start(query_id) + paf_query_end(line)
            } else {
                self.global_query_start(query_id) + paf_query_begin(line)
            },
        )
    }

    fn for_each_match<F>(self: &PafFile, line: &str, mut func: F)
        where
            F: FnMut(char, usize, bool, usize, usize),
    {
        let query_rev = paf_query_is_rev(line);
        let (x, y) = self.global_start(line, query_rev);
        let mut target_pos = x;
        let mut query_pos = y;
        // find and walk the cigar string
        //println!("{}", line);
        let mut cigars = line
            .split('\t')
            .skip_while(|s| !s.starts_with("cg:Z:"))
            .map(|s| s.strip_prefix("cg:Z:").unwrap())
            .collect::<Vec<&str>>();

        // Compute the fake CIGAR, in case the real one is missing in the line
        let query_len = paf_query_end(line) - paf_query_begin(line);
        let target_len = paf_target_end(line) - paf_target_begin(line);
        let fake_cigar = format!(
            "{}M",
            if query_len < target_len {
                query_len
            } else {
                target_len
            }
        );
        if cigars.is_empty() {
            cigars.push(fake_cigar.as_str());
        }

        for cigar in cigars {
            //println!("{}", cigar);

            let mut ungapped_alignment_len: usize = 0;
            let mut query_delta: usize = 0;
            let mut target_delta: usize = 0;

            let mut first: usize = 0;
            for (i, b) in cigar.bytes().enumerate() {
                let c = b as char;
                //println!("{} {}", i, b as char);
                match c {
                    'M' | '=' | 'X' => {
                        let n = cigar[first..i].parse::<usize>().unwrap() as usize;
                        func(c, query_pos, query_rev, target_pos, n);
                        query_pos += if query_rev { 0 - n } else { n };
                        target_pos += n;
                        first = i + 1;

                        if ungapped_alignment_len > 0 && (query_delta > 0 || target_delta > 0) {
                            println!("{}\t{}\t{}", ungapped_alignment_len, target_delta, query_delta);

                            ungapped_alignment_len = 0;
                            target_delta = 0;
                            query_delta = 0;
                        }

                        ungapped_alignment_len += n;
                    }
                    'D' => {
                        let n = cigar[first..i].parse::<usize>().unwrap();
                        target_pos += n;
                        first = i + 1;

                        target_delta += n;
                    }
                    'I' => {
                        let n = cigar[first..i].parse::<usize>().unwrap();
                        query_pos += if query_rev { 0 - n } else { n };
                        first = i + 1;

                        query_delta += n;
                    }
                    _ => {}
                }
            }

            if ungapped_alignment_len > 0 {
                println!("{}", ungapped_alignment_len);
            }
        }
    }

    pub(crate) fn for_each_match_in_file<F>(self: &PafFile, mut func: F)
        where
            F: FnMut(char, usize, bool, usize, usize),
    {
        for_each_line_in_file(&self.filename, |line: &str| {
            /*let (x, y) = self.global_start(line, paf_query_is_rev(line));
            println!(
                "{} {} {} {} {} {}",
                paf_query(line),
                paf_query_begin(line),
                paf_target(line),
                paf_target_begin(line),
                x,
                y
            );*/

            self.for_each_match(line, |c, x, r, y, d| func(c, x, r, y, d));
            //println!();
        });
    }
    fn get_axes(self: &PafFile, major_axis: usize) -> (usize, usize) {
        //let max_length = cmp::max(self.query_length, self.target_length) as f64;
        //println!("query = {} target = {}", self.query_length, self.target_length);
        if self.target_length > self.query_length {
            let ratio = self.query_length as f64 / self.target_length as f64;
            (major_axis, (major_axis as f64 * ratio) as usize)
        } else {
            let ratio = self.target_length as f64 / self.query_length as f64;
            ((major_axis as f64 * ratio) as usize, major_axis)
        }
    }
    fn get_axes_zoom(
        self: &PafFile,
        major_axis: usize,
        zoom: ((usize, usize), (usize, usize)),
    ) -> (usize, usize) {
        let t_length = zoom.0.1 - zoom.0.0;
        let q_length = zoom.1.1 - zoom.1.0;
        if t_length > q_length {
            let ratio = q_length as f64 / t_length as f64;
            (major_axis, (major_axis as f64 * ratio) as usize)
        } else {
            let ratio = t_length as f64 / q_length as f64;
            ((major_axis as f64 * ratio) as usize, major_axis)
        }
    }
    fn project_xy(self: &PafFile, x: usize, y: usize, axes: (usize, usize)) -> (f64, f64) {
        //println!("axes {} {}", axes.0, axes.1);
        (
            axes.0 as f64 * (x as f64 / self.target_length as f64),
            axes.1 as f64 * (y as f64 / self.query_length as f64),
        )
    }
    fn project_xy_zoom(
        self: &PafFile,
        x: usize,
        y: usize,
        axes: (usize, usize),
        zoom: ((usize, usize), (usize, usize)),
    ) -> (f64, f64) {
        //println!("axes {} {}", axes.0, axes.1);
        (
            axes.0 as f64 * (((x as f64) - (zoom.0.0 as f64)) / ((zoom.0.1 - zoom.0.0) as f64)),
            axes.1 as f64 * (((y as f64) - (zoom.1.0 as f64)) / ((zoom.1.1 - zoom.1.0) as f64)),
        )
    }
    /*
        fn get_pixel(self: &PafFile, x: i64, y: i64, axes: (usize, usize)) -> usize {
            let (q, t) = axes;
            x * (self.query_length / q as f64).round() as usize * q + y * (self.target_length / t as f64).round() as usize
        }
    */
}
