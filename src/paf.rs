use boomphf::Mphf;
use itertools::Itertools;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use flate2::read::MultiGzDecoder;

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

pub fn paf_query(line: &str) -> String {
    line.split('\t').next().unwrap().into()
}

pub fn paf_query_length(line: &str) -> usize {
    line.split('\t').nth(1).unwrap().parse::<usize>().unwrap()
}

pub fn paf_query_begin(line: &str) -> usize {
    line.split('\t').nth(2).unwrap().parse::<usize>().unwrap()
}

pub fn paf_query_end(line: &str) -> usize {
    line.split('\t').nth(3).unwrap().parse::<usize>().unwrap()
}

pub fn paf_query_is_rev(line: &str) -> bool {
    line.split('\t').nth(4).unwrap() == "-"
}

pub fn paf_target(line: &str) -> String {
    line.split('\t').nth(5).unwrap().into()
}

pub fn paf_target_length(line: &str) -> usize {
    line.split('\t').nth(6).unwrap().parse::<usize>().unwrap()
}

pub fn paf_target_begin(line: &str) -> usize {
    line.split('\t').nth(7).unwrap().parse::<usize>().unwrap()
}

pub fn paf_target_end(line: &str) -> usize {
    line.split('\t').nth(8).unwrap().parse::<usize>().unwrap()
}

fn _for_each_line_in_file(paf_filename: &str, mut callback: impl FnMut(&str)) {
    let file = File::open(paf_filename).unwrap();
    // Determine if the file is gzipped based on its extension
    let gzipped = Path::new(paf_filename)
        .extension()
        .map_or(false, |ext| ext == "gz");

    // Create a dynamic reader based on the file type
    let box_reader: Box<dyn Read> = if gzipped {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let reader = BufReader::new(box_reader);
    for line in reader.lines() {
        callback(&line.unwrap());
    }
}

impl PafFile {
    pub fn new(filename: &str) -> Self {
        let mut query_names: Vec<String> = Vec::new();
        let mut target_names: Vec<String> = Vec::new();
        _for_each_line_in_file(filename, |l: &str| {
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
        _for_each_line_in_file(filename, |l: &str| {
            let query_name: String = paf_query(l);
            let query_id = query_mphf.hash(&query_name) as usize;
            if !seen_queries[query_id] {
                seen_queries[query_id] = true;
                let query = &mut queries[query_id];
                query.name = query_name;
                query.length = paf_query_length(l);
            }
            let target_name: String = paf_target(l);
            let target_id = target_mphf.hash(&target_name) as usize;
            if !seen_targets[target_id] {
                seen_targets[target_id] = true;
                let target = &mut targets[target_id];
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
            let target = &mut targets[target_id];
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
            let query = &mut queries[query_id];
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
    pub(crate) fn global_start(self: &PafFile, line: &str, query_rev: bool) -> (usize, usize) {
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

    pub(crate) fn for_each_line_in_file<F>(self: &PafFile, mut func: F)
        where
            F: FnMut(&str),
    {
        _for_each_line_in_file(&self.filename, |line: &str| {
            func(line);
        });
    }
}
