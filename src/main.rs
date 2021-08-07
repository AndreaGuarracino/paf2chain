use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;
//use std::cmp;

use boomphf::*;

use rgb::*;
extern crate line_drawing;
use line_drawing::XiaolinWu;

use itertools::Itertools;

extern crate clap;
use clap::{App, Arg}; //, SubCommand};

mod paf;
use crate::paf::PafFile;

fn main() {
    let matches = App::new("paf2chain")
        .version("0.1.0")
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
        .arg(
            Arg::with_name("OUTPUT")
                .takes_value(true)
                .short("o")
                .long("output")
                .help("save the CHAIN format output to this file."),
        )
        .get_matches();

    let filename = matches.value_of("INPUT").unwrap();
    let paf = PafFile::new(filename);

    let default_output = format!("{}.chain", filename);
    let output_chain = matches.value_of("output").unwrap_or(&default_output);

    // // draw our grid
    // paf.targets.iter().for_each(|target| {
    //     if target.offset > 0 {
    //         let start = get_coords(target.offset, 0);
    //         let end = get_coords(target.offset, paf.query_length as usize);
    //         for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
    //             if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
    //                 let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
    //                 pixels[i] = get_color(val * 0.2);
    //             }
    //         }
    //     }
    // });
    // paf.queries.iter().for_each(|query| {
    //     if query.offset > 0 {
    //         let start = get_coords(0, query.offset);
    //         let end = get_coords(paf.target_length as usize, query.offset);
    //         for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
    //             if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
    //                 let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
    //                 pixels[i] = get_color(val * 0.2);
    //             }
    //         }
    //     }
    // });

    // for each match, we draw a line on our raster using Xiaolin Wu's antialiased line algorithm
    let draw_match = |_c, x: usize, rev: bool, y: usize, len: usize| {
        // //println!("draw_match {} {} {} {} {}", _c, x, rev, y, len);
        // let start = get_coords(x, y);
        // let end = get_coords(x + if rev { 0 - len } else { len }, y + len);
        // /*
        // println!(
        //     "start and end ({} {}) ({} {})",
        //     start.0, start.1, end.0, end.1
        // );
        //  */
        // for ((j, i), val) in XiaolinWu::<f64, i64>::new(start, end) {
        //     //println!("checking pixel {} {} {}", i, j, val);
        //     if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
        //         //println!("drawing pixel {} {} {}", i, j, val);
        //         let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
        //         pixels[i] = get_color(val);
        //     }
        // }
    };
    paf.for_each_match_in_file(draw_match);

    let path = &Path::new(output_chain);

    // encode_file takes the path to the image, a u8 array,
    // the width, the height, the color mode, and the bit depth
    // if let Err(e) = lodepng::encode_file(path, &raw, axes.0, axes.1, lodepng::ColorType::RGB, 8) {
    //     panic!("failed to write png: {:?}", e);
    // }
}
