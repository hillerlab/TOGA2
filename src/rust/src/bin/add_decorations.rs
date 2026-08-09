use clap::Parser;
use cubiculum::extract::extract::parse_bed;
use cubiculum::structs::structs::{BedEntry, Coordinates};
use fxhash::FxHashMap;
use phf::phf_map;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, stdout};
use std::path::Path;

use rust::read;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Given a Bed file and a TOGA2 mutation report file, prepare 
/// a decorator for the UCSC BigBed file\n
/// 
/// 
struct Args {

    /// path to the Bed12 file to decorate; 
    /// make sure that the contents correspond to those of the decorated BigBed file
    #[arg(long, short = 'b')]
    bed_file: String,

    ///  path to mutation file
    #[arg(long, short = 'm')]
    mutation_file: String,

    /// path to exon metadata table
    #[arg(long, short = 'e')]
    exon_meta: String,

    /// path to chromosome (contig, scaffold, etc.) size table
    #[arg(long, short = 'c')]
    chrom_sizes: String,

    /// path to the output file, if not set, the results are printed to standard output
    #[arg(long, short = 'o', default_value_t = String::from("stdout"))]
    output: String
}

const MUT2SHAPE: phf::Map<&'static str, &'static str> = phf_map!{
    "FS_DEL" => "Triangle",
    "FS_INS" => "InvTriangle",
    "STOP" => "Star",
    "SSMA" => "Square",
    "SSMD" => "Square"
};
const BLOCK_MUTS: (&str, &str) = ("BIG_DEL", "BIG_INS");

const ACCEPTOR: &str = "SSMA";
const DONOR: &str = "SSMD";
const PROJECTION: &str = "projection";

const REG_DEL: &str = "FS_DEL";
const REG_INS: &str = "FS_INS";
const BIG_DEL: &str = "BIG_DEL";
const BIG_INS: &str = "BIG_INS";

const DELETED: &str = "Deleted exon";
const MISSING: &str = "Missing exon";

const DELETED_REASON: &str = "Exon is deleted";
const MISSING_REASON : &str = "Exon is missing";

const MASKED: &str = "211,211,211,255";
const NOT_MASKED: &str = "198,0,9,255";
const GLYPH: &str = "glyph";
// static MUT2SHAPE: FxHashMap<&str, String> = FxHashMap::from_iter([
//     ("", String::from("")),
// ]);

enum MutationType {
    Interval,
    Point
}

#[derive(Clone)]
struct BlockMutation {
    chrom: String,
    start: u64,
    end: u64,
    exon: String,
    mut_type: String,
    color: &'static str
}

#[derive(Clone)]
struct PointMutation {
    chrom: String,
    start: u64,
    end: u64,
    exon: String,
    mut_type: String,
    color: &'static str
}

fn main() {
    let args = Args::parse();
    /// here go the sanity checks
    let mut proj2muts: FxHashMap<String, Vec<PointMutation>> = FxHashMap::default();
    // read the mutation file first
    let mutations = {
        let path = File::open(args.mutation_file).unwrap();
        Box::new(BufReader::new(path)) as Box<dyn BufRead>
    };
    for (i, line_) in mutations.lines().enumerate() {
        if let Ok(line) = line_ {
            let mut_comps: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
            if mut_comps.len() == 0 {continue};
            let proj: &str = mut_comps[0];
            // ignore fragmented projections; due to the current naming scheme, 
            // mutations cannot be attributed to fragments unambiguously
            if proj.contains(',') {continue};
            let mut_type: &str = mut_comps[7];
            // ignore mutations outside of the point mutation list§
            if !MUT2SHAPE.contains_key(mut_type) {continue};
            let mask_reason: &str = mut_comps[10];
            // ignore mutations corresponding to deleted and missing exons 
            if mask_reason == DELETED_REASON || mask_reason == MISSING_REASON {continue};
            let chrom: &str = mut_comps[4];
            // let start: u64 = mut_comps[5]
            //     .parse::<u64>()
            //     .expect(&format!(
            //         "Improper mutation file formatting at line {}: \"start\" field is not a valie integer", i+1
            //     ));
            let start: u64 = match mut_comps[5].parse::<u64>() {
                Ok(s) => {s},
                Err(_) => {
                    eprintln!("Improper mutation file formatting at line {}: \"start\" field is not a valid integer", i+1);
                    continue
                }
            };
            // let end: u64 = mut_comps[6]
            //     .parse::<u64>()
            //     .expect(&format!(
            //         "Improper mutation file formatting at line {}: \"end\" field is not a valie integer", i+1
            //     ));
            let end: u64 = match mut_comps[6].parse::<u64>() {
                Ok(e) => {e},
                Err(_) => {
                    eprintln!("Improper mutation file formatting at line {}: \"end\" field is not a valied integer", i+1);
                    continue
                }
            };
            let color: &str = match mut_comps[9] {
                "MASKED" => {MASKED},
                "NOT_MASKED" => {NOT_MASKED},
                _ => {
                    panic!(
                        "Improper mutation file formatting at line {}: \"is_masked\" field contains illegal value", i+1
                    )
                }
            };
            let exon: &str = mut_comps[1];
            let mut_entry: PointMutation = PointMutation { 
                chrom: chrom.to_string(), 
                start: start, 
                end: end, 
                exon: exon.to_string(),
                mut_type: mut_type.to_string(), 
                color: color 
            };
            proj2muts.entry(proj.to_string())
                .and_modify(|x| x.push(mut_entry.clone()))
                .or_insert(vec![mut_entry]);
        }
    }

    // parse the exon meta file
    let mut proj2exon2status: FxHashMap<String, FxHashMap<u16, bool>> = FxHashMap::default();
    for (i, line_) in read(args.exon_meta).lines().enumerate() {
        if let Ok(line) = line_ {
            let comps: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
            if comps[0] == PROJECTION {continue}
            let proj = comps[0];
            let exon_num = comps[1].parse::<u16>().expect(
                &format!(
                    "Improper formatting at exon meta file line {}; \"{}\" is not a valid exon number",
                    i + 1, comps[1]
                )
            );
            let status = comps[7] == "I";
            proj2exon2status
                .entry(proj.to_string())
                .and_modify(|x| {
                    x.insert(exon_num, status);
                })
                .or_insert(FxHashMap::from_iter(vec![(exon_num, status)]));
        }
    }

    // parse the chromosome size table
    let mut chrom_sizes: FxHashMap<String, u64> = FxHashMap::default();
    let chrom_sizes_file = {
        let path = File::open(args.chrom_sizes).unwrap();
        Box::new(BufReader::new(path)) as Box<dyn BufRead>
    };
    for (i, line_) in chrom_sizes_file.lines().enumerate() {
        if let Ok(line) = line_ {
            let chrom_sizes_comps: Vec<&str> = line.split('\t').collect::<Vec<&str>>();
            if chrom_sizes_comps.len() < 2 {
                panic!(
                    "Invalid chrom size file formatting at line {}; expected 2 fields, got {}",
                    i+1, chrom_sizes_comps.len()
                )
            }
            let size: u64 = chrom_sizes_comps[1]
                .parse::<u64>()
                .expect(
                    &format!(
                        "Invalid chrom size file formatting at line {}; value in the second field is not a valid integer",
                        chrom_sizes_comps[1]
                    )
                );
            chrom_sizes.insert(chrom_sizes_comps[0].to_string(), size);
        }
    }

    // create an output file handle
    let mut output_file = match args.output.as_str() {
        "stdout" => {Box::new(stdout()) as Box<dyn Write>},
        _ => {
            let path = Path::new(&args.output);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };

    // then, read the input Bed file; for each transcript encountered, 
    // retrieve mutation entries if present, prepare decoration entries,
    // and write them to the output file
    let bed_file = {
        let path = File::open(args.bed_file).unwrap();
        Box::new(BufReader::new(path)) as Box<dyn BufRead>
    };
    for line_ in bed_file.lines() {
        if let Ok(line) = line_ {
            let bed_entry: BedEntry = match parse_bed(line, 12, true) {
                Some(entry) => {entry},
                None => {continue}
            };
            let proj: &str = bed_entry.name().unwrap();
            // skip fragmented projections - TOGA2 currently does not support them
            if proj.contains(',') {continue};
            if let Some(muts) = proj2muts.get(proj) {
                // extract projection values needed for the mutation Bed record
                let strand: bool = bed_entry.strand().unwrap();
                let strand_lit: char = if strand {'+'} else {'-'};
                let tr_chrom: &str = bed_entry.chrom().unwrap();
                let chrom_size: u64 = *chrom_sizes
                    .get(tr_chrom)
                    .expect(&format!("Chromsome {} is missing from the chromosome size file", tr_chrom));
                let tr_start: &u64 = bed_entry.start().unwrap();
                let tr_end: &u64 = bed_entry.end().unwrap();
                let tr_coords: String = format!("{}:{}-{}:{}", tr_chrom, tr_start, tr_end, proj);
                'outer: for mutation in muts {
                    // ignore non-present exons
                    match mutation.exon.contains('_') {
                        true => {
                            let affected_exons = mutation
                                .exon
                                .split('_')
                                .map(|x| 
                                    x.parse::<u16>().expect(&format!("Invalid mutation exon value: {}", mutation.exon))
                                )
                                .collect::<Vec<u16>>();
                            for i in affected_exons[0]..*affected_exons.last().unwrap() {
                                let status = *proj2exon2status
                                    .get(proj)
                                    .expect(
                                        &format!("Unknown presence status for {}, exon {}", proj, i)
                                    )
                                    .get(&i)
                                    .expect(
                                        &format!("Unknown presence status for {}, exon {}", proj, i)
                                    );
                                if !status {continue 'outer;}
                            }
                        },
                        false => {
                            let exon = mutation.exon
                                .parse::<u16>()
                                .expect(&format!("Invalid mutation exon value: {}", mutation.exon));
                            let status = *proj2exon2status
                                .get(proj)
                                .expect(
                                    &format!("Unknown presence status for {}, exon {}", proj, exon)
                                )
                                .get(&exon)
                                .expect(
                                    &format!("Unknown presence status for {}, exon {}", proj, exon)
                                );
                            if !status {continue}
                        }
                    };
                    // point mutations have a total length of 1; 
                    // meanwhile, mutation files provided affected codon's coordinate;
                    // therefore, one has to first infer the exact coordinate
                    // let mut first_coord: u64  = if strand {mutation.end} else {mutation.start};
                    let mut first_coord: u64 = if &mutation.mut_type == ACCEPTOR {
                        if strand {mutation.end} else {mutation.start}
                    } else {
                        if strand {mutation.start} else {mutation.end}
                    };
                    // let mut first_coord: u64 = match mutation.exon.contains('_') {
                        // true => {
                        //     // mutation spans over multiple exons;
                        //     // acceptors should be mapped to mutation's end,
                        //     // the rest go to the first coordinate as always
                        //     if &mutation.mut_type == ACCEPTOR {
                        //         if strand {mutation.end} else {mutation.start}
                        //     } else {
                        //         if strand {mutation.start} else {mutation.end}
                        //     }
                        // },
                        // false => {
                        //     // mutation occurs within one exon;
                        //     // simply place it at the beginning of the codon
                        //     // NOTE: before v2.0.9e, donor mutations glyphs 
                        //     // were placed at the second position, creating an awkward gap
                        //     // if &mutation.mut_type == DONOR {
                        //     //     if strand {mutation.end} else {mutation.start}
                        //     // } else {
                        //     //     if strand {mutation.start} else {mutation.end}
                        //     // }
                        //     if strand {mutation.start} else {mutation.end}
                        // }
                    // };
                    let mut last_coord: u64 = first_coord + 1;
                    if last_coord >= chrom_size {
                        last_coord = chrom_size - 1;
                        if first_coord >= last_coord {first_coord = last_coord - 1}
                    }
                    let shape = MUT2SHAPE.get(&mutation.mut_type).unwrap();
                    let mut_name: String = format!(
                        "{}:{}-{}:{}:{}", 
                        tr_chrom, first_coord, last_coord, mutation.exon, mutation.mut_type
                    );
                    let decorator_line: String = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        mutation.chrom, first_coord, last_coord, mut_name, 0, strand_lit,
                        first_coord, last_coord, mutation.color, 1, 1, 0,
                        tr_coords, GLYPH, mutation.color, shape
                    );
                    if let Err(e) = output_file.write(decorator_line.as_bytes()) {
                        eprintln!("Failed to write the line: {}", e);
                    }
                };
            };
        }
    }
}