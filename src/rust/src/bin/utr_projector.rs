use chaintools::io::reader::Reader;
use clap::Parser;
use cubiculum::structs::structs::{BedEntry, Coordinates, Named, UtrSide};
use cubiculum::extract::extract::{BedFractionMode, extract_fraction, parse_bed, to_line};
use flate2::read::MultiGzDecoder;
use fxhash::{FxHashMap, FxHashSet};
use std::cmp::min;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::num::ParseIntError;
// use std::ops::Deref;
use std::path::Path;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Annotates untranslated regions (UTRs) by projecting reference annotation via genome alignment chains 
/// and adds the resulting UTR sequences to the TOGA2 output annotation
/// 
struct Args {
    /// Query annotation produced by TOGA2, in BED12 format
    #[arg(long, short = 'q')]
    query_annotation: String,

    /// Reference annotation used as TOGA2 input, in BED12 format. Note that UTR annotations 
    /// must be preserved for the module to annotate anything
    #[arg(long, short = 'r')]
    reference_annotation: String,

    /// Genome alignment chain file; might be gzip-compressed; for best performance, consider indexing it
    /// with chaintools indexer
    #[arg(long, short = 'c')]
    chain_file: String,

    /// Exon metadata file; used to infer loss status of the terminal exons
    #[arg(long, short = 'e')]
    exon_meta: String,

    /// A path to the output file; if set to `stdout` or not provided at all, 
    /// the results are printed to standard output
    #[arg(long, short = 'o', default_value_t = String::from("stdout"))]
    output: String,

    /// Absolute elongation threshold, in base pairs; unaligned UTR exon ends will not be extended farther 
    /// than this value
    #[arg(long, short = 'a', default_value_t = 3000)]
    abs_threshold: u64,

    /// Relative elongation threshold, times the exon length; unaligned UTR exon ends will not be extended farther 
    /// than this value
    #[arg(long, short = 'R', default_value_t = 2.5)]
    rel_threshold: f64,

    /// Deprecates unaligned sequence extrapolation
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    no_extrapolation: bool,

    /// If set, unprojected CDS-adjacent untranslated regions are added to the query projections
    /// based on their 
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    deduce_adjacent_regions: bool,

    /// If set, emulates the --deduce-adjacent-regions, except for 
    /// fixed value equal to --abs-threshold for the grafted UTR sequence 
    /// instead of the reference counterpart's length
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = false)]
    fixed_adjacent_regions: bool
}

fn read<T>(file: T) -> Box<dyn BufRead> 
where 
    T: AsRef<Path>
{
    if file.as_ref().extension().unwrap() == "gz" {
        let path = File::open(file).unwrap();
        return Box::new(BufReader::new(MultiGzDecoder::new(path))) as Box<dyn BufRead>
    } else {
        let path = File::open(file).unwrap();
        return Box::new(BufReader::new(path)) as Box<dyn BufRead>
    }
}

// fn read_gz<T>(file: T) -> Box<dyn BufRead>
// where
//     T: AsRef<Path>
// {
//     // let mut decoder = MultiGzDecoder::new();
//     //     decoder
//     //         .read_to_end(&mut args.exon_meta)
//     //         .with_context(|| format!("Failed to read file {:?}", file))?;
//     let path = File::open(file).unwrap();
//     Box::new(BufReader::new(MultiGzDecoder::new(path))) as Box<dyn BufRead>
// }

fn get_tr(proj: &str) -> String {
    let proj_comps = proj.split('#').collect::<Vec<&str>>();
    format!("{}#{}", proj_comps[0], proj_comps[1])
}

fn get_chain_id(proj: &str) -> Result<u32, ParseIntError> {
    proj.split('#').last().unwrap().parse::<u32>()
}

fn main() {
    let args = Args::parse();

    // create the hashmaps for further use
    // chain-to-transcript storage; for each chain, UTRs for the following transcripts will be further annotated
    let mut chain2trs: FxHashMap<u32, Vec<String>> = FxHashMap::default();
    // transcript-to-utr-exons; likely will be easier to store them in memory rather than recalulcate them every time
    let mut ref_tr2utrs: FxHashMap<String, Vec<BedEntry>> = FxHashMap::default();
    // transcript-to-strand storage; useful for quick data access when grafting
    // let mut ref_tr2strand: FxHashMap<String, bool> = FxHashMap::default();
    // UTR exon-to-maximal-relative-size storage; speeds up 
    let mut utr2rel_size: FxHashMap<String, u64> = FxHashMap::default();
    // projection-to-bed-entry; stores the projection data for annotation and output
    let mut query_tr2bed: FxHashMap<String, BedEntry> = FxHashMap::default();
    // feature storage; likely better than UtrBlock storage
    let mut ref_utr2features: FxHashMap<String, (bool, UtrSide)> = FxHashMap::default();
    // query chrom size dict for fast access at the final step of the pipeline
    let mut query_chr2size: FxHashMap<u32, u64> = FxHashMap::default();

    // open the output file
    let mut output_file = match args.output.as_str() {
        "stdout" => {Box::new(io::stdout()) as Box<dyn Write>},
        _ => {
            let path = Path::new(&args.output);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };

    // create a storage of CDS-adjacent UTR portions; if those are not projected properly,
    // substitutes equal in length to the reference exons will be added to the projections
    let mut tr2adj3: FxHashMap<String, u64> = FxHashMap::default();
    let mut tr2adj5: FxHashMap<String, u64> = FxHashMap::default();
    let mut adj3_appended: FxHashSet<String> = FxHashSet::default();
    let mut adj5_appended: FxHashSet<String> = FxHashSet::default();


    // read the reference BED file, infer the UTRs, and populate the transcript-to-UTRs dictionary
    let r_bed_start = Instant::now();
    for line in read(args.reference_annotation).lines() {
        if let Ok(x) = line {
            // turn record into a BedEntry item
            let bed_record = parse_bed(x, 12, false).unwrap();
            let tr = bed_record.name().unwrap();
            // get the UTR blocks // TODO: Add a line->fraction extractor function to cubiculum
            let utr = extract_fraction(
                &bed_record, BedFractionMode::Utr, false
            ).expect(&format!("UTR extraction from the entry {} failed", tr));
            if let None = utr {continue}
            // ref_tr2strand.insert(tr.to_string(), bed_record.strand().unwrap());
            let utr_blocks = utr.unwrap().to_blocks().unwrap();
            let cds_start = bed_record.thick_start().unwrap();
            let cds_end = bed_record.thick_end().unwrap();
            // convert each UTR block into a UTR object // TODO: Add a direct UTR extractor function to cubiculum
            for (i, mut utr_record) in utr_blocks.into_iter().enumerate() {
                let utr_start = utr_record.start().expect("Failed to infer start coordinate for a UTR block");
                let utr_end = utr_record.end().expect("Failed to infer end coordinate for a UTR block");
                let mut side: UtrSide = UtrSide::FivePrime;
                let mut is_adjacent = false;
                if *utr_end <= cds_start {
                    side = match utr_record.strand() {
                        Some(true) => {UtrSide::FivePrime},
                        Some(false) => {UtrSide::ThreePrime},
                        None => {continue} // TODO: Temporary workaround
                    };
                    // utr_record.set_side(side);
                    if *utr_end == cds_start {
                        is_adjacent = true;
                        if args.deduce_adjacent_regions || args.fixed_adjacent_regions {
                            // record that the transcript has an upstream-adjacent UTR block
                            let adj_utr_len = utr_record.length().unwrap();
                            match side {
                                UtrSide::ThreePrime => {tr2adj3.insert(tr.clone(), adj_utr_len);},
                                UtrSide::FivePrime => {tr2adj5.insert(tr.clone(), adj_utr_len);}
                            }
                        }
                        // utr_record.set_adjacency(is_adjacent);
                    }
                } else if *utr_start >= cds_end {
                    side = match utr_record.strand() {
                        Some(true) => {UtrSide::ThreePrime},
                        Some(false) => {UtrSide::FivePrime}
                        None => {continue} // TODO: Temporary workaround
                    };
                    if *utr_start == cds_end {
                        // utr_record.set_adjacency(is_adjacent);
                        is_adjacent = true;
                        if args.deduce_adjacent_regions || args.fixed_adjacent_regions {
                            // record that the transcript has an upstream-adjacent UTR block
                            let adj_utr_len = utr_record.length().unwrap();
                            match side {
                                UtrSide::ThreePrime => {tr2adj3.insert(tr.clone(), adj_utr_len);},
                                UtrSide::FivePrime => {tr2adj5.insert(tr.clone(), adj_utr_len);}
                            }
                        }
                    }
                }
                // utr_record.set_adjacency(is_adjacent);
                let utr_exon_name = format!("{}|{}", bed_record.name().unwrap(), i);
                ref_utr2features.insert(
                    utr_exon_name.clone(),
                    (is_adjacent, side)
                );
                // utr_record.set_side(side);
                utr_record.update_name(&format!("{}|{}", bed_record.name().unwrap(), i));
                ref_tr2utrs
                    .entry( tr.clone())
                    .and_modify(|x| x.push(utr_record.clone()))
                    .or_insert(vec![utr_record.clone()]);
                let rel_size = utr_record
                    .length()
                    .expect("Failed to infer UTR exon size") as f64 * args.rel_threshold;
                utr2rel_size.insert(utr_exon_name.clone(), rel_size as u64);
            }
        }
    } 
    println!("Reference bed parsing: time: {:?}", r_bed_start.elapsed());

    // read the query BED file, populate the chain-to-transcript and projection-to-bed dictionaries
    let q_bed_start = Instant::now();
    for line in read(args.query_annotation).lines() {
        if let Ok(x) = line {
            let query_bed = parse_bed(
                x, 12, false
            ).expect("Failed to parse the BED record in the query annotation");
            let comps: Vec<&str> = query_bed.name().unwrap().split('#').collect();
            let last_field: usize = comps.len() - 1;
            let chain_field: usize = if comps[last_field] == "retro" {last_field - 1} else {last_field};
            let chain_num = match comps[chain_field].parse::<u32>(){ // NOTE: picking the third rather than the last item to accommodate for `#retro` and other postfixes
                Ok(num) => {num},
                Err(_) => {
                    if let Err(e) = writeln!(output_file, "{}", to_line(&query_bed, 12).unwrap()) {
                        eprintln!("Failed to write the line: {}", e);
                    }
                    continue;
                }
            };
            let tr_name = comps[..chain_field].join("#");
            let proj_name = comps[..chain_field+1].join("#");
            chain2trs
                .entry(chain_num)
                .and_modify(|x| x.push(tr_name.clone()))
                .or_insert(vec![tr_name.clone()]);
            query_tr2bed.insert(proj_name, query_bed);
        }
    } 
    // println!("{:#?}", query_tr2bed);
    println!("Query bed parsing: time: {:?}", q_bed_start.elapsed());

    // read the exon meta file
    let mut proj2first_exon: FxHashMap<String, bool> = FxHashMap::default();
    let mut proj2last_exon: FxHashMap<String, bool> = FxHashMap::default();
    let exon_meta_start = Instant::now();
    {
        let mut proj2last_num: FxHashMap<String, u16> = FxHashMap::default();
        for line in read(args.exon_meta).lines() {
            if let Ok(parsed_line) = line {
                let exon_meta_data: Vec<&str> = parsed_line.split('\t').collect();
                if exon_meta_data.len() == 0 {continue}
                if exon_meta_data.len() < 8 {
                    panic!("Exon metadata file contains less than eight columns")
                }
                let proj = exon_meta_data[0];
                if proj == "projection" {continue}
                let exon_status = exon_meta_data[7]
                    .chars()
                    .nth(0)
                    .expect("Failed to parse exon status") == 'I';
                let exon_num = exon_meta_data[1].parse::<u16>().expect("Failed to parse exon number");
                if exon_num == 1 {
                    proj2first_exon.insert(proj.to_string(), exon_status);
                }
                let curr_last_exon = *proj2last_num.get(&proj.to_string()).unwrap_or(&0);
                if exon_num > curr_last_exon {
                    proj2last_num.insert(proj.to_string(), exon_num);
                    proj2last_exon.insert(proj.to_string(), exon_status);
                }
            }
        }
    }

    println!("Exon metadata file parsing time: {:?}", exon_meta_start.elapsed());

    // read the chain file; likely only a narrow selection of chains is required
    // TODO: Ideally read one chain at a time;
    let chain_start = Instant::now();
    let chains = Reader::from_file(args.chain_file).expect("Failed to read the chain file");
    println!("Chain parsing time: {:?}", chain_start.elapsed());


    // showtime
    // for each chain, get all the UTRs to project
    let chain_iteration_start = Instant::now();
    for (chain_id, chain) in chains.iter() {
        // if *chain_id != 120648 {continue;}
        // if *chain_id != 17199 {continue}
        // if *chain_id != 261 {continue}
        let transcripts = match chain2trs.get(chain_id) {
            Some(x) => {x},
            None => {continue}
        };
        query_chr2size.insert(*chain_id, chain.query.size);
        let mut utrs_for_the_chain = Vec::new();
        for tr in transcripts {
            match ref_tr2utrs.get(tr) {
                Some(x) => utrs_for_the_chain.extend(x),
                None => {continue}
            };
        }
        if utrs_for_the_chain.len() == 0 {continue}
         // and, well, project them
        let projected_utrs = chain.map_through_(
            &mut utrs_for_the_chain,
            // 0,
            // 0.0,
            // if args.no_extrapolation {0} else{args.abs_threshold}, 
            // if args.no_extrapolation {0.0} else {args.rel_threshold}, 
            !args.no_extrapolation,
            true
        ).expect("UTR projection failure!");
        // let mut sorted_intervals: Vec<Interval>  = projected_utrs.into_values().collect();
        // if sorted_intervals.len() == 0 {continue}
        // sorted_intervals.sort_by(
        //     |a, b| if a.start().unwrap() == b.start().unwrap() {
        //         a.end().unwrap().cmp(&b.end().unwrap())
        //     } else {
        //         a.start().unwrap().cmp(&b.start().unwrap())
        //     }
        // );
        for (utr_name, utr_proj) in projected_utrs {
        // for utr_proj in sorted_intervals {
            if let (None, None) = (utr_proj.start(), utr_proj.end()) {continue}
            let rel_size = utr2rel_size
                .get(utr_name)
                .expect(&format!("Unknown reference UTR exon size for {}", utr_name));
            let exon_len = utr_proj
                .length()
                .expect(
                    &format!("Unknown query UTR exon size for {}, chain {}", utr_name, chain_id)
                );
            if exon_len > *rel_size && exon_len > args.abs_threshold {continue}
            let tr_name = utr_name.split('|')
                .into_iter()
                .rev()
                .last()
                .unwrap(); // TODO: Does not look safe
            let proj_name = &format!("{}#{}", tr_name, chain_id);
            if proj_name == "XM_059874414#LOC132342188#8" {
                println!("utr_proj={:#?}", utr_proj);
            }
            let query_strand = query_tr2bed.get(proj_name).unwrap().strand().unwrap();
            let (is_adjacent, side) = ref_utr2features.get(utr_name).unwrap();
            // do not add 5'-UTRs if the first exon was not annotated
            let first_exon_present = proj2first_exon
                .get(proj_name)
                .expect(&format!("Projection {} is missing from exon metadata file", proj_name));
            if *side == UtrSide::FivePrime && !*first_exon_present {continue}
            // likewise, do not add 3'-UTRs if the last exon was not annotated
            let last_exon_present = proj2last_exon
                .get(proj_name)
                .expect(&format!("Projection {} is missing from exon metadata file", proj_name));
            if *side == UtrSide::ThreePrime && !*last_exon_present {continue}
            let (append_upstream, append_downstream) = match (is_adjacent, side, query_strand) {
                (false, _, _) => {(false, false)},
                (true, UtrSide::FivePrime, true) =>  {(true, false)},
                (true, UtrSide::FivePrime, false) => {(false, true)},
                (true, UtrSide::ThreePrime, true) => {(false, true)},
                (true, UtrSide::ThreePrime, false) => {(true, false)}
            };
            // println!("utr_proj={:#?}, append_upstream={}, append_downstream={}", utr_proj, append_upstream, append_downstream);
            // last size check, this time for adjacent UTR exons; append their expected size
            // first, for upstream-appended UTRs
            if append_upstream {
                let upstream_extend = query_tr2bed
                    .get(proj_name)
                    .unwrap()
                    .start()
                    .unwrap()
                    .checked_sub(*utr_proj.end().unwrap())
                    .unwrap_or(0);
                if upstream_extend > args.abs_threshold && upstream_extend > *rel_size { continue}
            }
            // same for downstream-appended exons
            if append_downstream {
                let downstream_extend = utr_proj
                    .start()
                    .unwrap()
                    .checked_sub(
                        *query_tr2bed.get(proj_name).unwrap().end().unwrap()
                    )
                    .unwrap_or(0);
                if downstream_extend > args.abs_threshold && downstream_extend > *rel_size {continue}
            }
            if proj_name == "XM_059874414#LOC132342188#8" {
                println!("append_upstream={:?}, append_downstream={:?}", append_upstream, append_downstream);
            }
            query_tr2bed
                .entry(proj_name.to_string())
                .and_modify(|x| {
                    x.graft(
                        utr_proj, 
                        true, 
                        true, 
                        true, 
                        false, 
                        append_upstream, 
                        append_downstream
                    );
                }
            );
            // println!("after_grafting={:#?}", to_line(
            //     query_tr2bed.get(proj_name).unwrap(), 12
            // ).unwrap()
            // );
            // record the adjacent UTR portion being projected and successfully grafted
            if *is_adjacent {
                match side {
                    UtrSide::ThreePrime => {adj3_appended.insert(proj_name.to_string());},
                    UtrSide::FivePrime => {adj5_appended.insert(proj_name.to_string());}
                }
            }
        }
    }

    // let out_path = Path::new(&args.output);
    // let output = Box::new(File::create(&out_path).unwrap()) as Box<dyn Write>;
    for (mut name, mut val) in query_tr2bed {
        // if unconditional adjacent UTR portion deductino has been enabled, 
        // check if the respective UTR portions
        if name.contains("#retro") {
            name = name.replace("#retro", "");
        }
        if args.deduce_adjacent_regions || args.fixed_adjacent_regions {
            let tr_name = get_tr(&name);
            let chain_id = get_chain_id(&name);
            let query_size = match chain_id {
                Ok(x) => {*query_chr2size.get(&x).unwrap()},
                Err(..) => {0}
            };
            if tr2adj3.contains_key(&tr_name) && !adj3_appended.contains(&name) {
                println!("Adding an unprojected 3'-adjacent UTR to projection {}", name);
                let utr_size = match args.deduce_adjacent_regions {
                    true => {
                        tr2adj3
                            .get(&tr_name)
                            .expect(&format!("Unrecorded 3'-adjacent UTR size for transcript {}", tr_name))
                        },
                    false => {&args.abs_threshold}
                };
                let strand = val.strand().unwrap();
                let (start, end) = match strand {
                    true => {
                        let thick_end = val.thick_end().unwrap();
                        (thick_end, min(thick_end + *utr_size, query_size))
                    },
                    false => {
                        let thick_start = val.thick_start().unwrap();
                        (thick_start.checked_sub(*utr_size).unwrap_or(0), thick_start)
                    }
                };
                let terminal_exon_present = match strand {
                    true => {proj2last_exon.get(&name)},
                    false => {proj2first_exon.get(&name)}
                }.expect(&format!("Projection {} is missing from exon metadata file", name));
                if *terminal_exon_present {
                    let adjacent_graft = BedEntry::bed3(String::new(), start, end);
                    if name == "XM_059874414#LOC132342188#8" {
                        println!("adjacent_graft_3'={:#?}", adjacent_graft);
                    }
                    val.graft(
                        adjacent_graft, 
                        true, 
                        false, 
                        true, 
                        false, 
                        !strand, 
                        strand
                    );
                }
            }
            if tr2adj5.contains_key(&tr_name) && !adj5_appended.contains(&name) {
                println!("Adding an unprojected 5'-adjacent UTR to projection {}", name);
                // let utr_size = tr2adj5
                //     .get(&tr_name)
                //     .expect(&format!("Unrecorded 5'-adjacent UTR size for transcript {}", tr_name));
                let utr_size = match args.deduce_adjacent_regions {
                    true => {
                        tr2adj5
                            .get(&tr_name)
                            .expect(&format!("Unrecorded 5'-adjacent UTR size for transcript {}", tr_name))
                        },
                    false => {&args.abs_threshold}
                };
                let strand = val.strand().unwrap();
                let (start, end) = match strand {
                    true => {
                        let thick_start = val.thick_start().unwrap();
                        (thick_start.checked_sub(*utr_size).unwrap_or(0), thick_start)
                    },
                    false => {
                        let thick_end = val.thick_end().unwrap();
                        (thick_end, min(thick_end + *utr_size, query_size))
                    }
                };
                let terminal_exon_present = match strand {
                    true => {proj2first_exon.get(&name)},
                    false => {proj2last_exon.get(&name)}
                }.expect(&format!("Projection {} is missing from exon metadata file", name));
                if *terminal_exon_present {
                    let adjacent_graft = BedEntry::bed3(String::new(), start, end);
                    val.graft(
                        adjacent_graft, 
                        true, 
                        false, 
                        true, 
                        false, 
                        strand, 
                        !strand
                    );
                }
            }
        }
        let result_line = to_line(&val, 12).expect(
            &format!("Failed to convert entry for {} into a BED line", name)
        );
        if let Err(e) = writeln!(output_file, "{}", result_line) {
            eprintln!("Failed to write the line: {}", e);
        }
    }
    println!("Chain dict iteration time: {:?}", chain_iteration_start.elapsed());
}