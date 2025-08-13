use anyhow::Context;
use clap::Parser;
use cubiculum::structs::structs::{BedEntry, Coordinates, Interval, Named};
use cubiculum::extract::extract::{BedFractionMode, extract_fraction, parse_bed};
use cubiculum::merge::merge::{discrete_interval_map, intersection};
use rust::get_cds_boundaries_;
use fxhash::{FxHashMap, FxHashSet};
use std::cmp::min;
use std::io::{self, BufRead, BufReader, Write};
use std::fs::File;
use std::path::Path;
use std::time::Instant;

const HEADER_LINE: &str = "transcript\tgene_overs\tchain\tsynt\tgl_score\tgl_exo\tchain_len\texon_qlen\t\
loc_exo\texon_cover\tintr_cover\tgene_len\tex_num\tex_fract\tintr_fract\tflank_cov\t\
clipped_exon_qlen\tclipped_intr_cover\n";


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Extracts feature from the chain-transcript pairs for future classification
/// 
struct Args {
    /// Reference annotation used as TOGA2 input, in BED12 format. Note that UTR annotations 
    /// must be preserved for precise feature estimation
    #[arg(long, short = 'b', default_value_t = String::from("/beegfs/projects/project-ymalovichko/toga_extension/duplication_tracing/TOGA2.0_tests/TREE_TESTS/mm10_iqtree_reduced_50/tmp/input_data/reference_annotation.bed"))]//String::from("/projects/hillerlab/genome/gbdb-HL/hg38/TOGA2/currentAnnotation/hg38.toga.transcripts.bed"))]
    bed_file: String,

    /// Genome alignment chain file; might be gzip-compressed; for best performance, consider indexing it
    /// with chaintools indexer (TODO: Chain indexing currently does not work properly)
    #[arg(long, short = 'c', default_value_t = String::from("/beegfs/projects/project-ymalovichko/toga_extension/duplication_tracing/TOGA2.0_tests/TREE_TESTS/mm10_iqtree_reduced_50/tmp/input_data/genome_alignment.chain"))]
    chain_file: String,

    /// [Optional] A two-column tab-separated file containing gene-to-transcript mapping, 
    /// in `long data` format. Used for synteny calculation
    #[arg(long, short = 'i')]
    isoform_file: Option<String>,

    /// A path to the output file; if set to `stdout` or not provided at all, 
    /// the results are printed to standard output
    #[arg(long, short = 'o', default_value_t = String::from("stdout"))]
    output: String,
}

fn main() {
    let args = Args::parse();
    // dummy values
    let chain_file = args.chain_file;
    // let chain_file = "/beegfs/projects/project-ymalovichko/toga_extension/duplication_tracing/TOGA2.0_tests/TREE_TESTS/HLanoCau1_iqtree_reduced_50/tmp/input_data/genome_alignment.chain";
    let bed_file = args.bed_file;

    // open the output file
    let mut output_file = match args.output.as_str() {
        "stdout" => {Box::new(io::stdout()) as Box<dyn Write>},
        _ => {
            let path = Path::new(&args.output);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };

    let program_start = Instant::now();

    // workflow

    // 1) parse bed
    let mut chrom2trs: FxHashMap<String, Vec<BedEntry>> = FxHashMap::default();
    let input_file = {
        let path = File::open(bed_file).unwrap();
        Box::new(BufReader::new(path)) as Box<dyn BufRead>
    };

    let now = Instant::now();
    for line_ in input_file.lines() {
        if let Ok(line) = line_ {
            let record: Option<BedEntry> = parse_bed(line, 12, false);
            if let Some(fraction) = record {
                chrom2trs
                    .entry(fraction.chrom().unwrap().clone())
                    .or_insert(Vec::new())
                    .push(fraction);
            }
        } else {println!("")}
    }
    let elapsed = now.elapsed();
    println!("BED parsing duration time: {:?}", elapsed);
    
    // 2) aggregate genes by chromosomes, then sort by coding sequence coordinates
    // since the transcripts are stored chromosome-wise, create a shortcut {transcript:index} dictionary
    let mut tr2index: FxHashMap<String, usize> = FxHashMap::default();
    let now = Instant::now();
    for (_, trs) in chrom2trs.iter_mut() {
        trs.sort_by(
            |a, b| if a.thick_start().unwrap() == b.thick_start().unwrap() {
                a.thick_end().unwrap().cmp(&b.thick_end().unwrap())
            } else {
                a.thick_start().unwrap().cmp(&b.thick_start().unwrap())
            }
            );
        for (i, tr) in trs.iter().enumerate() {
            tr2index.insert(
                tr.name().unwrap().clone(),
                i
            );
        }
    }
    let tr2index = tr2index;
    let elapsed = now.elapsed();
    println!("BED sorting duration time: {:?}", elapsed);

    // 3) read the chain file
    let now = Instant::now();
    // let chainmap = chaintools::io::reader::Reader::extract_ix(chain_file, None).with_context(||
    //     {"Failed to read the chain file"}
    // ); // TODO: Fails for certain transcripts!
    let chainmap = chaintools::io::reader::Reader::from_file(chain_file).with_context(||
        {"Failed to read the chain file"}
    ).unwrap();
    let elapsed = now.elapsed();
    println!("Chain reading duration time: {:?}", elapsed);

    // 4) if provided, parse the transcript-to-chain mapping
    let mut tr2gene: FxHashMap<String, String> = FxHashMap::default();
    if let Some(x) = &args.isoform_file {
        let isoform_file = {
            let path = File::open(x).unwrap();
            Box::new(BufReader::new(path)) as Box<dyn BufRead>
        };
        for (i, line_) in isoform_file.lines().enumerate() {
            if let Ok(line) = line_ {
                let contents = line.split("\t").collect::<Vec<&str>>();
                if contents.len() == 0 {continue}
                if contents.len() != 2 {
                    panic!(
                        "Wrong format encountered in the isoforms file at line {}; expected 2 columns, got {}",
                        i + 1, contents.len()
                    )
                }
                tr2gene.insert(
                    contents[1].to_string(),
                    contents[0].to_string()
                );
            }
        }
    }

    // write the header line to the output file
    if let Err(e) = output_file.write(HEADER_LINE.as_bytes()) {
        eprintln!("Failed to write the header line: {}", e);
    }


    // create the storage objects for transcript features that will be potentially utilized  more than once
    let mut transcript_length: FxHashMap<String, u64> = FxHashMap::default();
    let mut exon_number: FxHashMap<String, u16> = FxHashMap::default();
    let mut cds_length: FxHashMap<String, u64> = FxHashMap::default();
    let mut intron_length: FxHashMap<String, u64> = FxHashMap::default();
    let mut cds_intron_length: FxHashMap<String, u64> = FxHashMap::default();
    // let utr_length: FxHashMap<&str, u64> = FxHashMap::default();

    let now = Instant::now();
    let mut i: u64 = 0;
    let mut tr_counter: usize = 0;
    let mut tr_set: FxHashSet<String> = FxHashSet::default();
    for chain_id in chainmap.keys() {
        // if *chain_id != 1 {continue}
        let chain = chainmap
            .map
            .get(chain_id).
            unwrap();
        let transcripts = match chrom2trs.get_mut(&chain.refs.chr) {
            Some(x) => {x},
            None => continue
        };
        // define which transcripts are intersected by the chain
        let mut intersected = chain.intersect_to_cds_vector(transcripts, true);
        // let elapsed = now.elapsed();
        // println!("Transcript intersection duration time: {:?}", elapsed);
        if *chain_id <= 10 {
            println!("Intersected by chain {}: {}", chain_id, intersected.len());
        }
        // if chain.refs.chr == "chr9" {
        //     println!("Number of transcritps for chr9: {}", transcripts.len());
        // }
        // proceed if the chain has no transcripts intersected
        if intersected.len() == 0 {continue};

        // otherwise, start gathering the stats
        let tr_num = intersected.len();
        let chain_aln_score = chain.score;
        let ref_chain_len = chain.refs.end - chain.refs.start;
        let query_chain_len = (chain.query.end - chain.query.start) as f64;
        let query_chain_start: u64 = if chain.query.strand == '+' {chain.query.start} else {
            chain.query.size - chain.query.end
        };
        let query_chain_end: u64 = if chain.query.strand == '+' {chain.query.end} else {
            chain.query.size - chain.query.start
        };
        // get Aligning fraction (A); this is total aligned sequence minus the aligned UTR fraction
        let aln_block_sum = chain.alignment_sum();

        // create counters for chain-level coverage stats
        let mut cds_coverage: u64 = 0;
        let mut utr_coverage: u64 = 0;
        let mut intron_coverage: u64 = 0;
        // let mut clipped_cds_coverage: u64 = 0;
        // let mut clipped_intron_coverage: u64 = 0;

        // for synteny estimation, save the names of the transcrips 
        // which have at least one of their coding exons covered
        let mut syntenic_transcripts: FxHashSet<String> = FxHashSet::default();

        // record the query chromosome size; will come in handy later
        let ref_chrom_size = chain.refs.size;
        let ref_chrom = chain.refs.chr.clone();

        // create the transcriptwise storage hashmaps
        let mut tr2cds_cov: FxHashMap<String, u64> = FxHashMap::default();
        let mut tr2intron_cov: FxHashMap<String, u64> = FxHashMap::default();
        let mut tr2clipped_intron_cov: FxHashMap<String, u64> = FxHashMap::default();
        let mut tr2flank_cov: FxHashMap<String, u64> = FxHashMap::default();
        let mut tr2utr_cov: FxHashMap<String, u64> = FxHashMap::default();
        // the following two storages collections apply to 'clipped' features only
        let mut tr2clipped_cov: FxHashMap<String, u64> = FxHashMap::default();
        // let mut tr2clipped_exons: FxHashMap<String, u16> = FxHashMap::default();

        // for each transcript, infer the intervals to intersect
        // prepare the list of coordinate intervals for intersection
        let mut all_intervals = Vec::new();
        // for processed pseudogene classification, track the 'chain-clipped' statistics
        let mut cds_intervals: Vec<BedEntry> = Vec::new();
        let mut clipped_transcripts: FxHashSet<String> = FxHashSet::default();
        // let mut has_partially_uncovered: bool = false;
        let mut overlaps_cds_exon = false;
        let mut multi_exon_chain: bool = false;
        let mut spans_over_cds = false;

        for transcript in &mut intersected{
            tr_set.insert(transcript.name().unwrap().clone());
            if transcript.thin_start().unwrap() >= chain.refs.start && 
                transcript.thin_end().unwrap() <= chain.refs.end {
                    spans_over_cds = true;
                }
            let tr_name: String = transcript.name().unwrap().clone();
            // let mut has_introns = false;
            // record the total transcript length
            if !transcript_length.contains_key::<str>(&tr_name) {
                transcript_length.insert(
                    tr_name.clone(),
                    transcript.length().unwrap()
                );
            }
            //
            if !exon_number.contains_key::<str>(&tr_name) {
                exon_number.insert(
                    tr_name.clone(),
                    transcript.exon_num().unwrap()
                );
            }
            // add CDS exon blocks to the intersection vector
            let mut exon_counter = 0;
            let mut all_cds_exons: Vec<BedEntry> = transcript
                .to_cds(false)
                .unwrap()// TODO: Currently assumes that all transcripts in the annotation are coding
                .to_blocks()
                .unwrap()
                .into_iter()
                .map(|mut x|
                    {
                        x.update_name(&format!("{}$cds${}", x.name().unwrap(), exon_counter));
                        exon_counter += 1;
                        x
                    }
                )
                .collect();
            if !overlaps_cds_exon {
                let num_exons_covered = all_cds_exons
                    .iter()
                    .map(|x|
                        (intersection(
                            *x.start().unwrap(), *x.end().unwrap(), chain.refs.start, chain.refs.end
                        ).unwrap_or(0) > 0) as u8
                    ).sum::<u8>();
                    overlaps_cds_exon = num_exons_covered > 0;
                    // if !multi_exon_chain {
                    //     multi_exon_chain = num_exons_covered > 1;
                    // }
            }
            if !cds_length.contains_key::<str>(&tr_name) {
                cds_length.insert(
                    tr_name.clone(), 
                    all_cds_exons
                        .iter()
                        .map(|x| x.length().unwrap())
                        .sum()
                );
            }
            cds_intervals.append(&mut all_cds_exons.clone());
            all_intervals.append(&mut all_cds_exons);

            // push the intron intervals
            let all_introns = extract_fraction(
                transcript, 
                BedFractionMode::All, 
                true
            )
                .expect(&format!("Intron extraction failed for transcript {}", tr_name));
            if let Some(x) = all_introns {
                all_intervals.append(
                    &mut x.to_blocks()
                        .unwrap()
                        .into_iter()
                        .map(|mut y|
                            {
                                y.update_name(&format!("{}$intron", x.name().unwrap()));
                                y
                            }
                        )
                        .collect()
                );
                // let clipped_introns = extract_fraction(
                //     transcript, 
                //     BedFractionMode::Cds,
                //     false
                // )
                //     .expect(&format!("Coding sequence intron extraction failed for transcript {}", tr_name))
                //     .unwrap();
                // if let Some(x) = clipped_introns {
                //     all_intervals.append(
                //         &mut x.to_blocks()
                //             .unwrap()
                //             .into_iter()
                //             .map(|mut y|
                //                 {
                //                     y.update_name(&format!("{}$coding_intron", x.name().unwrap()));
                //                     y
                //                 }
                //             )
                //             .collect()
                //     );
                // }
                if !intron_length.contains_key::<str>(&tr_name) {
                    intron_length.insert(tr_name.clone(), x.block_length());
                }
            }

            // extract CDS introns

            // push the flank intervals
            let flank_name = format!("{}$flank", transcript.name().unwrap().clone());
            let left_flank_start = match (*transcript.start().unwrap()).checked_sub(10000) {
                Some(x) => {x},
                None => {0}
            };
            let left_flank_end = *transcript.start().unwrap();//min(transcript.end().unwrap() + 10000, query_chrom_size);
            let left_flank = BedEntry::bed4(
                ref_chrom.clone(), left_flank_start, left_flank_end, flank_name.clone()
            );
            all_intervals.push(left_flank);
            let right_flank_start = *transcript.end().unwrap();
            let right_flank_end = min(right_flank_start + 10000, ref_chrom_size);
            let right_flank = BedEntry::bed4(
                ref_chrom.clone(), right_flank_start,right_flank_end, flank_name.clone()
            );
            all_intervals.push(right_flank);

            // finally, push the UTRs
            // UTR block coverage will be `excluded` from chainwise statistics
            // but accounted for when computing coverage-to-span statistics
            let all_utrs = extract_fraction(
                transcript, 
                BedFractionMode::Utr, 
                false
            )
                .expect(&format!("UTR extraction failed for transcript {}", tr_name));
            if let Some(x) = all_utrs {
                all_intervals.append(
                    &mut x.to_blocks()
                        .unwrap()
                        .into_iter()
                        .map(|mut y|
                            {
                                y.update_name(&format!("{}$utr", x.name().unwrap()));
                                y
                            }
                        )
                        .collect()
                );
            }
            // if coding sequence is only partially covered by the chain, add the 'clipped' intron coverage
            if (transcript.thick_start().unwrap() < chain.refs.start) || (
                transcript.thick_end().unwrap() >= chain.refs.end
            ) {
                // has_partially_uncovered = true;
                clipped_transcripts.insert(tr_name.clone());
                // println!("chain={}, chain start={}, chain end={}", chain_id, chain.refs.start, chain.refs.end);
                let chainclipped_entry = transcript
                    .to_cds(false)
                    .unwrap()
                    .clip_by(Some(chain.refs.start), Some(chain.refs.end), false)
                    .expect(&format!(
                        "Clipping entry {} by chain coordinates ({}-{}) failed:", tr_name, query_chain_start, query_chain_end
                    ));
                let clipped_block_num = chainclipped_entry.exon_num().unwrap();
                if clipped_block_num == 0 {continue}
                let clipped_introns = extract_fraction(
                    &chainclipped_entry, 
                    BedFractionMode::Cds, 
                    true
                ).expect(&format!("Extracting introns for chain-clipped transcript {} at chain {} failed: ", tr_name, chain_id));
                if let Some(x) = clipped_introns {
                    all_intervals.append(
                        &mut x.to_blocks()
                            .unwrap()
                            .into_iter()
                            .map(|mut y|
                                {
                                    y.update_name(&format!("{}$clipped_intron", x.name().unwrap()));
                                    y
                                }
                            )
                            .collect()
                    );
                }
            } else {
                // otherwise, push the coding sequence introns
                // record the coding sequence introns' total length
                let cds_introns = extract_fraction(
                    transcript, 
                    BedFractionMode::Cds,
                    true
                ).expect(&format!("Coding sequence intron extraction failed for transcript {}", tr_name));
                let mut cds_intron_len = 0;
                if let Some(x) = cds_introns {
                    cds_intron_len = x.block_length();
                    all_intervals.append(
                        &mut x.to_blocks()
                            .unwrap()
                            .into_iter()
                            .map(|mut y|
                                {
                                    y.update_name(&format!("{}$clipped_intron", x.name().unwrap()));
                                    y
                                }
                            )
                            .collect()
                    );
                }
                if !cds_intron_length.contains_key::<str>(&tr_name) {
                    cds_intron_length.insert(
                        tr_name.clone(),
                        cds_intron_len
                    );
                }
            }
        }

        // skip chains which do not align at least one coding exon
        // unless the span over the complete cds
        if !overlaps_cds_exon && !spans_over_cds {continue}
        // if *chain_id == 770539 {
        //     println!("overlaps_cds_exon={}, spans_over_cds={}", overlaps_cds_exon, spans_over_cds)
        // }

        // discretize the contents of `all_intervals`
        let all_interval_num: usize = all_intervals.len();
        let (mut discrete_intervals, unit2interval) = discrete_interval_map(&mut all_intervals); 
        if *chain_id <= 10 {
            println!("Chain {} will be intersected to {} discrete intervals  ({} intervals before discretion)", chain_id, discrete_intervals.len(), all_interval_num)
        }
        let (mut discrete_cds_intervals, _) = discrete_interval_map(&mut cds_intervals);
        

        // intersect the discrete intervals to the aligned chain blocks
        let cov_res = chain.alignment_cov_(
            &mut discrete_intervals
        ).expect("Failed exon-to-chain overlap!");

        // extract the clipped coordinates
        let (clipped_coords_ref, clipped_coords_query) = match get_cds_boundaries_(
            chain, &mut discrete_cds_intervals
        ).expect(&format!("Failed to project CDS boundaries via chain {}", chain_id)) {
            (Some(x), Some(y)) => {(x, y)},
            _ => {(Interval::new(), Interval::new())}
        };
        let ref_clipped_start = *clipped_coords_ref.start().unwrap_or(&0);
        let ref_clipped_end = *clipped_coords_ref.end().unwrap_or(&0);
        let non_coding_chain: bool = (ref_clipped_end - ref_clipped_start) == 0;
        let ref_clipped_start = Some(ref_clipped_start);
        let ref_clipped_end = Some(ref_clipped_end);
        
        // if *chain_id <= 10 {
        //     println!("clipped_coords_ref={:#?}, clipped_coords_query={:#?}", clipped_coords_ref, clipped_coords_query);
        // }

        let mut tr2cov_exons: FxHashMap<String, &str> = FxHashMap::default();
        // TODO: Instead of two exons aligned, check for at least two exons spanned in the same transcript
        for (inter_id, intersection) in cov_res.iter() {
            // for (header, intersection) in cov_res.iter() {
            let mut block_cds_coverage: u64 = 0;
            let mut block_utr_coverage: u64 = 0;
            let mut block_intron_coverage: u64 = 0;
            for header in unit2interval.get(*inter_id).unwrap() {
                let tr = header.split('$').nth(0).unwrap().to_string();
                if header.contains("$cds") {
                    // record coverage in a provisional variable
                    block_cds_coverage += *intersection;
                    tr2cds_cov
                        .entry(tr.clone())
                        .and_modify(|x| *x += *intersection)
                        .or_insert(*intersection);
                    if *intersection > 0 {
                        // mark the transcript as counting for the synteny calculation
                        syntenic_transcripts.insert(tr.clone());
                        // count the exon as covered by the chain; if at least two coding exons for the same transcript 
                        // are aligned, the chain is eligible for processed pseudogene classification
                        let ex_num = header.split('$').last().unwrap();
                        if !multi_exon_chain {
                            match tr2cov_exons.get(&*tr) {
                                Some(x) => {
                                    if !(*x == ex_num) {
                                        multi_exon_chain = true
                                    }
                                },
                                None => {
                                    tr2cov_exons.insert(tr.clone(), ex_num);
                                }
                            };
                        }
                    }
                } else if header.contains("$intron") {
                    block_intron_coverage += *intersection;
                    // intron_coverage += *intersection;
                    tr2intron_cov
                        .entry(tr.clone())
                        .and_modify(|x| *x += *intersection)
                        .or_insert(*intersection);
                } else if header.contains("$clipped_intron") {
                    tr2clipped_intron_cov
                        .entry(tr.clone())
                        .and_modify(|x| *x += *intersection)
                        .or_insert(*intersection);
                } else if header.contains("$flank") {
                    tr2flank_cov
                        .entry(tr.clone())
                        .and_modify(|x| *x += *intersection)
                        .or_insert(*intersection);
                } else if header.contains("$utr") {
                    block_utr_coverage += *intersection;
                    tr2utr_cov
                        .entry(tr.clone())
                        .and_modify(|x| *x += *intersection)
                        .or_insert(*intersection);
                } else if header.contains("$clipped_exon") {
                    tr2clipped_cov
                        .entry(tr.clone())
                        .and_modify(|x| *x += *intersection)
                        .or_insert(*intersection);
                    // clipped_cds_coverage += *intersection;
                }  else if header.contains("$clipped_intron") {
                    tr2clipped_cov
                        .entry(tr.clone())
                        .and_modify(|x| *x += *intersection)
                        .or_insert(*intersection);
                    // clipped_intron_coverage += *intersection;
                }
            }
            // if a discrete interval belongs to at least one UTR exon, exclude it from the 
            // global CDS coverage, otherwise add it to the total CDS coverage counter
            if block_utr_coverage > 0 && block_cds_coverage == 0 && block_intron_coverage == 0 {
                utr_coverage += *intersection
            } else {
                if block_cds_coverage > 0 {cds_coverage += *intersection}
                if block_intron_coverage > 0 {intron_coverage += *intersection}
            }
        }

        tr_counter += cov_res.len();

        // a safeguard filter against 'partial' spanning chains
        if syntenic_transcripts.len() == 0 && !spans_over_cds {continue}
        // if syntenic_transcripts.len() == 0 {multi_exon_chain = false}

        // at this point, all the necessary stats must have been obtained
        // first, summarize the remaining chain features
        // first, get synteny value
        let synt = match args.isoform_file {
            Some(_) => {
                let mut syntenic_genes: FxHashSet<&str> = FxHashSet::default();
                // TODO: Currently ignoring the coding sequence requirement
                // to comply with the T1 model architecture
                // for synt_tr in syntenic_transcripts { 
                    for synt_tr in intersected.iter() {
                    // TODO: See above
                    // let gene = tr2gene.get(&synt_tr);
                    let gene = tr2gene.get(synt_tr.name().unwrap());
                    if let Some(x) = gene {
                        syntenic_genes.insert(x);
                    }
                };
                syntenic_genes.len()
            },
            None => {syntenic_transcripts.len()}
        };
        // then get global CDS fraction as C / A
        let global_cds_fraction = cds_coverage as f64 / aln_block_sum as f64;
        // TODO: Corresponds to the T1 code but contradicts the original paper's M&Ms
        // TODO: The original paper doubles down on UTR blocks's exclusion but the code uses bulk cov as denominator
        // let global_cds_fraction = cds_coverage as f64 / (aln_block_sum - utr_coverage) as f64; // ORIG
        // get exon_qlen parameter; for that, divide the sum of CDS-covered based by Q, UTR exons excluded
        let exon_qlen = if multi_exon_chain {cds_coverage as f64 / (query_chain_len - utr_coverage as f64)} else {0.0};
        // TODO: Corresponds to the T1 code but contradicts the original paper's M&Ms
        // TODO: The original paper does not exclude the UTRs yet the code does 
        // let exon_qlen = if multi_exon_chain {cds_coverage as f64 / query_chain_len} else {0.0}; // ORIG

        // for the latter, compute also the 'clipped' analogue; 
        // numerator is 'clipped' CDS + intron coverage, denominator is the CDS-clipped span in the query
        let clipped_exon_qlen = if multi_exon_chain {
            if clipped_coords_query.length().unwrap() == 0 {0.0} else {cds_coverage as f64 / clipped_coords_query.length().unwrap() as f64}
        } else {0.0};
        if *chain_id == 585326 {
            println!("chain_id={}, multi_exon_chain={}", chain_id, multi_exon_chain);
        }

        // then, summarize transcript features and write the resulting line to output
        for transcript in intersected{
            let tr_name: &str = transcript.name().unwrap().as_ref();
            // coding sequence coverage (c); total number of CDS exon bases aligned for this transcript
            let local_cds_coverage = *tr2cds_cov
                .get(tr_name)
                .unwrap_or(&0);
            // intron coverage (i); total number of intron (CDS + UTR) bases aligned for this transcript
            let local_intron_coverage = *tr2intron_cov
                .get(tr_name)
                .unwrap_or(&0);
            // coding sequence intron fraction; needed for processed pseudogdne classification
            let local_clipped_intron_coverage = *tr2clipped_intron_cov
                .get(tr_name)
                .unwrap_or(&0);
            // local aligned fraction (a); total aligned fraction for the focal transcript except for UTR
            let local_aligned_fraction = local_cds_coverage + local_intron_coverage;
            // local CDS fraction; c / a
            let local_cds_fraction = match local_aligned_fraction {
                0 => {0.0},
                _ => {local_cds_coverage as f64 / local_aligned_fraction as f64}
            };
            // retrieve total coding sequence length for the transcript (coding bases only); must be non-zero
            let total_cds_len = *cds_length
                .get(tr_name)
                .expect(&format!("No total coding sequence length recorded for transcript {}", tr_name));
            // same for coding sequence's total intron length; this can be equal to zero for single-exon transcripts
            let total_intron_len = *intron_length
                .get(tr_name)
                .unwrap_or(&0);
                // .expect(&format!("No total intron sequence length recorded for transcript {}", tr_name));

            // get total UTR length for the transcript; can equal zero
            // let utr_cov = *tr2utr_cov
            //     .get(tr_name)
            //     .unwrap_or(&0) as f64;
            
            // retrieve total flank coverage and divide it by total flank length (2 * 10000)
            let flank_cov: f64 = *tr2flank_cov
                .get(transcript.name().unwrap()).unwrap_or(&0) as f64 / 20000.0;

            // retrieve the total transcript length (including UTRs, ~gene length)
            let tr_length = *transcript_length
                .get(tr_name)
                .expect(&format!("No total length recorded for transcript {}", tr_name));
            // and the total exon number
            let exon_num = exon_number
                .get(tr_name)
                .expect(&format!("No exon number recorded for transcript {}", tr_name));

            // last, record the 'clipped' features
            // let clipped_exon_cov: f64;
            let clipped_intron_cov: f64;
            // for transcripts which are only partially covered by the chain, 
            // the clipped intron stat must be recorded with respect to the covered portion's span
            if clipped_transcripts.contains(tr_name) {
                if non_coding_chain {
                    // clipped_exon_cov = 0.0;
                    clipped_intron_cov = 0.0;
                } else {
                    // TODO: Separate coverage tracker for CDS introns' coverage
                    // get the original BedEntry
                    let index = tr2index
                        .get(&*tr_name)
                        .expect(&format!("Missing entry index for transcript {}", tr_name));
                    let orig_tr = chrom2trs
                        .get(&ref_chrom)
                        .unwrap()[*index]
                        .clone() // TODO: Devise an exteranl clip_by() function that does not require input mutability to avoid unnecessary cloning 
                        .clip_by(
                            ref_clipped_start, 
                            ref_clipped_end, 
                            false
                        ).unwrap();
                    if *chain_id == 359705 {
                        println!("ref_clipped_start={}, ref_clipped_end={}, orig_tr={:#?}", ref_clipped_start.unwrap(), ref_clipped_end.unwrap(), orig_tr);
                    }
                    // let clipped_exon_span: f64;
                    let clipped_intron_span: f64;
                    if orig_tr.exon_num().unwrap() == 0 {
                        // clipped_exon_span = 0.0;
                        clipped_intron_span = 0.0;
                    } else {
                        // clipped_exon_span = orig_tr.to_cds(false).unwrap().block_length() as f64;
                        clipped_intron_span = match extract_fraction(&orig_tr, BedFractionMode::Cds, true).unwrap() {
                            Some(x) => {
                                // if *chain_id == 33406 {println!("{:#?}", x);}
                                x.block_length() as f64
                            },
                            None => {0.0}
                        };
                    }
                    // clipped_exon_cov = if clipped_exon_span == 0.0 {0.0} else {local_cds_coverage as f64 / clipped_exon_span};
                    clipped_intron_cov = if clipped_intron_span == 0.0 {-1.0} else {local_clipped_intron_coverage as f64 / clipped_intron_span};
                }
            } else {
                // otherwise, the whole coding sequence intron fraction is considered
                let total_cds_intron = *cds_intron_length.get(&*tr_name).unwrap_or(&0);
                // clipped_exon_cov = local_cds_coverage as f64 / total_cds_len as f64;
                clipped_intron_cov = if total_cds_intron == 0 {0.0} else {local_clipped_intron_coverage as f64 / total_cds_intron as f64}; // TODO: Only CDS intron coverage should be considered
            }

            // file all the features as a single output line
            let out_line = format!(
                "{}\t{}\t{}\t{}\t{}\t{:.16}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                tr_name, tr_num, chain_id, 
                synt, chain_aln_score, global_cds_fraction, 
                ref_chain_len, exon_qlen,
                local_cds_fraction, local_cds_coverage, local_intron_coverage,
                tr_length, exon_num,
                total_cds_len, total_intron_len, flank_cov,
                clipped_exon_qlen, clipped_intron_cov //clipped_exon_cov, 
                // clipped_cds, clipped_utr
            );
            if let Err(e) = output_file.write(out_line.as_bytes()) {
                eprintln!("Failed to write the line: {}", e);
            }
        }

        i += 1;
    }
    let elapsed = now.elapsed();
    println!("Iterating over all chains; elapsed time: {:?}", elapsed);
    println!("Processed {} chains, with a total of {} non-redundant transcripts and {} redundant intervals processed", i, tr_set.len(), tr_counter);

    println!("Total time elapsed: {:?}", program_start.elapsed());
}


// ## LOGICAL 'OR' FOR UTRs