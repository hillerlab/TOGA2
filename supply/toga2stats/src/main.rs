// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::Parser;
use dashmap::DashMap;
use rayon::prelude::*;
use simple_logger::init_with_level;
use toga2_stats::{
    cli::{Args, Mode, Query},
    utils::reader,
};

use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
    vec,
};

/// Root of the HillerLab genome database tree.
const GBDB_HL: &str = "/projects/hillerlab/genome/gbdb-HL";

/// Name of the loss-summary file inside each TOGA2 output directory.
const LOSS_SUMMARY: &str = "loss_summary.tsv";

/// Name of the TOGA2 subdirectory under each reference genome.
const TOGA2: &str = "TOGA2";

/// Default ancestral gene set file when `--ancestral` is not provided.
const ANCESTRAL: &str =
    "/projects/hillerlab/genome/gbdb-HL/hg38/TOGA2/AncestralGeneSet/Ancestral_placental.txt";

/// TSV header line printed when `--no-header` is not set.
const HEADER: &str = "ASSEMBLY\tGENE_INTACT\tGENE_INACT\tGENE_MISSING\tTRANSCRIPT_INTACT\tTRANSCRIPT_INACT\tTRANSCRIPT_MISSING\ttANCESTRAL_INTACT\tANCESTRAL_INACT\tANCESTRAL_MISSING";

/// Entry point: parses CLI arguments, selects mode, processes assemblies,
/// and prints a TSV summary line (or lines) to stdout.
fn main() {
    let start = std::time::Instant::now();
    let args = Args::parse();

    init_with_level(args.level).unwrap();

    let mode = Mode::from(args.hillerlab);

    let accumulator = ParallelAccumulator::default();
    let ancestral = if let Some(ancestral) = &args.ancestral {
        get_ancestral(ancestral)
    } else {
        get_ancestral(&PathBuf::from(ANCESTRAL))
    };

    if !args.no_header {
        println!("{HEADER}");
    }

    match mode {
        Mode::HillerLab => {
            log::debug!("INFO: Running in HillerLab mode!");
            let query_type = Args::get_query_type(&args);

            match query_type {
                Query::All | Query::Group => {
                    log::debug!("INFO: Mode choosen -> All/Group queries");
                    __distribute(&args, &ancestral, query_type);
                }
                Query::Single => {
                    log::debug!("INFO: Mode choosen -> Single query");

                    let reference = &args.reference.unwrap_or_else(|| {
                        log::error!("ERROR: --ref is required in HillerLab mode!");
                        std::process::exit(1);
                    });
                    let query = args.query.unwrap_or_else(|| {
                        log::error!("ERROR: --query is required in HillerLab mode!");
                        std::process::exit(1);
                    });

                    let assembly = PathBuf::from(GBDB_HL)
                        .join(reference)
                        .join(TOGA2)
                        .join(format!("vs_{}", &query));

                    log::debug!("INFO: Target assembly -> {assembly:?}");
                    __get_projections(&assembly, &accumulator, &ancestral);

                    let line = accumulator
                        .to_line(ancestral.len())
                        .replace("ASSEMBLY", &query);
                    println!("{}", line);
                }
                Query::None => unreachable!(),
            }
        }
        Mode::Public => {
            log::debug!("INFO: Running TOGA2 directory mode!");

            let target = &args.togadir.unwrap_or_else(|| {
                log::error!("ERROR: --togadir is required in TOGA2 directory mode!");
                std::process::exit(1);
            });

            __get_projections(target, &accumulator, &ancestral);

            let line = accumulator
                .to_line(ancestral.len())
                .replace("ASSEMBLY", target.file_stem().unwrap().to_str().unwrap());

            println!("{}", line);
        }
    }

    let elapsed = start.elapsed();
    log::debug!("Elapsed time: {:.3?}", elapsed);
}

/// Processes multiple query assemblies in HillerLab mode.
///
/// Handles both [`Query::All`] (scans the TOGA2 directory for `vs_*` subdirs)
/// and [`Query::Group`] (reads a file listing query names). Each assembly is
/// processed in parallel using rayon. Optionally preserves input order via an
/// `RwLock<HashMap>` when `--preserve-order` is set.
fn __distribute(args: &Args, ancestral: &HashSet<String>, query_type: Query) {
    let reference = args.reference.as_ref().unwrap_or_else(|| {
        log::error!("ERROR: --ref is required in HillerLab mode!");
        std::process::exit(1);
    });
    let origin = PathBuf::from(GBDB_HL).join(reference).join(TOGA2);

    if !origin.exists() {
        log::error!("ERROR: Reference path does not exist: {origin:?}");
        std::process::exit(1);
    }

    match query_type {
        Query::All => {
            let entries = std::fs::read_dir(&origin).unwrap_or_else(|e| {
                log::error!("ERROR: Could not read TOGA2 directory: {e}");
                std::process::exit(1);
            });

            let assemblies: Vec<PathBuf> = entries
                .filter_map(|entry| entry.ok())
                .map(|entry| entry.path())
                .filter(|path| {
                    path.is_dir()
                        && path
                            .file_name()
                            .unwrap()
                            .to_str()
                            .unwrap()
                            .starts_with("vs_")
                })
                .collect();

            let assembly_map: Option<std::sync::RwLock<HashMap<PathBuf, String>>> =
                if args.preserve_order {
                    let mut mapper = HashMap::new();

                    for assembly in assemblies.iter() {
                        mapper.insert(assembly.clone(), String::new());
                    }
                    Some(std::sync::RwLock::new(mapper))
                } else {
                    None
                };

            log::debug!("INFO: Found {} assemblies", assemblies.len());

            assemblies.par_iter().for_each(|assembly| {
                let local_accumulator = ParallelAccumulator::default();

                log::debug!(
                    "DEBUG: Processing assembly {:?} in {:?}",
                    assembly.file_name().unwrap(),
                    assembly
                );
                __get_projections(assembly, &local_accumulator, ancestral);

                let line = local_accumulator
                    .to_line(ancestral.len())
                    .replace("ASSEMBLY", assembly.file_name().unwrap().to_str().unwrap());

                if !args.preserve_order {
                    println!("{}", line);
                } else {
                    let mut handle = assembly_map.as_ref().unwrap().write().unwrap_or_else(|e| {
                        log::error!("ERROR: Could not lock assembly map: {e}");
                        std::process::exit(1);
                    });
                    if let Some(name) = handle.get_mut(assembly) {
                        *name = line;
                    }
                }
            });

            if args.preserve_order {
                let mapper = assembly_map
                    .unwrap_or_else(|| panic!("ERROR: Could not create assembly map!"))
                    .into_inner()
                    .unwrap_or_else(|e| panic!("ERROR: Could not get assembly map! -> {e}"));

                for assembly in assemblies {
                    let line = mapper.get(&assembly).unwrap_or_else(|| {
                        panic!("ERROR: Could not get assembly map! -> {assembly:?}")
                    });

                    println!("{}", line);
                }
            }
        }
        Query::Group => {
            // INFO: need to read file and get list of queries
            let queries_file = args.queries.as_ref().unwrap_or_else(|| {
                log::error!(
                    "ERROR: --queries is required in HillerLab mode when --all is not set. Seems that there is nothing here!"
                );
                std::process::exit(1);
            });

            let queries = reader(queries_file).unwrap_or_else(|e| {
                log::error!("ERROR: Could not read queries file: {e}");
                std::process::exit(1);
            });

            let query_list: Vec<PathBuf> = queries
                .lines()
                .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
                .map(|line| origin.join(format!("vs_{}", line.trim())))
                .collect();

            log::debug!("INFO: Found {} queries", query_list.len());

            let assembly_map: Option<std::sync::RwLock<HashMap<PathBuf, String>>> =
                if args.preserve_order {
                    let mut mapper = HashMap::new();

                    for assembly in query_list.iter() {
                        mapper.insert(assembly.clone(), String::new());
                    }
                    Some(std::sync::RwLock::new(mapper))
                } else {
                    None
                };

            query_list.par_iter().for_each(|assembly| {
                let local_accumulator = ParallelAccumulator::default();

                if !assembly.exists() {
                    log::warn!("WARN: Assembly path does not exist: {assembly:?}, skipping...");
                    return;
                }

                log::debug!(
                    "DEBUG: Processing assembly {:?} in {:?}",
                    assembly.file_name().unwrap(),
                    assembly
                );
                __get_projections(assembly, &local_accumulator, ancestral);

                let line = local_accumulator
                    .to_line(ancestral.len())
                    .replace("ASSEMBLY", assembly.file_name().unwrap().to_str().unwrap());

                if !args.preserve_order {
                    println!("{}", line);
                } else {
                    let mut handle = assembly_map.as_ref().unwrap().write().unwrap_or_else(|e| {
                        log::error!("ERROR: Could not lock assembly map: {e}");
                        std::process::exit(1);
                    });
                    if let Some(name) = handle.get_mut(assembly) {
                        *name = line;
                    }
                }
            });

            if args.preserve_order {
                let mapper = assembly_map
                    .unwrap_or_else(|| panic!("ERROR: Could not create assembly map!"))
                    .into_inner()
                    .unwrap_or_else(|e| panic!("ERROR: Could not get assembly map! -> {e}"));

                for assembly in query_list {
                    let line = mapper.get(&assembly).unwrap_or_else(|| {
                        panic!("ERROR: Could not get assembly map! -> {assembly:?}")
                    });

                    println!("{}", line);
                }
            }
        }
        _ => unreachable!(),
    }
}

/// Thread-safe counters for gene, transcript, and ancestral gene statuses.
///
/// Backed by `DashMap` for lock-free concurrent access from rayon-parallel
/// iterations. Each map uses string keys (`"intact"`, `"inact_mut"`,
/// `"missing_seq"`) to categorise counts.
///
/// # Fields
///
/// * `genes` — Counts for gene-level projections
/// * `transcripts` — Counts for transcript-level projections
/// * `ancestral` — Counts restricted to genes present in the ancestral set
pub struct ParallelAccumulator {
    pub genes: DashMap<String, usize>,
    pub transcripts: DashMap<String, usize>,
    pub ancestral: DashMap<String, usize>,
}

impl Default for ParallelAccumulator {
    /// Creates an empty [`ParallelAccumulator`] with all DashMaps initialised.
    fn default() -> Self {
        ParallelAccumulator {
            genes: DashMap::new(),
            transcripts: DashMap::new(),
            ancestral: DashMap::new(),
        }
    }
}

impl ParallelAccumulator {
    /// Formats the accumulated counts into a single TSV line.
    ///
    /// The line contains 10 tab-separated fields:
    /// `ASSEMBLY`, `GENE_INTACT`, `GENE_INACT`, `GENE_MISSING`,
    /// `TRANSCRIPT_INTACT`, `TRANSCRIPT_INACT`, `TRANSCRIPT_MISSING`,
    /// `ANCESTRAL_INTACT`, `ANCESTRAL_INACT`, `ANCESTRAL_MISSING`.
    ///
    /// The placeholder `"ASSEMBLY"` is replaced by the caller with the actual
    /// assembly name. `ANCESTRAL_MISSING` is computed as
    /// `total_ancestral - (intact + inact_mut)`.
    ///
    /// # Arguments
    ///
    /// * `ancestral_len` — Total number of genes in the ancestral set
    pub fn to_line(&self, ancestral_len: usize) -> String {
        let mut line = String::from("ASSEMBLY\t");

        for category in vec![&self.genes, &self.transcripts] {
            let intact = category.get("intact").map(|v| *v).unwrap_or(0);
            let missing = category.get("missing_seq").map(|v| *v).unwrap_or(0);
            let inact = category.get("inact_mut").map(|v| *v).unwrap_or(0);

            line.push_str(&format!("{intact}\t{inact}\t{missing}\t",));
        }

        let intact = &mut self.ancestral.get("intact").map(|v| *v).unwrap_or(0);
        let missing = &mut self.ancestral.get("missing_seq").map(|v| *v).unwrap_or(0);
        let inact = &mut self.ancestral.get("inact_mut").map(|v| *v).unwrap_or(0);

        let absent = ancestral_len - (*intact + *missing + *inact);
        *missing += absent;

        line.push_str(&format!("{intact}\t{inact}\t{missing}\t",));

        line
    }
}

/// Reads the ancestral gene set file into a `HashSet<String>`.
///
/// Lines starting with `#` are skipped. The first whitespace-delimited token
/// of each remaining line is added to the set. Uses rayon for parallel
/// fold–reduce.
///
/// # Arguments
///
/// * `ancestral` — Path to the ancestral gene set file
fn get_ancestral(ancestral: &PathBuf) -> HashSet<String> {
    let ancestral = reader(ancestral).unwrap_or_else(|e| {
        log::error!("ERROR: Could not read projections file: {e}");
        std::process::exit(1);
    });

    let tracks = ancestral
        .par_lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.split_whitespace().next())
        .fold(
            || HashSet::new(),
            |mut acc: HashSet<String>, record| {
                acc.insert(record.to_string());
                acc
            },
        )
        .reduce(
            || HashSet::new(),
            |mut acc, map| {
                acc.extend(map);
                acc
            },
        );

    log::debug!("INFO: Read {} ancestral genes", tracks.len());
    tracks
}

/// Parses the `loss_summary.tsv` inside an assembly directory and updates the
/// accumulator.
///
/// Skips comment lines (starting with `#`) and the header line (starting with
/// `"level"`). Each remaining line is parsed via [`Projection::parse`].
///
/// # Arguments
///
/// * `assembly` — Path to the assembly directory containing `loss_summary.tsv`
/// * `accumulator` — Target accumulator for parsed counts
/// * `ancestral` — Set of ancestral gene names for ancestral tracking
fn __get_projections(
    assembly: &PathBuf,
    accumulator: &ParallelAccumulator,
    ancestral: &HashSet<String>,
) {
    if assembly.join(LOSS_SUMMARY).exists() {
        let projections = reader(assembly.join(LOSS_SUMMARY)).unwrap_or_else(|e| {
            log::error!("ERROR: Could not read projections file: {e}");
            std::process::exit(1);
        });

        projections
            .par_lines()
            .filter(|line| !line.starts_with('#') && !line.starts_with("level"))
            .for_each(|line| {
                if let Err(e) = Projection::parse(line, accumulator, ancestral) {
                    log::trace!("WARN: Error parsing {}: {}", line, e);
                }
            });
    } else {
        log::debug!("ERROR: loss_summary does not exist for TOGA2 dir -> {assembly:?}")
    }
}

/// A single line from a `loss_summary.tsv` file.
///
/// Records the name, hierarchical level, and projection status of a gene or
/// transcript. The optional `class` field can be set after parsing for
/// downstream filtering.
///
/// # Fields
///
/// * `name` — Gene or transcript identifier
/// * `level` — Whether this is a gene, transcript, or projection entry
/// * `status` — Projection status (intact, lost, missing, etc.)
/// * `class` — Optional classification label
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Projection {
    pub name: String,
    pub level: ProjectionLevel,
    pub status: Status,
    pub class: Option<String>,
}

/// Hierarchical level of a projection entry.
///
/// Each line in `loss_summary.tsv` belongs to one of these levels.
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub enum ProjectionLevel {
    /// Gene-level projection.
    Gene,
    /// Transcript-level projection.
    Transcript,
    /// Raw projection (skipped during processing).
    Projection,
}

impl From<&str> for ProjectionLevel {
    /// Parses a `ProjectionLevel` from its string representation.
    ///
    /// `"GENE"`, `"TRANSCRIPT"`, and `"PROJECTION"` are valid inputs. Any
    /// other value logs an error and exits.
    fn from(s: &str) -> Self {
        match s {
            "GENE" => ProjectionLevel::Gene,
            "TRANSCRIPT" => ProjectionLevel::Transcript,
            "PROJECTION" => ProjectionLevel::Projection,
            _ => {
                log::error!("ERROR: Invalid projection level: {}", s);
                std::process::exit(1);
            }
        }
    }
}

impl Projection {
    /// Parses a tab-separated line from `loss_summary.tsv`.
    ///
    /// Fields (tab-separated):
    ///   1. Level (`GENE`, `TRANSCRIPT`, or `PROJECTION`)
    ///   2. Name (gene or transcript identifier)
    ///   3. Status code (`FI`, `I`, `PI`, `M`, `PM`, `PG`, `L`, `UL`, `N`)
    ///
    /// `PROJECTION`-level lines are skipped. For `GENE` lines, the function
    /// also updates the ancestral counter when the gene name appears in the
    /// ancestral set.
    ///
    /// # Arguments
    ///
    /// * `line` — A single line from `loss_summary.tsv`
    /// * `accumulator` — Target accumulator for parsed counts
    /// * `ancestral` — Set of ancestral gene names
    pub fn parse(
        line: &str,
        accumulator: &ParallelAccumulator,
        ancestral: &HashSet<String>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut data = line.split('\t');

        let level = ProjectionLevel::from(data.next().unwrap_or_else(|| {
            log::error!("ERROR: Invalid line, missing level field: {}", line);
            std::process::exit(1);
        }));

        match level {
            ProjectionLevel::Projection => return Err("Skipping PROJECTION level".into()),
            _ => {}
        }

        let projection = data.next().unwrap_or_else(|| {
            log::error!("ERROR: Invalid line, missing name field: {}", line);
            std::process::exit(1);
        });

        let status = Status::from(data.next().unwrap_or_else(|| {
            log::error!("ERROR: Invalid line, missing status field: {}", line);
            std::process::exit(1);
        }));

        match level {
            ProjectionLevel::Gene => match status {
                Status::FullyIntact | Status::Intact => {
                    accumulator
                        .genes
                        .entry("intact".to_string())
                        .and_modify(|e| *e += 1)
                        .or_insert(1);

                    if ancestral.contains(projection) {
                        accumulator
                            .ancestral
                            .entry("intact".to_string())
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                }
                Status::PartiallyIntact
                | Status::Missing
                | Status::PartiallyMissing
                | Status::Pseudogene
                | Status::Absent => {
                    accumulator
                        .genes
                        .entry("missing_seq".to_string())
                        .and_modify(|e| *e += 1)
                        .or_insert(1);

                    if ancestral.contains(projection) {
                        accumulator
                            .ancestral
                            .entry("missing_seq".to_string())
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                }
                Status::Lost | Status::UncertainLoss => {
                    accumulator
                        .genes
                        .entry("inact_mut".to_string())
                        .and_modify(|e| *e += 1)
                        .or_insert(1);

                    if ancestral.contains(projection) {
                        accumulator
                            .ancestral
                            .entry("inact_mut".to_string())
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    }
                }
            },
            ProjectionLevel::Transcript => match status {
                Status::FullyIntact | Status::Intact => {
                    accumulator
                        .transcripts
                        .entry("intact".to_string())
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                }
                Status::PartiallyIntact
                | Status::Missing
                | Status::PartiallyMissing
                | Status::Pseudogene
                | Status::Absent => {
                    accumulator
                        .transcripts
                        .entry("missing_seq".to_string())
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                }
                Status::Lost | Status::UncertainLoss => {
                    accumulator
                        .transcripts
                        .entry("inact_mut".to_string())
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                }
            },
            _ => unreachable!(),
        }

        Ok(())
    }

    /// Assigns a classification label to this projection.
    ///
    /// # Arguments
    ///
    /// * `class` — Classification string
    pub fn set_class(&mut self, class: String) {
        self.class = Some(class);
    }
}

/// Projection status of a gene or transcript.
///
/// Mapped from two-letter codes found in `loss_summary.tsv`.
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub enum Status {
    /// Fully intact projection (`FI`).
    FullyIntact,
    /// Intact projection (`I`).
    Intact,
    /// Partially intact projection (`PI`).
    PartiallyIntact,
    /// Missing sequence (`M`).
    Missing,
    /// Partially missing sequence (`PM`).
    PartiallyMissing,
    /// Pseudogene (`PG`).
    Pseudogene,
    /// Absent / not found (`N` or any unknown code).
    Absent,
    /// Lost projection (`L`).
    Lost,
    /// Uncertain loss (`UL`).
    UncertainLoss,
}

impl From<&str> for Status {
    /// Parses a [`Status`] from its two-letter code.
    ///
    /// Known codes: `FI`, `I`, `PI`, `M`, `PM`, `PG`, `L`, `UL`. Any other
    /// value (including `N`) maps to [`Status::Absent`].
    fn from(s: &str) -> Self {
        match s {
            "FI" => Status::FullyIntact,
            "I" => Status::Intact,
            "PI" => Status::PartiallyIntact,
            "M" => Status::Missing,
            "PM" => Status::PartiallyMissing,
            "PG" => Status::Pseudogene,
            "L" => Status::Lost,
            "UL" => Status::UncertainLoss,
            "N" | _ => Status::Absent,
        }
    }
}
