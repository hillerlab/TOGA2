// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! CLI argument parsing and mode dispatch.
//!
//! Defines the [`Args`] struct parsed from command-line flags, the [`Mode`]
//! enum selecting input layout, and the [`Query`] enum selecting the scope of
//! assemblies to process.

use clap::{self, ArgGroup, Parser};
use log;

use std::path::PathBuf;

/// Command-line arguments for `toga2_stats`.
///
/// Two mutually exclusive modes control the input source:
///
/// * `--hillerlab` — hard-coded base paths under `/projects/hillerlab/genome/gbdb-HL`
/// * `--togadir` — direct path to a TOGA2 output directory
///
/// Within HillerLab mode, exactly one of `--query`, `--queries`, or `--all`
/// must be provided.
///
/// # Fields
///
/// * `togadir` — Path to a TOGA2 output directory (Public mode)
/// * `reference` — Reference assembly name (required in HillerLab mode)
/// * `query` — Single query assembly name
/// * `queries` — File listing query assembly names, one per line
/// * `hillerlab` — Enable HillerLab-specific base paths
/// * `no_header` — Suppress the TSV header line in output
/// * `preserve_order` — Preserve input order of query assemblies in output
/// * `all` — Process every `vs_*` subdirectory under the reference
/// * `ancestral` — Path to ancestral gene set file
/// * `level` — Logging verbosity (default: Info)
/// * `threads` — Worker thread count (default: number of CPU cores)
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(group(
    ArgGroup::new("mode")
        .args(&["hillerlab", "togadir"])
        .required(true)
        .multiple(false)
))]
#[command(group(
    ArgGroup::new("q_args")
        .args(&["query", "queries", "all"])
        .required(false)
        .multiple(false)
))]
pub struct Args {
    #[arg(
        short = 't',
        long = "togadir",
        help = "Path to TOGA2 output directory",
        value_name = "PATH",
        required = false
    )]
    pub togadir: Option<PathBuf>,

    #[arg(
        short = 'r',
        long = "ref",
        help = "Reference assembly name to process",
        value_name = "REFERENCE",
        requires_if("true", "hillerlab")
    )]
    pub reference: Option<String>,

    #[arg(
        short = 'q',
        long = "query",
        help = "Query assembly name to process (assumes TOGA2 output)",
        value_name = "QUERY",
        requires_if("true", "hillerlab"),
        conflicts_with = "all",
        conflicts_with = "queries"
    )]
    pub query: Option<String>,

    #[arg(
        short = 'Q',
        long = "queries",
        help = "Query assemblies to read from a file (assumes TOGA2 output)",
        value_name = "FILE",
        requires_if("true", "hillerlab"),
        conflicts_with = "all",
        conflicts_with = "query"
    )]
    pub queries: Option<String>,

    #[arg(
        short = 'H',
        long = "hillerlab",
        help = "Activate Hillerlab specific paths",
        value_name = "FLAG",
        action = clap::ArgAction::SetTrue,
    )]
    pub hillerlab: bool,

    #[arg(
        short = 'N',
        long = "no-header",
        help = "Do not print header",
        value_name = "FLAG",
        action = clap::ArgAction::SetTrue,
    )]
    pub no_header: bool,

    #[arg(
        short = 'P',
        long = "preserve-order",
        help = "Preserve order of query assemblies",
        value_name = "FLAG",
        action = clap::ArgAction::SetTrue,
    )]
    pub preserve_order: bool,

    #[arg(
        short = 'A',
        long = "all",
        help = "Process all samples in the input directory",
        requires_if("true", "hillerlab"),
        value_name = "FLAG",
        conflicts_with = "query",
        conflicts_with = "queries",
        action = clap::ArgAction::SetTrue,
    )]
    pub all: bool,

    #[arg(
        short = 'a',
        long = "ancestral",
        help = "Path to ANCESTRAL file directory",
        value_name = "PATH",
        required = false
    )]
    pub ancestral: Option<PathBuf>,

    #[arg(
        short = 'L',
        long = "level",
        help = "Logging level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: log::Level,

    #[arg(
        short = 'T',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,
}

/// Operating mode selection.
///
/// Determines whether the tool uses HillerLab's hard-coded directory layout
/// or a user-supplied TOGA2 output directory.
#[derive(Clone, Debug)]
pub enum Mode {
    /// Use hard-coded paths under `/projects/hillerlab/genome/gbdb-HL`.
    HillerLab,
    /// Use a user-supplied TOGA2 output directory.
    Public,
}

impl From<bool> for Mode {
    /// Converts a boolean to a [`Mode`].
    ///
    /// `true` maps to [`Mode::HillerLab`], `false` to [`Mode::Public`].
    fn from(value: bool) -> Self {
        if value { Mode::HillerLab } else { Mode::Public }
    }
}

impl Args {
    /// Determines which [`Query`] variant applies based on the parsed flags.
    ///
    /// Checks `--all`, `--query`, and `--queries` in priority order and
    /// returns the corresponding variant, or [`Query::None`] if none is set.
    pub fn get_query_type(&self) -> Query {
        if self.all {
            Query::All
        } else if self.query.is_some() {
            Query::Single
        } else if self.queries.is_some() {
            Query::Group
        } else {
            Query::None
        }
    }
}

/// Query-assembly selection scope.
///
/// Controls which assemblies are processed within HillerLab mode.
#[derive(Debug, Clone, Copy)]
pub enum Query {
    /// Process every `vs_*` subdirectory under the reference.
    All,
    /// Process a single query assembly.
    Single,
    /// Process a list of query assemblies from a file.
    Group,
    /// No query selection (placeholder).
    None,
}
