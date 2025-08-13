use clap::Parser;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::io::prelude::*;
use std::path::Path;

use cubiculum::extract::extract::{bed_to_fraction, BedFractionMode};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// A combination of HillerLab bed12ToCDSOnly and UCSC bed12ToBed6 utilities
struct Args {
    /// Input file; if set to 'stdin', expects data to come from standard output stream
    #[arg(long, short = 'i', default_value_t = String::from("stdin"))]
    input: String,

    /// Output file; if set to 'stdout', will write the output data to standard output stream
    #[arg(long, short = 'o', default_value_t = String::from("stdout"))]
    output: String,

    /// Fraction to report; possible values are:
    /// cds - for coding sequence;
    /// utr - for all untranslated sequence;
    /// 5utr - for 5'-side untranslated region;
    /// 3utr - for 3'-side untranslated region
    #[arg(long, short='m', default_value_t = String::from("all"))]
    mode: String,

    /// If set, intron intervals will be reported instead of exons as BED12 blocks
    #[arg(long, short='n', action)]
    intron: bool,

    /// If set, output format is switched to BED6; each interval will be reported 
    /// as a separate BED6 entry
    #[arg(long, short='b', action)]
    bed6: bool

}

fn main() {
    let args = Args::parse();

    match args.mode.as_str() {
        "all" => { BedFractionMode::All },
        "cds" => { BedFractionMode::Cds },
        "utr" => { BedFractionMode::Utr },
        "5utr" => { BedFractionMode::Utr5 },
        "3utr" => { BedFractionMode::Utr3 },
        _ => {
            panic!("Invalid 'mode' has been provided: {}. Valid modes are: all, cds, utr, 3utr, 5utr", args.mode)
        }
    };

    // let mut output_file = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     // .append(true)
    //     .open(args.output)
    //     .unwrap();

    let input_file = match args.input.as_str() {
        "stdin" => {Box::new(io::stdin().lock()) as Box<dyn BufRead>},
        _ => {
            let path = File::open(args.input).unwrap();
            Box::new(BufReader::new(path)) as Box<dyn BufRead>
        }
    };

    let mut output_file = match args.output.as_str() {
        "stdout" => {Box::new(io::stdout()) as Box<dyn Write>},
        _ => {
            let path = Path::new(&args.output);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };

    for line_ in input_file.lines() {
        if let Ok(line) = line_ {
                let result: Option<String> = bed_to_fraction(line, &args.mode, args.intron, args.bed6);
                if let Some(fraction) = result {
                // println!("{}", fraction);
                if let Err(e) = writeln!(output_file, "{}", fraction) {
                    eprintln!("Failed to write the line: {}", e);
                }
            }
        }
    }

}
