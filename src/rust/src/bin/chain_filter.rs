use clap::Parser;
use flate2::read::MultiGzDecoder;
use memchr::memchr;
use core::str;
use std::path::Path;
use std::fs::File;
use std::io::{self, Read, Write};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Given a chain file and a minimal score value X, picks the chains with score of at least X 
/// and writes them to the output destination
struct Args {
    /// Input chain file; can be compressed into gzip format as long as the file has the ".gz" extension
    #[arg(long, short = 'i')]
    chain_file: String,

    /// Output file; if not value is provided or value is set to "stdout", 
    /// the results are written to standard output
    #[arg(long, short = 'o', default_value_t = String::from("stdout"))]
    output: String,

    /// Minimal chain score; chain scoring below this value are filtered out from output
    #[arg(long, short='s', default_value_t = 0)]
    min_score: u64
}

fn main() {
    let args = Args::parse();
    let start_time = Instant::now();
    // initialize the output file
    let mut output_file = match args.output.as_str() {
        "stdout" => {Box::new(io::stdout()) as Box<dyn Write>},
        _ => {
            let path = Path::new(&args.output);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };

    // open the file as suggested in chaintools::Reader::open()
    let stat = std::fs::metadata(&args.chain_file)
        .expect(&format!("Failed to get metadata for {:?}", args.chain_file));
    let mut data = vec![];
    data.reserve(stat.len() as usize + 1);

    let mut f = File::open(&args.chain_file)
        .expect(&format!("Failed to open file {:?}", args.chain_file));

    let input_path = Path::new(&args.chain_file);
    if input_path.extension().unwrap() == "gz" {
        let mut decoder = MultiGzDecoder::new(f);
        decoder
            .read_to_end(&mut data)
            .expect( &format!("Failed to read file {:?}", args.chain_file));
    } else {
        f.read_to_end(&mut data)
            .expect( &format!("Failed to read file {:?}", args.chain_file));
    }

    // direct ripoff from the chaintools from_file() function
    let mut data = &data[..];
    loop {
        let sep = memchr(b'\n', &data).expect(
            &format!(
                "No newline symbols found"
            )
        );
        // println!("Block delimitation time: {:?}", block_delim_time.elapsed());
        let Some(end) = memchr(b'c', &data[sep..]) else {
            let header = &data[..sep];
            // check the score value
            let score_start = memchr(b' ', header).expect(
                "Improperly formatted header"
            );
            let score_end = memchr(b' ', &header[score_start+1..]).expect(
                "Improperly formatted header"
            );
            let score_field = str::from_utf8(&header[score_start+1..score_end+score_start+1])
                .unwrap()
                .parse::<u64>();
            let score = match score_field {
                Ok(x) => {x},
                Err(_) => {
                    if header[score_start+1] == b'-' {
                        eprintln!("Negative score value encountered; skipping");
                        break
                    } else {
                        panic!("Failed to parse the score line: {:?}", String::from_utf8_lossy(header))
                    }
                }
            };
            if score >= args.min_score {
                // write the contents to the file
                if let Err(e) = output_file.write(data) {
                    eprintln!("Failed to write the chain data to output: {:?}", e);
                };
            }
            break;
        };
        let header = &data[..sep];
        // check the score value
        let score_start = memchr(b' ', header).expect(
            "Improperly formatted header"
        );
        let score_end = memchr(b' ', &header[score_start+1..]).expect(
            "Improperly formatted header"
        );
        let score_field = str::from_utf8(&header[score_start+1..score_end+score_start+1])
            .unwrap()
            .parse::<u64>();
        let score = match score_field {
            Ok(x) => {x},
            Err(_) => {
                if header[score_start+1] == b'-' {
                    eprintln!("Negative score value encountered; skipping");
                    data = &data[sep + end..];
                    continue
                } else {
                    panic!("Failed to parse the score line: {:?}", String::from_utf8_lossy(header))
                }
            }
        };
        if score < args.min_score {
            data = &data[sep + end..];
            continue
        }
        // write the valid chain instance
        if let Err(e) = output_file.write(&data[..sep + end]) {
            eprintln!("Failed to write the chain data to output: {:?}", e);
        };

        data = &data[sep + end..];
    }
    println!("Chain filtering time: {:?}", start_time.elapsed());
}
