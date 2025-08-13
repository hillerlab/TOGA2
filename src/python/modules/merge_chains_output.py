#!/usr/bin/env python3
"""Parse raw chain runner output.

Chain features extraction steps results in numerous files.
This script merges these files and then
builds s table containing chain features.
"""
import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

import argparse
from dataclasses import dataclass
from datetime import datetime as dt
from collections import defaultdict
from modules.shared import CommandLineManager, CONTEXT_SETTINGS
from typing import Dict, List, Optional, TextIO, Union
# from version import __version__

import click

__author__ = ''
__credits__ = ("Bogdan M. Kirilenko")

GENES: str = 'genes'
CHAIN: str = 'chain'
TIMESTAMP: str = '#estimated'
HEADER: str = '\t'.join(
    (
        'transcript', 'gene_overs', 'chain', 'synt', 'gl_score', 'gl_exo', 'chain_len', 
        'exon_qlen', 'loc_exo', 'exon_cover', 'intr_cover', 'gene_len', 'ex_num', 
        'ex_fract', 'intr_fract', 'flank_cov', 'clipped_exon_qlen', 'clipped_intr_cover'
    )
)
OK: str = 'ok'
# t0 = dt.now()

@dataclass
class TranscriptFeatures:
    __slots__ = ('gene_len', 'exon_fraction', 'intron_fraction', 'exons_num')
    gene_len: int
    exon_fraction: float
    intron_fraction: float
    exons_num: int

@dataclass
class ChainFeatures:
    __slots__ = (
        'synteny', 'global_score', 
        'global_exons', 'chain_q_len', 'chain_len',
        'local_exons', 'coverages', 'introns', 'flank_cov',
        'clipped_chain_qlen', 'clipped_exons', 'clipped_introns'
    )
    synteny: str
    global_score: str
    global_exons: str
    chain_q_len: str
    chain_len: str
    local_exons: Dict[str, str]
    coverages: Dict[str, str]
    introns: Dict[str, str]
    flank_cov: Dict[str, str]
    clipped_chain_qlen: str
    clipped_exons: Dict[str, str]
    clipped_introns: Dict[str, str]


# def parse_args():
#     """Read args, check."""
#     app = argparse.ArgumentParser()
#     app.add_argument("results_dir", type=str, help="Directory containing the results.")
#     app.add_argument(
#         "bed_file", type=str, help="Bed file containing annotations for genes analyzed."
#     )
#     app.add_argument("output", type=str, help="Save output here.")
#     app.add_argument(
#         "--isoforms",
#         "-i",
#         type=str,
#         default=None,
#         help="File with isoforms. Each line means "
#         "different ids of a gene, like gene_id<tab>[,-separated list of transcripts].",
#     )
#     app.add_argument(
#         "--verbose", "-v", action="store_true", dest="verbose", help="Verbose messages."
#     )
#     app.add_argument(
#         "--exon_cov_chains",
#         "-e",
#         action="store_true",
#         dest="exon_cov_chains",
#         help="If yes, computing the gene_overs parameter we consider only "
#         "those chains that cover at least 1 base in exons. If no, "
#         "consider all the chains overlapping and introns only any.",
#     )
#     app.add_argument("--log_file", type=str, help="Log file")
#     # print help if there are no args
#     if len(sys.argv) < 2:
#         app.print_help()
#         sys.exit(0)
#     args = app.parse_args()
#     return args


# def read_bed_data(bed_file):
#     """Get the necessary data from the bed file."""
#     result = {}  # return this dictionary
#     to_log(f"Reading {bed_file}")
#     f = open(bed_file, "r")
#     for line in f:
#         # parse tab-separated bed file
#         all_bed_info = line.rstrip().split("\t")
#         cds_track = make_cds_track(line)  # we need CDS only
#         cds_bed_info = cds_track.rstrip().split("\t")

#         if len(all_bed_info) != 12 or len(cds_bed_info) != 12:
#             # if there are not 12 fields - no guarantee that we parse what we want
#             die(f"Error! Bed12 file {bed_file} is corrupted!")

#         # extract fields that we need
#         chromStart = int(all_bed_info[1])  # gene start
#         chromEnd = int(all_bed_info[2])  # and end
#         # blocks represent exons
#         all_blockSizes = [int(x) for x in all_bed_info[10].split(",") if x != ""]
#         cds_blockSizes = [int(x) for x in cds_bed_info[10].split(",") if x != ""]
#         # data to save
#         gene_len = abs(chromStart - chromEnd)

#         # for decision tree I will need number of exons
#         # and number of bases in exonic and intronic fractions
#         exons_num = len(all_blockSizes)
#         exon_fraction = sum(all_blockSizes)  # including UTR
#         cds_fraction = sum(cds_blockSizes)  # CDS only
#         intron_fraction = gene_len - exon_fraction
#         gene_name = all_bed_info[3]

#         # save the data
#         result[gene_name] = {
#             "gene_len": gene_len,
#             "exon_fraction": cds_fraction,
#             "intron_fraction": intron_fraction,
#             "exons_num": exons_num,
#         }
#     f.close()
#     to_log(f"merge_chains_output: got data for {len(result.keys())} transcripts")
#     return result


# def process_gene_line(gene_line):
#     """Parse gene line."""
#     # contain [gene_id]=[intersected_chain_id] blocks
#     # remove "gene" and "", get (gene. chain) tuples
#     data_fields = [x for x in gene_line if "=" in x]
#     data = [(x.split("=")[1], x.split("=")[0]) for x in data_fields]
#     chain = data[0][0]  # chain is equal everywhere
#     genes = [x[1] for x in data]
#     return chain, genes


# def process_chain_line(chain_line):
#     """Parse chain line."""
#     chain_id = chain_line[1]
#     # collect the data in this dictionary
#     data_fields = {chain_id: {}}
#     # chain_id synteny global_score global_exo [local_exos] [coverage] [introns] [chain_len]

#     # simple values:
#     target = data_fields[chain_id]
#     target["synteny"] = chain_line[2]
#     target["global_score"] = chain_line[3]
#     target["global_exo"] = chain_line[4]
#     target["chain_Q_len"] = chain_line[5]
#     target["chain_len"] = chain_line[10]

#     # these values are more complex
#     # there are comma-separated pairs of values (TRANSCRIPT, VALUE),
#     target["local_exos"] = parse_pairs(chain_line[6])
#     target["coverages"] = parse_pairs(chain_line[7])
#     target["introns"] = parse_pairs(chain_line[8])
#     target["flank_cov"] = parse_pairs(chain_line[9])
#     return data_fields


# def load_results(results_dir):
#     """Load and sort the chain feature extractor results."""
#     to_log("merge_chains_output: Loading the results...")
#     results_files = os.listdir(results_dir)
#     to_log(f"merge_chains_output: There are {len(results_files)} result files to combine")

#     # to hold data from fields "genes":
#     chain_genes_data = defaultdict(list)
#     # to hold data from "chains" field:
#     chain_raw_data = {}
#     # read file-by-file, otherwise it takes too much place
#     genes_counter, chain_counter = 0, 0  # count chain and genes lines

#     for results_file in results_files:
#         # there are N files: read them one-by-one
#         path = os.path.join(results_dir, results_file)
#         f = open(path, "r")
#         for line in f:
#             # read file line-by-line, all fields are tab-separated
#             line_data = line.rstrip().split("\t")
#             # define the class of this line
#             # a line could be either gene or chain-related
#             if line_data[0] == "genes":
#                 # process as a gene line
#                 chain, genes = process_gene_line(line_data)
#                 chain_genes_data[chain].extend(genes)
#                 genes_counter += 1
#             elif line_data[0] == "chain":
#                 # chain related data
#                 the_chain_related = process_chain_line(line_data)
#                 # add this chain-related dict to the global one
#                 chain_raw_data.update(the_chain_related)
#                 chain_counter += 1
#         # do not forget to close the file
#         f.close()

#     to_log(f"merge_chains_output: got {len(chain_genes_data)} keys in chain_genes_data")
#     to_log(f"merge_chains_output: got {len(chain_raw_data)} keys in chain_raw_data")
#     to_log(f"merge_chains_output: There were {genes_counter} transcript lines and {chain_counter} chain lines")
#     # actually, these values must be equal
#     # just a sanity check
#     if not genes_counter == chain_counter:
#         err_msg = (
#             f"merge_chains_output: ERROR!\ngenes_counter and chain_counter hold different "
#             f"values:\n{genes_counter} and {chain_counter} respectively"
#         )
#         to_log(err_msg)
#         die("Some features extracting jobs died!")
#     return chain_genes_data, chain_raw_data


# def revert_dict(dct):
#     """Revert {a: [list of b's]} to {b: [list of a's]}."""
#     reverted = defaultdict(list)
#     for k, value in dct.items():
#         for v in value:
#             reverted[v].append(k)
#     to_log(f"merge_chains_output: chain_genes_data dict reverted, there are {len(reverted.keys())} keys now")
#     return reverted


# def make_synteny(genes, isoforms):
#     """Return synteny for a list of genes and dict of isoforms."""
#     return len(list(set([isoforms.get(gene) for gene in genes])))


# def combine(bed_data, chain_data, genes_data, exon_cov, isoforms):
#     """Combine chain and bed data and return gene-oriented table."""
#     to_log(f"merge_chains_output: Combining the data...")
#     combined = defaultdict(list)  # {gene: [related lines]}

#     for chain, data in chain_data.items():
#         # data dict contain the following keys:
#         # 'synteny', 'global_score', 'global_exo', 'chain_len', 'local_exos', 'coverages', 'introns'
#         # get the genes intersected this chain
#         genes = list(data["local_exos"].keys())
#         synteny = data["synteny"] if not isoforms else make_synteny(genes, isoforms)

#         for gene in genes:  # iterate gene-by-gene
#             # build a gene-oriented line
#             gene_feat = bed_data.get(gene)
#             if not gene_feat:
#                 continue  # it should not happen but...
#             # extract gene-related features from this chain
#             local_exo = data["local_exos"][gene]
#             exon_coverage = str(data["coverages"][gene])
#             intron_coverage = str(data["introns"][gene])
#             flank_cov = str(data["flank_cov"][gene])
#             # get a number of chains that covert this gene
#             # else: if you need the chains that overlap EXONS
#             gene_overs = (
#                 len(genes_data[gene])
#                 if not exon_cov
#                 else len([x for x in genes_data[gene] if x != "None"])
#             )
#             # fill the line
#             line_data = [  # gene and chain overlap information
#                 gene,
#                 gene_overs,
#                 chain,
#                 synteny,
#                 # global chain features
#                 data["global_score"],
#                 data["global_exo"],
#                 data["chain_len"],
#                 data["chain_Q_len"],
#                 # chain to gene local features
#                 local_exo,
#                 exon_coverage,
#                 intron_coverage,
#                 # gene features
#                 gene_feat["gene_len"],
#                 gene_feat["exons_num"],
#                 gene_feat["exon_fraction"],
#                 gene_feat["intron_fraction"],
#                 flank_cov,
#             ]

#             line = "\t".join([str(x) for x in line_data]) + "\n"
#             combined[gene].append(line)
#     # got all the data --> return back
#     to_log(f"merge_chains_output: got combined dict with {len(combined.keys())} keys")
#     return combined


# def save(data, output):
#     """Save the data into the file."""
#     # make the header
#     header_fields = (
#         "gene gene_overs chain synt gl_score gl_exo chain_len exon_qlen loc_exo exon_cover "
#         "intr_cover gene_len ex_num ex_fract intr_fract flank_cov".split()
#     )
#     header = "\t".join(header_fields) + "\n"
#     # define the stream to write the data
#     to_log(f"merge_chains_output: Writing output to {output}")
#     f = open(output, "w") if output != "stdout" else sys.stdout
#     f.write(header)
#     for lines in data.values():
#         f.write("".join(lines))
#     f.close()


# def merge_chains_output(
#     bed_file, isoforms_file, results_dir, output, exon_cov_chains=False
# ):
#     """Chains output merger core function."""
#     # read bed file, get gene features
#     bed_data = read_bed_data(bed_file)
#     # load isoforms data if provided
#     # returns 3 values, we keep only isoform to gene
#     if isoforms_file:
#         _, isoforms, _ = read_isoforms_file(isoforms_file)
#     else:
#         isoforms = None
#     # read result files from unit
#     chain_genes_data, chain_raw_data = load_results(results_dir)
#     # I need this dict reverted actually
#     # not chain-genes-data but gene-chains-data
#     genes_data = revert_dict(chain_genes_data)

#     # combine all the data into one gene-oriented dictionary
#     combined_data = combine(
#         bed_data, chain_raw_data, genes_data, exon_cov_chains, isoforms
#     )
#     # save this data
#     save(combined_data, output)
#     # finish the program
#     to_log(f"merge_chains_output: total runtime: {format(dt.now() - t0)}")

def parse_pairs(pairs: str) -> Dict[str, str]:
    """Parse lines like X=5,Y=56, and returns a dict."""
    return {k.strip():v.strip() for k,v in (x.split('=') for x in pairs.split(',') if x)}

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'results_dir',
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    metavar='RESULTS_DIR'
)
@click.argument(
    'bed_file',
    type=click.File('r', lazy=True),
    metavar='BED_FILE'
)
@click.argument(
    'output',
    type=click.File('w', lazy=True)
)
@click.option(
    '--isoforms',
    '-i',
    type=click.File('r', lazy=True),
    metavar='FILE',
    default=None,
    show_default=True,
    help=(
        'A file containing reference gene-to-isoform mapping'
    )
)
@click.option(
    '--exon_cov_chains',
    '-e',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, only the projections with at least one reference base covered '
        'are considered for further analysis'
    )
)
@click.option(
    '--log_name',
    '-ln',
    type=str,
    metavar='STR',
    default=None,
    show_default=True,
    help='Logger name to use; relevant only upon main class import'
)
@click.option(
    '--verbose',
    '-v',
    metavar='FLAG',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls the execution verbosity'
)

class FeatureAggregator(CommandLineManager):
    """
    Merges output files from the projection feature extraction step. Arguments are:\n
    * RESULTS_DIR is a directory containing output files from the feature extraction step;\n
    * BED_FILE is a reference annotation in BED format;\n
    * OUTPUT is a path to the output file
    """
    __slots__ = ('ref_bed', 'isoform2gene', 'gene2chains', 'chain_feature_data', 'only_covering')

    def __init__(
        self,
        results_dir: click.Path,
        bed_file: click.File,
        output: click.File,
        isoforms: Optional[Union[click.Path, None]],
        exon_cov_chains: Optional[bool],
        log_name: Optional[str],
        verbose: Optional[bool]
    ) -> None:
        self.v: bool = verbose
        self.set_logging(log_name)

        self.only_covering: bool = exon_cov_chains
        self.ref_bed: Dict[str, TranscriptFeatures] = {}
        self._to_log('Parsing reference BED file')
        self.parse_ref_bed(bed_file)
        self.isoform2gene: Dict[str, List[str]] = {}
        self.parse_ref_isoforms(isoforms)
        self.gene2chains: Dict[str, str] = defaultdict(list)
        self.chain_feature_data: Dict[str, List[str]] = {}
        self.load_features(results_dir)
        self._to_log('Feature data successfully aggregated')
        self._to_log('Writing output data')
        self.write_combined_output(output)


    def parse_ref_bed(self, file: TextIO) -> None:
        """Retrieves the necessary features from the reference BED file"""
        """Get the necessary data from the bed file."""
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 12:
                self._die('Reference BED file contains less than twelve fiels at line %i' % i)
            if data[1] == data[2] or data[6] == data[7]:
                self._to_log(
                    'Reference BED file contains an entry with no coding sequence at line %i' % i,
                    'warning'
                )
            thin_start: int = int(data[1])
            thin_end: int = int(data[2])
            name: str = data[3]
            cds_start: int = int(data[6])
            cds_end: int = int(data[7])
            exon_num: int = 0
            exon_sum: int = 0
            cds_sum: int = 0
            block_sizes: List[int] = [int(x) for x in data[10].split(',') if x]
            block_starts: List[int] = [int(x) for x in data[11].split(',') if x]
            if len(block_sizes) != len(block_starts):
                self._die('Number of block starts does not equal the number of block sizes at line %i' % i)
            for i in range(len(block_sizes)):
                block_start: int = thin_start + block_starts[i]
                block_end: int = block_start + block_sizes[i]
                exon_sum += block_sizes[i]
                exon_num += 1
                if block_end <= cds_start or block_start >= cds_end:
                    continue
                block_start = max(block_start, cds_start)
                block_end = min(block_end, cds_end)
                cds_sum += (block_end - block_start)
            gene_len: int = thin_end - thin_start
            intron_sum: int = gene_len - exon_sum
            self.ref_bed[name] = TranscriptFeatures(gene_len, cds_sum, intron_sum, exon_num)
        self._to_log('Extracted classification data for %i reference transcripts' % len(self.ref_bed))

    def parse_ref_isoforms(self, file: Union[TextIO, None]) -> None:
        """Parses reference isoform file"""
        if file is None:
            self._to_log('No isoform mapping file was provided')
            return
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split('\t')
            if not data or not data[0]:
                continue
            if len(data) < 2:
                self._die('Isoforms file contained less than two columns at line %i' % i)
            gene: str = data[0]
            tr: str = data[1]
            self.isoform2gene[tr] = gene

    def load_features(self, results_dir: click.Path) -> None:
        """Aggregates the results over multiple chain_runner.py runs stored in a single directory"""
        for file in os.listdir(results_dir):
            ## ignore the successful execution stamps
            if file.split('_')[-1] == OK:
                continue
            filepath: str = os.path.join(results_dir, file)
            with open(filepath, 'r') as h:
                for i, line in enumerate(h, start=1):
                    data: List[str] = line.rstrip().split('\t')
                    if not data or not data[0]:
                        continue
                    if data[0] == GENES:
                        split_gene_data: List[str] = [x.split('=') for x in data[1:] if '=' in x]
                        for gene, chain in split_gene_data:
                            self.gene2chains[gene].append(chain)
                    elif data[0] == CHAIN:
                        chain_id: str = data[1]
                        synteny: str = data[2]
                        global_score: str = data[3]
                        global_exon: str = data[4]
                        chain_q_len: str = data[5]
                        chain_len: str = data[10]
                        local_exons: Dict[str, str] = parse_pairs(data[6])
                        coverages: Dict[str, str] = parse_pairs(data[7])
                        introns: Dict[str, str] = parse_pairs(data[8])
                        flank_cov: Dict[str, str] = parse_pairs(data[9])
                        clipped_chain_qlen: str = data[11]
                        clipped_exons: Dict[str, str] = parse_pairs(data[12])
                        clipped_introns: Dict[str, str] = parse_pairs(data[13])
                        self.chain_feature_data[chain_id] = ChainFeatures(
                            synteny, global_score, global_exon,
                            chain_q_len, chain_len, local_exons, 
                            coverages, introns, flank_cov,
                            clipped_chain_qlen, clipped_exons, clipped_introns
                        )
                    elif TIMESTAMP in line:
                        continue
                    else:
                        self._die('Erroneous formatting at line %i in file %s' % (i, filepath))

    def write_combined_output(self, output: TextIO) -> None:
        """Combines the loaded data into output suitable for TOGA2 use"""
        output.write(HEADER + '\n')
        for chain, chain_features in self.chain_feature_data.items():
            trs: List[str] = chain_features.local_exons.keys()
            synteny: str = (
                chain_features.synteny if not self.isoform2gene else str(self._get_synteny(trs))
            )
            for tr in trs:
                tr_features: Union[TranscriptFeatures, None] = self.ref_bed.get(tr, None)
                if tr_features is None:
                    self._to_log('Transcript %s is missing from the reference annotation data' % tr)
                    continue 
                local_exon_score: str = chain_features.local_exons[tr]
                exon_coverage: str = str(chain_features.coverages[tr])
                intron_coverage: str = str(chain_features.introns[tr])
                flank_cov: str = str(chain_features.flank_cov[tr])
                # get a number of chains that covert this gene
                # else: if you need the chains that overlap EXONS
                gene_overs = (
                    len(self.gene2chains[tr])
                    if not self.only_covering
                    else len([x for x in self.gene2chains[tr] if x != 'None'])
                )
                ## get additional features for processed pseudogene classification
                clipped_introns: str = chain_features.clipped_introns[tr]
                line_data: List[str] = [
                    tr,
                    gene_overs,
                    chain,
                    synteny,
                    chain_features.global_score,
                    chain_features.global_exons,
                    chain_features.chain_len,
                    chain_features.chain_q_len,
                    local_exon_score,
                    exon_coverage,
                    intron_coverage,
                    tr_features.gene_len,
                    tr_features.exons_num,
                    tr_features.exon_fraction,
                    tr_features.intron_fraction,
                    flank_cov,
                    chain_features.clipped_chain_qlen,
                    clipped_introns
                ]
                output.write(
                    '\t'.join(map(str, line_data)) + '\n'
                )

    def _get_synteny(self, isoforms: List[str]) -> int:
        """Returns the number of unique isoforms for a given set of genes"""
        return len({self.isoform2gene[i] for i in isoforms if i in self.isoform2gene})


# def main():
#     """Entry point."""
#     args = parse_args()
#     setup_logger(args.log_file)
#     merge_chains_output(
#         args.bed_file,
#         args.isoforms,
#         args.results_dir,
#         args.output,
#         args.exon_cov_chains,
#     )


if __name__ == "__main__":
    # main()
    FeatureAggregator()
