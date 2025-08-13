#!/usr/bin/env python3
"""Script to run chain classification job.

Extract features from each chain to gene intersection.
Loads a list of chain: genes tasks and calls
modules.processor.unit for each chain: genes task.
"""
# import argparse
import os
import sys

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
sys.path.extend([LOCATION, PARENT])

from datetime import datetime as dt
from modules.overlap_select import overlap_select
from modules.common import bed_extract_id
from modules.common import make_cds_track
from modules.common import die
from modules.common import load_chain_dict
from modules.common import setup_logger
from modules.common import to_log
from modules.shared import CommandLineManager, CONTEXT_SETTINGS
from typing import Any, Dict, List, Optional, Tuple, Union
# from version import __version__

import click

__author__ = ('Bogdan M. Kirilenko', 'Yury V. Malovichko')

FLANK_SIZE: int = 10000  # gene flank size -> for flank ali feature
COMBINED_BED_ID: str = "COMBINED"  # placeholder gene name for intermediate tracks
COMBINED_CLIPPED_BED_ID: str = "COMBINED_CLIPPED"
ALL_EXONS_COMBINED: str = "ALL_EXONS_COMBINED"


def intersect(range_1, range_2):
    """Return intersection size."""
    return min(range_1[1], range_2[1]) - max(range_1[0], range_2[0])


def merge_ranges(range_1, range_2):
    """Return merged range."""
    return min(range_1[0], range_2[0]), max(range_1[1], range_2[1])


def extract_chain(chain_file, chain_dict, chain):
    """Extract chain string.

    We have: chain file, chain_id, start byte and offset.
    """
    f = open(chain_file, "rb")
    start, offset = chain_dict.get(int(chain))
    f.seek(start)  # jump to start_byte_position
    chain = f.read(offset).decode("utf-8")  # read OFFSET bytes
    f.close()
    return chain


def bed12_to_ranges(bed):
    """Convert bed-12 file to set of sorted ranges."""
    ranges_unsort, chrom = [], None
    for line in bed.split("\n"):
        # parse line and extract blocks
        line_info = line.split("\t")
        if not line_info or not line_info[0]:
            continue
        chrom = line_info[0]
        glob_start = int(line_info[1])
        blocks_num = int(line_info[9])
        block_sizes = [int(x) for x in line_info[10].split(",") if x != ""]
        block_starts = [
            glob_start + int(x) for x in line_info[11].split(",") if x != ""
        ]
        block_ends = [block_starts[i] + block_sizes[i] for i in range(blocks_num)]
        for i in range(blocks_num):  # save the range for each exon
            ranges_unsort.append((block_starts[i], block_ends[i]))
    # return sorted ranges
    die("(bed12_to_ranges) error, cannot read bed properly") if not chrom else None
    return chrom, sorted(ranges_unsort, key=lambda x: x[0])


def bedcov_ranges(ranges, chrom, clipped: bool = False):
    """Return a set of exons without overlaps.

    Python re-implementation of bedCov (kent) functionality.
    """
    ranges_filtered, pointer = [ranges[0]], 0  # initial values for filter
    gene = COMBINED_BED_ID if not clipped else COMBINED_CLIPPED_BED_ID  # if there is a mixture of genes - no ID anyway
    nested = False  # default value
    for i in range(1, len(ranges)):  # ranges are sorted so we can
        # compare each with only the next one
        if intersect(ranges[i], ranges_filtered[pointer]) <= 0:
            pointer += 1  # we have no intersection
            ranges_filtered.append(ranges[i])
        else:  # intersect - add merged range to the pointer
            # pointer = pointer - don't move it
            # replace the last range with merged last + new one
            nested = True  # at least one pair intersected
            ranges_filtered[pointer] = merge_ranges(ranges_filtered[pointer], ranges[i])

    # chr | start | end | gene | - bed4 structure
    # now make bed4 file
    exons, template = [], "{0}\t{1}\t{2}\t{3}\n"
    for grange in ranges_filtered:
        exons.append(template.format(chrom, grange[0], grange[1], gene))
    return exons, nested


def check_nest(work_data, cds_bed, clipped: bool = False):
    """Return True if genes are nested."""
    chrom, ranges = bed12_to_ranges(cds_bed)
    if not ranges:
        return False
    exons, nested = bedcov_ranges(ranges, chrom, clipped)
    work_data["exons" if not clipped else "clipped_exons"] = exons
    return nested


def get_tot_exons_track(work_data):
    """Get all exons including UTR and collapse them."""
    chrom, ranges = bed12_to_ranges(work_data["bed"])
    exons, _ = bedcov_ranges(ranges, chrom)
    # print(f'{ranges=}')
    # print(f'{exons=}')
    bed_template = "{0}\t{1}\t{2}\t{6}\t1000\t+\t{1}\t{2}\t0,0,0\t{3}\t{4}\t{5}"
    gene = ALL_EXONS_COMBINED
    blocks_uns = [(int(x.split("\t")[1]), int(x.split("\t")[2])) for x in exons]
    blocks = sorted(
        blocks_uns, key=lambda x: x[0]
    )  # no guarantee that it is sorted initially
    bed_12_start = min([x[0] for x in blocks])
    bed_12_end = max(x[1] for x in blocks)
    block_starts = ",".join([str(x[0] - bed_12_start) for x in blocks]) + ","
    block_sizes = ",".join([str(x[1] - x[0]) for x in blocks]) + ","
    bed_12 = bed_template.format(
        chrom, bed_12_start, bed_12_end, len(exons), block_sizes, block_starts, gene
    )
    return bed_12


def collapse_exons(work_data, clipped: bool = False):
    """Compensate nested genes."""
    # how bed12 looks like:
    # chr15	19964665	19988117	O	1000	+	19964665	19988117
    # 0,0,0	3	307,46,313,	0,390,22990,
    # I need to fill it with chrom, start and end and blocks info
    exon_key: str = "exons" if not clipped else "clipped_exons"
    if exon_key not in work_data:
        return
    bed_template = "{0}\t{1}\t{2}\t{6}\t1000\t+\t{1}\t{2}\t0,0,0\t{3}\t{4}\t{5}"
    chrom, _, _, gene = work_data[exon_key][0][:-1].split("\t")
    blocks_uns = [
        (int(x.split("\t")[1]), int(x.split("\t")[2])) for x in work_data[exon_key]
    ]
    # TODO: fix duplicated code fragment
    blocks = sorted(
        blocks_uns, key=lambda x: x[0]
    )  # no guarantee that it is sorted initially
    bed_12_start = min([x[0] for x in blocks])
    bed_12_end = max(x[1] for x in blocks)
    block_starts = ",".join([str(x[0] - bed_12_start) for x in blocks]) + ","
    block_sizes = ",".join([str(x[1] - x[0]) for x in blocks]) + ","
    bed_12 = bed_template.format(
        chrom,
        bed_12_start,
        bed_12_end,
        len(work_data[exon_key]),
        block_sizes,
        block_starts,
        gene,
    )
    nested_key: str = 'nested' if not clipped else 'clipped_nested'
    work_data[nested_key] = bed_12


def make_chain_clipped_track(
    bed_line: str, chain_start: int, chain_end: int
) -> Tuple[str, int]:
    """
    """
    data: List[str] = bed_line.rstrip().split('\t')
    chain_start, chain_end = sorted((chain_start, chain_end))
    init_start: int = int(data[1])
    init_end: int = int(data[2])
    new_start: int = init_start
    new_end: int = init_end
    data[3] = f'{data[3].replace("_CDS", "")}_chainclip'
    if init_start >= chain_start and init_end <= chain_end:
        exon_counter: int = int(data[9])
        return '\t'.join(data), exon_counter
    chrom_starts: List[int] = [int(x) for x in data[11].split(',') if x]
    chrom_sizes: List[int] = [int(x) for x in data[10].split(',') if x]
    new_exon_starts: str = ''
    new_exon_sizes: str = ''
    new_exon_counter: int = 0
    new_start_encountered: bool = False
    for i in range(len(chrom_starts)):
        start: int = init_start + chrom_starts[i]
        if start >= chain_end:
            break
        end: int = start + chrom_sizes[i]
        if end <= chain_start:
            continue
        new_exon_counter += 1
        start = max(start, chain_start)
        end = min(end, chain_end)
        if not new_start_encountered:
            new_start = start
            new_start_encountered = True
        new_exon_starts += f'{start - new_start},'
        new_exon_sizes += f'{end - start},'
    new_end = end
    if not new_exon_starts:
        new_exon_starts = '0,'
        new_exon_sizes = '0,'
        new_end = new_start
    # print(f'{data[3]=}, {init_start=}, {init_end=}, {chain_start=}, {chain_end=}, {new_start=}, {new_end=}')
    data[1] = str(new_start)
    data[6] = str(new_start)
    data[2] = str(new_end)
    data[7] = str(new_end)
    data[9] = str(new_exon_counter)
    data[10] = new_exon_sizes
    data[11] = new_exon_starts
    return '\t'.join(data), new_exon_counter


def make_grange_track(bed_line: str) -> str:
    """Converts BED12 track line into a 'genomic range' entry"""
    grange_track = bed_line.split("\t")  # tab-separated file
    # create the second track for the genomic region of the same gene
    # also known as "gene body"
    grange_track[3] = (
        grange_track[3] + "_grange"
    )  # I add _grange for the gene name, mark it
    grange_track[11] = "0"  # one block --> one start, starts from 0
    # size of block == size of the gene
    grange_track[10] = str(int(grange_track[2]) - int(grange_track[1]))
    grange_track[9] = "1"  # it means that it will be only one block
    return '\t'.join(grange_track)

def extend_bed_lines(bed_lines: str, chain_start: Optional[int] = None, chain_end: Optional[int] = None):
    """Create bed tracks for overlapSelect."""
    bed_lines_extended = ""  # init the variable to store the extended bed lines
    cds2covered_exons: Dict[str, int] = {}
    for line in bed_lines.split("\n")[:-1]:
        # verbose(f"Extending line:\n{line}")
        bed_lines_extended += line + "\n"  # first, I add the original bed line
        grange_track: str = make_grange_track(line)
        bed_lines_extended += grange_track + "\n"

        # create a separate track for FLANKS
        flanks_track = line.split("\t")
        name: str = flanks_track[3]
        flanks_track[3] = name + "_flanks"
        flanks_track[11] = "0"
        flanks_track[9] = "1"
        # need to avoid negative value here!
        # used 1 just for robustness
        flank_start = int(flanks_track[1]) - FLANK_SIZE
        flank_start = 1 if flank_start < 1 else flank_start
        flanks_track[1] = str(flank_start)
        flanks_track[6] = flanks_track[1]
        # no need to care about flank_end > chrSize
        # if not exists -> will be no blocks
        flanks_track[2] = str(int(flanks_track[2]) + FLANK_SIZE)
        flanks_track[7] = flanks_track[2]
        flanks_track[10] = str(int(flanks_track[2]) - int(flanks_track[1]))
        bed_lines_extended += "\t".join(flanks_track) + "\n"

        # add CDS track
        cds_track = make_cds_track(line)
        bed_lines_extended += cds_track + "\n"

        if chain_start is not None and chain_end is not None:
            clipped_cds, exon_counter = make_chain_clipped_track(cds_track, chain_start, chain_end)
            cds2covered_exons[name] = exon_counter
            bed_lines_extended += clipped_cds + '\n'
            clipped_cds_range: str = make_grange_track(clipped_cds)
            bed_lines_extended += clipped_cds_range + '\n'
    # print(bed_lines_extended)
    return bed_lines_extended, cds2covered_exons


# def cound_cds_exons(bed_lines_extended):
#     """Count CDS exons in each gene/transcript."""
#     bed_lines = [x.rstrip().split("\t") for x in bed_lines_extended.split("\n") if x != ""]
#     cds_lines = [x for x in bed_lines if x[3].endswith("_CDS")]
#     ret = {x[3][:-4]: int(x[9]) for x in cds_lines}
#     return ret


# def get_clipped_spans(bed_lines: str) -> Dict[str, Tuple[int, int]]:
#     """Gets total exon and intron lengths for chain-clipped BED records"""
#     out_dict: Dict[str, Tuple[int, int]] = {}
#     for line in bed_lines.split('\n'):
#         data: List[str] = line.rstrip().split('\t')
#         if not data or not data[0]:
#             continue
#         if not data[3].endswith('_chainclip'):
#             continue
#         cds_start: int = int(data[6])
#         cds_end: int = int(data[7])
#         exon_sum: int = sum(int(x) for x in data[10].split(',') if x)
#         intron_sum: int = (cds_end - cds_start) - exon_sum
#         gene: str = data[3].rstrip('_chainclip')
#         out_dict[gene] = (exon_sum, intron_sum)
#     return out_dict


def gene2clipped_spans(bed_lines: str, gene: str, start: int, stop: int) -> Tuple[int, int]:
    """
    Reports the total exon and intron length confined between the two coordinates
    """
    for line in bed_lines.split('\n'):
        data: List[str] = line.rstrip().split('\t')
        if not data:
            continue
        name: str = data[3]
        if name != f'{gene}_chainclip':
            continue
        # print(line.rstrip())
        cds_start: int = int(data[6])
        cds_end: int = int(data[7])
        start, stop = sorted((start, stop))
        total_span: int = min(stop, cds_end) - max(start, cds_start)
        exon_sizes: List[int] = [int(x) for x in data[10].split(',') if x]
        exon_starts: List[int] = [int(x) for x in data[11].split(',') if x]
        total_exon_sum: int = 0
        for i in range(len(exon_sizes)):
            exon_start: int = cds_start + exon_starts[i]
            exon_end: int = exon_start + exon_sizes[i]
            if exon_start > stop:
                break
            if exon_end < start:
                continue
            exon_start = max(exon_start, start)
            exon_end = min(exon_end, stop)
            total_exon_sum += (exon_end - exon_start)
        # print(f'{start=}, {stop=}, {cds_start=}, {cds_end=} {total_span=}')
        total_intron_sum: int = total_span - total_exon_sum
        return (total_exon_sum, total_intron_sum)


def get_features(
        work_data: Dict[str, Any], 
        result: Dict[str, Any], 
        bed_lines_extended: str,
        cds2cov_exons: Dict[str, int], 
        nested: bool = False
    ) -> None:
    """Compute local exon overlap score.

    For every line in the bed file (X - exonic base, I - intronic base)
    the new line will be created, represents the genomic region (called gene_grange).
    After the overlapSelect have been called, it returns the number of bases in chain blocks
    overlapping the gene exons:
    ----XXXXIIIIXXXXIIIIXXXX----- gene A - 5 overlapped bases / in exons and blocks
    ----XXXXXXXXXXXXXXXXXXXX----- gene A_grange - 9 overlapped bases / in all genomic region and blocks
    --bbbb----bbb-----bbbb------- chain N
    Here we consider the entire gene as a single exon
    In this toy example local_fractionExonOverlap score would be 5/9

    The raw overlapSelectOutput looks like:
    #inId   selectId        inOverlap       selectOverlap   overBases       similarity      inBases selectBases
    ENSG00000167232 chr17   0.112   1       399     0.201   3576    399
    ENSG00000167232_grange    chr17   0.0111  1       399     0.022   35952   399
    """
    # call overlap select
    chain_glob_bases, local_exo_dict, bed_cov_times, clipped_query_len, covered_cds_range = overlap_select(
        bed_lines_extended, work_data["chain"]
    )
    nums_of_cds_exons_covered = [
        len(v) for k, v in bed_cov_times.items() if k.endswith("_CDS")
    ]
    max_num_of_cds_exons_covered = (
        max(nums_of_cds_exons_covered) if len(nums_of_cds_exons_covered) > 0 else 0
    )
    # gene2clipped_spans: Dict[str, Tuple[int, int]] = get_clipped_spans(bed_lines_extended)
    # compute for each gene finally
    chain_cds_bases = 0  # summarize global set here
    clipped_chain_cds_bases: int = 0 ## summarize the same for 
    # print(f'{local_exo_dict=}')

    for gene in work_data["genes"]:
        # pick the data from overlap select table
        blocks_v_exons = local_exo_dict[gene]
        # print(f'{blocks_v_exons=}')
        blocks_v_cds = local_exo_dict[gene + "_CDS"]
        blocks_v_gene = local_exo_dict[gene + "_grange"]
        blocks_v_flanks_and_gene = local_exo_dict[gene + "_flanks"]
        # blocks_v_chain_clip: int = local_exo_dict.get(gene + "_chainclip", 0)
        blocks_v_chain_clip: int = local_exo_dict.get(gene + "_chainclip_grange", 0)
        blocks_v_exon_clip: int = local_exo_dict.get(gene + "_chainclip", 0)
        # print(f'{blocks_v_exons=}')
        # print(f'{blocks_v_cds=}')
        # print(f'{blocks_v_gene=}')
        # print(f'{bed_cov_times=}')
        # cds_exons_num = gene_to_cds_exons[gene]

        # all exons - CDS exons -> UTR exons
        blocks_v_utr_exons = blocks_v_exons - blocks_v_cds
        # gene blocks - UTR exons -> gene without UTR exons
        blocks_v_no_utr_exons = blocks_v_gene - blocks_v_utr_exons

        # if something like this happened -> there is a bug
        assert blocks_v_exons >= blocks_v_utr_exons
        assert blocks_v_exons >= blocks_v_cds

        # blocks gene + flanks - blocks gene -> blocks X flanks
        blocks_v_flanks = blocks_v_flanks_and_gene - blocks_v_gene
        # blocks gene - blocks exons -> blocks introns
        blocks_v_introns = blocks_v_gene - blocks_v_exons
        assert blocks_v_introns >= 0

        flank_feature = blocks_v_flanks / (FLANK_SIZE * 2)

        # global counters
        # CDS bases increase with blocks V cds in the gene
        chain_cds_bases += blocks_v_cds
        # increase number of UTR exons
        # chain_utr_exon_bases += blocks_v_utr_exons
        clipped_chain_cds_bases += blocks_v_exon_clip
        blocks_v_intron_clip = blocks_v_chain_clip - blocks_v_exon_clip

        # get local results
        result["gene_coverage"] += f"{gene}={blocks_v_cds},"
        result["gene_introns"] += f"{gene}={blocks_v_introns},"
        result["flanks_cov"] += f"{gene}={flank_feature},"

        # print(f'{covered_cds_range=}')
        clipped_cov_exon, clipped_cov_intron = gene2clipped_spans(
            bed_lines_extended, gene, *covered_cds_range
        )#gene2clipped_spans[gene]
        # print(
        #     f'{gene=}, {clipped_cov_exon=}, {clipped_cov_intron=}, {blocks_v_exon_clip=}, {blocks_v_intron_clip=}, {blocks_v_chain_clip=}, {covered_cds_range=}'
        # )
        clipped_exon_fraction: float = (
            blocks_v_exon_clip / clipped_cov_exon
        ) if clipped_cov_exon else 0
        if cds2cov_exons[gene] < 2:
            clipped_intron_fraction: float = -1.0
        else:
            clipped_intron_fraction: float = (
                blocks_v_intron_clip / clipped_cov_intron
            ) if clipped_cov_intron else 0
        # print(f'{clipped_query_len=}, {blocks_v_chain_clip=}, {blocks_v_exon_clip=}, {clipped_cov_exon=}, {blocks_v_intron_clip=}, {clipped_cov_intron=}')
        result["exon_cov_clipped"] += f"{gene}={clipped_exon_fraction},"
        result["intron_cov_clipped"] += f"{gene}={clipped_intron_fraction},"

        local_exo = (
            blocks_v_cds / blocks_v_no_utr_exons
            if blocks_v_no_utr_exons != 0.0
            else 0.0
        )
        assert local_exo >= 0
        assert local_exo <= 1
        result["local_exons"] += "{0}={1},".format(gene, local_exo)
        # increase synteny if > 0 CDS bases covered
        if blocks_v_cds > 0:
            result["chain_synteny"] += 1
            ov_block = f"{gene}={work_data['chain_id']}"
            result["gene_overlaps"].append(ov_block)
        else:
            # verbose(f"Chain don't overlap any exons in {gene}")
            # TO CHECK - it was like this here:
            result["gene_overlaps"].append(f"{gene}=None")
    # do not forget about global feature
    # chain_glob_bases -= chain_utr_exon_bases  # ignore UTR exons!
    # print(f'{chain_cds_bases=}, {local_exo_dict[COMBINED_BED_ID]=}')
    chain_v_all_exons = local_exo_dict[ALL_EXONS_COMBINED]
    chain_cds_bases = local_exo_dict[COMBINED_BED_ID] if nested else chain_cds_bases
    chain_v_utr_exons = chain_v_all_exons - chain_cds_bases
    # for x in work_data["exons"]:
    #     print(x)
    # for x in work_data["clipped_exons"]:
    #     print(x)
    # print(f'{clipped_query_len=}, {clipped_chain_cds_bases=}, {nested=}')
    # if COMBINED_BED_ID in local_exo_dict:
    #     print(f'{local_exo_dict[COMBINED_BED_ID]=}, {local_exo_dict[COMBINED_CLIPPED_BED_ID]=}')
    clipped_chain_cds_bases = local_exo_dict[COMBINED_BED_ID] if nested else clipped_chain_cds_bases
    q_len_corrected = work_data["chain_QLen"] - chain_v_utr_exons
    assert (
        q_len_corrected >= 0
    )  # chain length in query - blocks cover UTR cannot be a negative number
    result["global_exo"] = (
        chain_cds_bases / chain_glob_bases if chain_glob_bases != 0 else 0
    )
    # here we consider this transcript separately
    # nested genes do not affect this feature
    # chain_exon_bases = local_exo_dict[COMBINED_BED_ID] if nested else chain_exon_bases
    if max_num_of_cds_exons_covered > 1:
        result["Exlen_to_Qlen"] = (
            chain_cds_bases / q_len_corrected if q_len_corrected != 0 else 0
        )
        result["Exlen_to_Qlen_chainclip"] = (
            clipped_chain_cds_bases / clipped_query_len if clipped_query_len != 0 else 0
        )
    else:  # if chain covers at most 1 CDS exon -> this feature is not applicable
        result["Exlen_to_Qlen"] = 0
        result["Exlen_to_Qlen_chainclip"] = 0
    assert result["Exlen_to_Qlen"] <= 1  # blocklen / qlen cannot be > 1


def extract_cds_lines(all_bed_lines, clipped: bool = False):
    """Extract bed lines with names end with _CDS."""
    selected = []
    for line in all_bed_lines.split("\n"):
        if line == "":
            continue
        if line.split("\t")[3].endswith("_CDS" if not clipped else '_chainclip'):
            selected.append(line)
    if not selected:
        return ''
    return "\n".join(selected) + "\n"


def make_output(work_data, result, t0):
    """Arrange the output."""
    # verbose("Making the output...")
    chain_fields = [
        "chain",
        work_data["chain_id"],
        result["chain_synteny"],
        result["chain_global_score"],
        result["global_exo"],
        result["Exlen_to_Qlen"],
        result["local_exons"],
        result["gene_coverage"],
        result["gene_introns"],
        result["flanks_cov"],
        result["chain_len"],
        result["Exlen_to_Qlen_chainclip"],
        result["exon_cov_clipped"],
        result["intron_cov_clipped"]
    ]
    chain_output = "\t".join([str(x) for x in chain_fields]) + "\n"
    genes_output = "genes\t{0}\n".format("\t".join(result["gene_overlaps"]))
    time_output = f"#estimated time: {dt.now() - t0}\n"
    return chain_output, genes_output, time_output


def extended_output(result, t0):
    """Make human-readable output for small tests."""
    chain_output = "Chain-related features:\n"
    for key, value in result.items():
        if key == "gene_overlaps":
            continue
        chain_output += f'"{key}": {value}\n'
    genes_output = "These genes are overlapped by these chains:\n{0}".format(
        "\t".join(result["gene_overlaps"])
    )
    time_output = f"#estimated time: {dt.now() - t0}\n"
    return chain_output, genes_output, time_output


## TODO: Make it read a two-column file instead
@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'ref_hdf5',
    type=click.Path(exists=True),
    metavar='REF_HDF5'
)
@click.argument(
    'chain_file',
    type=click.Path(exists=True),
    metavar='CHAIN_FILE'
)
@click.option(
    '--input_file',
    '-i',
    type=click.File('r', lazy=True),
    metavar='INPUT_FILE',
    default=None,
    show_default=True,
    help=(
        'A tab-separated two-column file containing chain identifiers and '
        'respective comma-separated lists of projected transcripts. Deprecates '
        'the values provided with --chain_id and --transcripts options'
    )
)
@click.option(
    '--chain_id',
    '-c',
    type=str,
    metavar='CHAIN_ID',
    default=None,
    show_default=True,
    help=(
        'Id of the chain to analyse. Valid only if transcript list is provided '
        'instead of the input file'
    )
)
@click.option(
    '--transcripts',
    '-t',
    type=str,
    metavar='TRANSCRIPTS',
    default=None,
    show_default=True,
    help=(
        'A comma-separated list of transcript to analyse. Valid only if chain ID '
        'is provided instead of the input file'
    )
)
@click.option(
    '--extended',
    '-e',
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        'If set, reports output in human-readable format; '
        'recommended for standalone runs only'
    )
)
@click.option(
    '--output',
    '-o',
    type=click.File('a', lazy=True),
    metavar='OUTPUT_FILE',
    default=sys.stdout,
    show_default=False,
    help=(
        'A path to write the output to; the file is processed in the "add" mode '
        '[default: stdout]'
    )
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    default=False,
    show_default=True,
    help='Controls execution verbosity'
)


class ChainRunner(CommandLineManager):

    """
    A module for extracting projection features.\n
    This is a legacy code from TOGA 1.0 with minimal reflavouring. A proper
    update will be released after TOGA 2.0 pilot version is ready.\n\n

    Arguments are:\n
    * CHAIN_ID is a chain ID; only one chain is processed per each chain_runner
    instance;\n
    * TRANSCRIPTS is a comma-separated list of reference transcript IDs;\n
    * REF_HDF5 is an HDF5 file containing reference annotation data;\n
    * CHAIN_FILE is a local copy of alignment chain file. Note that the script
    uses the internal TOGA indices dumped in the same directory as the local copy;
    if, for any reason, you are running this script outside of the TOGA pipeline,
    make sure you are using the TOGA chain file located in tmp/input_data,
    otherwise index the input chain file manually.
    """

    __slots__ = [
        'ref_hdf5', 'chain_file', 'input_file', 'chain', 'transcripts',
        'chain2trs', 'index_file', 'extended', 'output', 'v',
        'chain_dict', 'logger'
    ]

    def __init__(
        self,
        ref_hdf5: click.Path,
        chain_file: click.Path,
        input_file: Union[click.File, None],
        chain_id: Union[str, None],
        transcripts: Union[str, None],
        extended: bool,
        output: click.File,
        verbose: bool
    ) -> None:
        self.v: bool = verbose
        self.set_logging()

        self.input_file: Union[click.File, None] = input_file
        self.chain: Union[str, None] = chain_id
        self.transcripts: Union[str, None] = transcripts
        self.chain2trs: Dict[str, List[str]] = {}
        self.check_input_consistency()
        if ref_hdf5[-4:] != 'hdf5':
            ref_hdf5 = '.'.join(ref_hdf5.split('.')[:-1]) + '.hdf5'
            if not os.path.isfile(ref_hdf5):
                self._die(
                    'Input reference annotation file is not an HDF5 instance '
                    'and does not have an explicitly named HDF5 counterpart'
                )
            self.ref_hdf5: str = ref_hdf5
        else:
            self.ref_hdf5: str = ref_hdf5
        ## CREATE A CHAIN ATTRIBUTE
        self.chain_file: str = chain_file
        index_file: str = chain_file + '_ID_position'
        if not os.path.isfile(index_file):
            self._die(
                'Chain file does not have a proper index or index file\'s name '
                'differs from the expected template'
            )
        self.index_file: str = index_file
        self.extended: bool = extended
        self.output: click.File = output
        self.run()

    def check_input_consistency(self) -> None:
        """
        Checks chain and transcript identifier arguments for consistency.
        If a (valid) input file is provided, its contents are used, overwriting
        any data provided via --chain and --transcripts.
        Otherwise, the method checks whether both chain ID and transcript names
        are provided.
        """
        if self.chain is None and self.transcripts is None and self.input_file is None:
            self._die(
                'No input was provided for chain_runner.py! Please provide either '
                '--input_file or --chain_id AND --transcripts options'
            )
        if self.input_file is not None:
            for line in self.input_file:
                data: List[str] = line.strip().split('\t')
                if not data:
                    continue
                if len(data) != 2:
                    self._die(
                        'Input file for chain_runner.py contains inconsistent '
                        'number of columns. Please make sure that argument for '
                        '--input_file contains two tab-separated columns'
                    )
                chain: str = data[0]
                trs: List[str] = [x for x in data[1].split(',') if x]
                if not trs:
                    self._die(
                        'Input file for chain_runner.py contains no '
                        f'transcripts for chain {chain}'
                    )
                self.chain2trs[chain] = trs
            return
        if self.chain is None or self.transcripts is None:
            self._die(
                'Either --chains or --transcripts were not provided. Please provide '
                'values for both options or set --input_file instead'
            )
        trs: List[str] = [x for x in self.transcripts.split(',') if x]
        self.chain2trs[chain] = trs

    def run(self) -> None:
        """
        Main running method
        """
        ## disable HDF5 file locking
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        ## read chain index file
        self.chain_dict: Dict[str, Tuple[int, int]] = load_chain_dict(self.index_file)

        ## run the extractor function
        ## NOTE: Code for chain extraction has been barely changed in TOGA 2.0;
        ## this function will likely be re-implemented as the main classe's method
        ## TODO: Rewrite as ChainRunner method
        for chain, trs in self.chain2trs.items():
            results: Tuple[str, str, str] = self.extract_chain_features(chain, trs)
            ## write the results
            # self.output.write(''.join(results[:2]))
            self.output.write(''.join(results))

    def extract_chain_features(
        self, chain: str, transcripts: List[str]
    ) -> Tuple[str, str]:
        """
        A method adaptation of extract_chain_features() function from TOGA 1.0
        """
        t0 = dt.now()
        ## LEGACY DATA STRUCTURES; TO BE MODIFIED
        work_data = {
            "bed": "",
            "chain": "",
            "nested": None,
            "clipped_nested": None,
            "chain_id": "",
            "chain_QLen": 0,
            "genes": [],
            "chain_len": 0,
            "chain_Tstarts": 0,
            "chain_Tends": 0,
            "chain_global_score": 0,
        }
        # structure to collect the results
        result = {
            "global_exo": 0.0,
            "flanks_cov": "",
            "gene_coverage": "",
            "gene_introns": "",
            "chain_synteny": 0,
            "local_exons": "",
            "gene_overlaps": [],
            "Exlen_to_Qlen": 0,
            "Exlen_to_Qlen_chainclip": 0,
            "exon_cov_clipped": "",
            "intron_cov_clipped": ""
        }

        ## populate :work_data: slots with respective values
        work_data['chain_id'] = chain
        work_data['bed'] = bed_extract_id(self.ref_hdf5, transcripts)
        work_data['genes'] = [
            x.split('\t')[3] for x in work_data['bed'].split('\n')[:-1]
        ] ## TODO: Outstandingly ugly; rewrite

        ## check if data for all transcripts required were successfully extracted
        if len(transcripts) != len(work_data['bed'].split("\n")[:-1]):
            self._echo('Warning. Not all the transcripts were found!')
            need_: int = len(transcripts)
            extracted_: int = len(work_data['bed'].split('\n')[:-1])
            self._echo(f'Expected {need_} transcripts, extracted {extracted_}')
            missing_genes: str = ','.join(
                [x for x in transcripts if x not in work_data['genes']]
            )
            self._echo(f'Missing transcripts:\n{missing_genes}')

        ## extract chain body from the file
        work_data['chain'] = extract_chain(self.chain_file, self.chain_dict, chain)

        ## parse chain header
        chain_header: List[str] = work_data['chain'].split('\n')[0].split()
        t_strand: bool = chain_header[4] == '+'
        q_start: int = int(chain_header[10])
        q_end: int = int(chain_header[11])
        t_start: int = int(chain_header[5])
        t_end: int = int(chain_header[6])
        # print(f'{t_start=}, {t_end=}, {t_strand=}')
        if not t_strand:
            t_size: int = int(chain_header[8])
            tstart: int = t_size - t_end
            t_end = t_size - t_start
            t_start = tstart
        q_len: int = abs(q_end - q_start)
        work_data['chain_QLen'] = q_len
        work_data['chain_Tstarts'] = int(chain_header[5])
        work_data['chain_Tends'] = int(chain_header[6])
        result['chain_global_score'] = int(chain_header[1])
        result['chain_len'] = work_data['chain_Tends'] - work_data['chain_Tstarts']

        ## computation part; to be modified
        bed_lines_extended, cds2cov_exons = extend_bed_lines(work_data['bed'], t_start, t_end)
        cds_bed_lines = extract_cds_lines(bed_lines_extended)
        tot_track: str = get_tot_exons_track(work_data)
        bed_lines_extended += f"{tot_track}\n"
        # print(bed_lines_extended)
        nested: bool = check_nest(work_data, cds_bed_lines)  # check if the genes are nested
        # print(bed_lines_extended)
        clipped_cds_lines: str = extract_cds_lines(bed_lines_extended, clipped=True)
        # print(f'{clipped_cds_lines=}')
        if clipped_cds_lines:
            _ = check_nest(work_data, clipped_cds_lines, clipped=True)

        if not nested:
            # 99% cases go here
            # there are no nested genes
            get_features(work_data, result, bed_lines_extended, cds2cov_exons)
        else:
            # another case, firstly need to make bed track with no intersections
            # and only after that call this function with flag NESTED for updated bed file
            collapse_exons(work_data)
            if work_data['nested'] is not None:
                bed_lines_extended += f"{work_data['nested']}\n"
            collapse_exons(work_data, clipped=True)
            if work_data['clipped_nested'] is not None:
                bed_lines_extended += f'{work_data["clipped_nested"]}\n'
            # print(bed_lines_extended)
            get_features(work_data, result, bed_lines_extended, cds2cov_exons, nested=True)
        if self.extended:
            # provide extended output
            # human-readable version
            output = extended_output(result, t0)
        else:
            # provide short version of output
            output = make_output(work_data, result, t0)
        return output





if __name__ == "__main__":
    ChainRunner()
