"""Python replacement for overlapSelect.

For a chain and a set of genes returns the following:
gene: how many bases this chain overlap in exons.
"""
from collections import defaultdict
from typing import Dict, List, Tuple
# from version import __version__

__author__ = "Bogdan Kirilenko, 2020."
__email__ = "bogdan.kirilenko@senckenberg.de"
__credits__ = ["Michael Hiller", "Virag Sharma", "David Jebb"]


def make_bed_ranges(bed_line):
    """Convert a bed line to a set of exon ranges."""
    line_info = bed_line.split("\t")
    # parse bed-12 line according to specification
    chrom_start = int(line_info[1])
    gene_name = line_info[3]
    # basically get exon absolute coordinates
    block_starts = [chrom_start + int(x) for x in line_info[11].split(",") if x != ""]
    block_sizes = [int(x) for x in line_info[10].split(",") if x != ""]
    blocks_num = int(line_info[9])
    block_ends = [block_starts[i] + block_sizes[i] for i in range(blocks_num)]
    # create and return (exon start, exon end, gene name) tuples
    genes = [gene_name for _ in range(blocks_num)]
    return list(zip(block_starts, block_ends, genes))


def parse_bed(bed_lines):
    """Return sorted genomic regions."""
    ranges = []
    for bed_line in bed_lines.split("\n"):
        if bed_line == "":
            continue
        gene_ranges = make_bed_ranges(bed_line)
        ranges.extend(gene_ranges)
    return list(sorted(ranges, key=lambda x: x[0]))


def chain_reader(chain):
    """Yield chain blocks one by one."""
    chain_data = chain.split("\n")
    chain_head = chain_data[0].split()
    del chain_data[0]  # keep blocks only in the chain_data list
    # define starting point
    progress = int(chain_head[5])
    query_strand: bool = chain_head[9] == '+'
    if query_strand:
        query_progress: int = int(chain_head[10])
    else:
        query_progress: int = int(chain_head[8]) - int(chain_head[10])
    blocks_num = len(chain_data)
    for i in range(blocks_num):
        # read block-by-block
        block_info = chain_data[i].split()
        if len(block_info) == 1:
            # only the last block has length == 1
            block_size = int(block_info[0])
            block_start = progress
            block_end = block_start + block_size
            query_block_start: int = query_progress
            query_block_end: int = query_block_start + (block_size if query_strand else - block_size)#block_size
            query_block_start, query_block_end = sorted((query_block_start, query_block_end))
            yield block_start, block_end, query_block_start, query_block_end
            break  # end the loop

        # get intermediate block data
        block_size = int(block_info[0])
        dt = int(block_info[1])
        block_start = progress
        block_end = block_start + block_size
        progress = block_end + dt
        dq: int = int(block_info[2])
        query_block_start: int = query_progress
        query_block_end: int = query_block_start + (block_size if query_strand else -block_size)#block_size
        query_progress = query_block_end + (dq if query_strand else - dq)#dq
        query_block_start, query_block_end = sorted((query_block_start, query_block_end))
        yield block_start, block_end, query_block_start, query_block_end


def intersect(ch_block, be_block):
    """Return intersection size."""
    return min(ch_block[1], be_block[1]) - max(ch_block[0], be_block[0])


def gene2strand_dict(bed: str) -> Dict[str, bool]:
    """Creates a sequence:strand dictionary"""
    out_dict: Dict[str, str] = {}
    for line in bed.split('\n'):
        data: List[str] = line.split('\t')
        if not data or not data[0]:
            continue
        name: str = data[3]
        strand: str = data[5]
        out_dict[name] = strand == '+'
    return out_dict

def overlap_select(bed, chain):
    """Python implementation of some overlapSelect (kent) functionality."""
    ranges = parse_bed(bed)  # list of exon ranges
    genes = [x[2] for x in ranges]
    gene2strand: Dict[str, str] = gene2strand_dict(bed)

    # bed overlaps: our results, count overlapped bases per gene
    bed_overlaps = {gene: 0 for gene in genes}
    bed_covered_times = defaultdict(
        set
    )  # for each bed track: how many exons intersected
    chain_len = 0  # sum of chain blocks
    start_with = 0  # bed tracks counter

    chain_header: List[str] = chain.split('\n')[0].split(' ')
    query_strand: bool = chain_header[9] == '+'
    # t_start: int = int(chain_header[5])
    # t_end: int = int(chain_header[6])
    # t_strand: bool = chain_header[4] == '+'
    # if not t_strand:
    #     t_size: int = int(chain_header[3])
    #     _t_start: int = t_size - t_end
    #     t_end = t_size - t_start
    #     t_start = _t_start
    # cds_start: int = min(
    #     [x[0] for x in cds_ranges]
    # )
    # cds_end: int = max([x[1] for x in cds_ranges])
    first_cds_block: Tuple[int, int, int, int] = None
    last_cds_block: Tuple[int, int, int, int] = None
    first_cds_exon: Tuple[int, int, str] = None
    last_cds_exon: Tuple[int, int, str] = None

    # init blocks generator
    for block in chain_reader(chain):
        # add to len
        chain_len += block[1] - block[0]  # block[1] - block[0] is block length
        FLAG = False  # was there an intersection or not?
        FIRST = True
        bed_num = 0
        while True:
            if FIRST:  # start with previous start, first iteration
                bed_num = start_with
                FIRST = False  # guarantee that this condition works ONCE per loop
            else:  # just increase the pointer
                bed_num += 1  # to avoid inf loop

            if bed_num >= len(ranges):
                # we intersected all beds, they are over
                break
            exon = ranges[bed_num]  # pick the exon
            if block[1] < exon[0]:  # too late
                break  # means that bed is "righter" than chain block

            # get intersection between the chain block and exon
            block_vs_exon = intersect(block, (exon[0], exon[1]))
            if block_vs_exon > 0:  # they intersect
                # if exon[2].endswith('_chainclip'):
                #     print(f'{block=}, {exon=}')
                if not FLAG:  # the FIRST intersection of this chain block
                    # shift the start: all exons left-side will never be reached

                    # guarantee that I will assign to starts with
                    # only the FIRST intersection (if it took place)
                    start_with = bed_num
                    FLAG = True  # otherwise starts_with will be preserved
                # add the intersection size
                bed_overlaps[exon[2]] += block_vs_exon
                bed_covered_times[exon[2]].add(exon[0])
                if exon[2].endswith('_chainclip'):
                    if first_cds_block is None:
                        first_cds_block = block
                        first_cds_exon = exon
                    last_cds_block = block
                    last_cds_exon = exon

            else:  # we recorded all the region with intersections
                if block[0] > exon[1]:  # too early
                    # in case like:
                    # gene A EEEEE----------------------------------------EEEEEE #
                    # gene B               EEEEEEEEEE                            #
                    # gene C                               EEEEEEEEE             #
                    # chain                                    ccccc             #
                    # at gene A I will get FLAG = True and NO intersection with gene B
                    # --> I would miss gene C in this case without this condition.
                    continue

                elif FLAG:  # this is not a nested gene
                    break  # and all intersections are saved --> proceed to the next chain

    ## define the coding boundaries in the query
    if first_cds_block is not None:
        cds_start_offset: int = first_cds_exon[0] - first_cds_block[0]
        cds_end_offset: int = last_cds_block[1] - last_cds_exon[1]
        if query_strand:#gene2strand[first_cds_exon[2]] == query_strand:
            cds_start_in_query: int = first_cds_block[2] + cds_start_offset
            cds_end_in_query: int = last_cds_block[3] - cds_end_offset
        else:
            cds_start_in_query: int = first_cds_block[3] - cds_start_offset
            cds_end_in_query: int = last_cds_block[2] + cds_end_offset
        # if query_strand:#gene2strand[last_cds_exon[2]] == query_strand:
        #     cds_end_in_query: int = last_cds_block[3] - cds_end_offset
        # else:
        #     cds_end_in_query: int = last_cds_block[2] + cds_end_offset
        # print(f'{first_cds_block=}, {last_cds_block=}')
        # print(f'{first_cds_exon=}, {last_cds_exon=}')
        # print(f'{cds_start_in_query=}, {cds_end_in_query=}')
        query_coding_len: int = abs(cds_end_in_query - cds_start_in_query)
        covered_cds_range: Tuple[int, int] = (
            max(first_cds_exon[0], first_cds_block[0]), 
            min(last_cds_exon[1], last_cds_block[1])
        )
    else:
        query_coding_len: int = 0
        covered_cds_range: int = (0, 0)
    # return the required values
    return chain_len, bed_overlaps, bed_covered_times, query_coding_len, covered_cds_range


# def get_coding_boundaries_in_query(
#     bed: str,
#     chain: str
# ) -> Tuple[int, int]:
#     ranges: List[Tuple[int, int, str]] = parse_bed(bed)
#     cds_ranges: List[Tuple[int, int, str]] = [
#         x for x in ranges if x[2].endswith('_CDS')
#     ]
#     chain_header: List[str] = chain.split('\n')[0].split(' ')
#     t_start: int = int(chain_header[5])
#     t_end: int = int(chain_header[6])
#     t_strand: bool = chain_header[4] == '+'
#     if not t_strand:
#         t_size: int = int(chain_header[3])
#         _t_start: int = t_size - t_end
#         t_end = t_size - t_start
#         t_start = _t_start
#     cds_start: int = min(
#         [x[0] for x in cds_ranges]
#     )
#     cds_end: int = max([x[1] for x in cds_ranges])
#     for x in 
