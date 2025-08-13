"""
The Cython source for the projection extraction script.
This functionale allows one to:
  1) find the chains overlapping the target search space
  2) extract the coordinates in the query corresponding to the target search space
  3) give the coordinates of all the exons covered by the search space, and
  4) give the coordinates of all the exon projections to the query
"""

import cython
cimport cython

__author__ = "Yury V. Malovichko"
__date__ = "2022"
__version__ = "0.1"
__email__ = "yury.malovichko@senckenberg.de"

## this is a dummy constant for the orthology indicator
## correc once you obtain a valid chain classification table
cdef double ORTHOLOGY_THRESHOLD = 0.5

cdef list getChainHeaders_(str file, list allowed_names):
    cdef str line
    cdef list chain_headers = []
    with open(file, 'r') as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('chain'):
                if allowed_names and line.split()[-1] not in allowed_names:
                    continue
                chain_headers.append(line)
#    chain_headers = sorted(chain_headers, key=lambda x: int(x.split()[2]))
    return chain_headers


def getChainHeaders(file: str, allowed_names: list = []) -> list:
    return getChainHeaders_(file, allowed_names)


cdef dict parseBed_(str file):
    cdef dict bed_annot = dict()
    cdef str line, chrom, name
    cdef int start, stop, exon_start, exon_size, strand
    cdef list line_split, exon_starts, exon_sizes
    cdef tuple transcript_record, exons

    with open(file, 'r') as handle:
      for line in handle:
        line_split = line.rstrip().split('\t')
        ## TODO: consider creating a user-defined class
        chrom = line_split[0]
        ## should I consider the whole gene or only the coding sequence?
        start = int(line_split[1])
        stop = int(line_split[2])
        transcript_id = line_split[3]
        strand = 1 if line_split[5] == '+' else 0
        exon_starts = [ int(x) for x in line_split[11].split(',')[:-1] ]
        exon_sizes = [
             int(x) for x in line_split[10].split(',')[:-1]
        ]
        exons = tuple(
            map(
                lambda x, y, z: (x, x + y, z),
                exon_starts,
                exon_sizes,
                [strand] * len(exon_starts)
            )
        )
        transcript_record = (transcript_id, start, stop, chrom, strand, exons)
        ## TODO: consider making the output dict nested, with the second level
        ## standing for the strand
        bed_annot[transcript_id] = transcript_record
    ## sort the resulting dictionary to easen the intersection search
    ## thank God the dictionaries are now ordered in Python
    bed_annot = dict(sorted(bed_annot.items(), key=lambda x: (x[1][3], x[1][1])))
    return bed_annot


def parseBed(file: str) -> dict:
    return parseBed_(file)


cdef list parseFocalSpaces_(str file):
      cdef str line
      cdef list focal_list = []
      with open(file, 'r') as handle:
          for line in handle:
              focal_list.append(line.strip())
      return focal_list


def parseFocalSpaces(file: str):
    return parseFocalSpaces_(file)


cdef dict parseClassdChains_(str file):
##    If the chain classification is provided,
##    store the list of the approved chains
    cdef dict trans_to_chain_map = dict()
    cdef tuple record
    cdef list line_info
    cdef str line, chain_id, transcript_id
    cdef float orthology_score
    with open(file, 'r') as handle:
        for line in handle:
            if line[:10] == 'transcript':
                continue
            line_info = line.strip().split('\t')
            orthology_score = float(line_info[2])
            if orthology_score >= ORTHOLOGY_THRESHOLD:
                transcript_id = line_info[0]
                chain_id = line_info[1]
                record = (chain_id, orthology_score)
                if transcript_id in trans_to_chain_map.keys():
                    trans_to_chain_map[transcript_id].append(record)
#                    if orthology_score > trans_to_chain_map[transcript_id][1]:
#                        trans_to_chain_map[transcript_id] = record
                else:
                    trans_to_chain_map[transcript_id] = [record] ## TODO: a very ugly if-else solution, consider rewriting
    return trans_to_chain_map


def parseClassifiedChainFile(file: str) -> dict:
    return parseClassdChains_(file)


cdef list parseChains(str chain_header, str chain_file, str source):
##    Captures the desired chain, then parses its coverage pattern
##    in the genome of interest
    cdef bint proceed, chain_space = True, False
    cdef str line, strand
    cdef list line_info
    cdef int block_ind, block_size, block_start, block_stop, chrom_size
    cdef list block_list = list()
    with open(chain_file, 'r') as handle:
        while proceed:
            line = next(handle).strip()
            if line == chain_header: # catch the desired chain
                chain_space = True # indicates the start of the chain body
                line_info = line.split()
                if source == 'target':
                    # set the chrom size, strand, chain start and the distance column
                    chrom_size = line_info[3]
                    strand = line_info[4]
                    block_start = int(line_info[5])
                    block_ind = 1
                elif source == 'query':
                    chrom_size = line_info[8]
                    strand = line_info[9]
                    block_start = int(line_info[10])
                    block_ind = 2
                if strand == "-":
                    block_start = chrom_size - block_start + 1
            elif chain_space: # start the chain body parsing
                line_info = list(map(int, line.split('\t')))
                block_size = line_info[0]
                if strand == "-":
                    block_stop = block_start - block_size + 1
                else:
                    block_stop = block_start + block_size
                block_list.append((block_start, block_stop))
                if len(line_info) > 1:
                    # if the distances to the next block are given,
                    # increase the start position
                    block_start += (block_size + line_info[block_ind])
                else:
                    # otherwise the chain body is over
                    proceed = False

    return block_list


cdef bint checkChainSpan(str chain_header, tuple transcript_record):
##    Checks whether the chain spans over the transcript
    cdef list chain_info = chain_header.split()
    cdef int chain_start = <int>chain_info[5]
    cdef int chain_end = <int>chain_info[6]
    return (min(chain_end, transcript_record[3]) - \
            max(chain_start, transcript_record[2])) > 0


cdef bint intersection(int start1, int end1, int start2, int end2):
##    An abstract version of the intersection checker
    return (min(end1, end2) - max(start1, start2)) > 0
    ## TODO: should the exon be completely covered by the aligned block?


cdef bint checkBlockCoverage(list chain_blocks, tuple transcript_record):
##    Checks whether the genomic entity is covered by
##    any aligned block of the chain
    cdef tuple block
    cdef int block_number = len(chain_blocks)
    cdef int exon_number = len(transcript_record[4])
    cdef int block_index, exon_index
    cdef int block_start, block_end
    cdef int exon_start, exon_end
    cdef int current_block = 0
    cdef list values = []
    ## assume that the exon coordinate tuples are properly sorted
    for exon_index in range(exon_number):
        exon_start, exon_end = transcript_record[4][exon_index]
        ## also, assume that the blocks are properly sorted
        for block_index in range(current_block, block_number):
            block_start, block_end = chain_blocks[block_index]
            ## once the intersecting block is captured, the upstream blocks
            ## pose no interest (provided the blocks are sorted)
            if intersection(block_start, block_end,
               exon_start, exon_end):
                values.append(True)
                current_block = block_index
                continue
            values.append(False)
    return all(values) ## TODO: should this line be Cythonized with a cycle?


cdef tuple chain_bed_intersect(list chains, dict bed, str chain_file):
##    For each transcript, assesses whether it is covered by any chain,
##    then test whether each of the transcript's exons is covered by the aligned
##    block.
    cdef tuple transcript
    cdef tuple map_entry
    cdef bint all_transcripts_covered
    cdef chain_dict = dict()
    cdef list c2b_mapping = list
    cdef list skipped = list()
    cdef set missing_chroms, present_chroms
    cdef str chain_header, chain_key, bed_key, chrom, missing_chrom, chain
    cdef int chrom_ind, trans_ind, chain_ind
    cdef int chain_anchor
    ## TODO: consider moving this block to the chain header parser
    for chain_header in chains:
        chrom = chain_header.split()[2]
        if chrom not in chain_dict.keys():
            chain_dict[chrom] = []
        chain_dict[chrom].append(chain_header)
    ## Bogdan's feature: track the transcripts from the chromosomes
    ## absent in either bed or chain file
    missing_chroms = set(chain_dict.keys()).difference(bed.keys())
    for missing_chrom in missing_chroms:
        for transcript in bed[missing_chrom]:
            skipped.append(transcript)
    present_chroms = set(chain_dict.keys()).difference(missing_chroms)
    ## now, find the spanning chain for each transcript
    for chrom_ind in range(len(present_chroms)):
        chrom = present_chroms[chrom_ind]
        chain_anchor = 0
        for trans_ind in range(len(bed[chrom])):
            transcript = bed[chrom][trans_ind]
            for chain_ind in range(chain_anchor, len(chain_dict[chrom])):
                chain = chain_dict[chrom]
                if checkChainSpan(chain, transcript):
                    ## we have found the spanning chain
                    ## now parse its block content
                    chain_blocks = parseChains(chain, chain_file, 'target')
                    all_transcripts_covered = checkBlockCoverage(chain_blocks,
                                                                 transcript
                    )
                    chain_id = chain.split()[-1]
                    map_entry = (transcript[1], chain_id, all_transcripts_covered)
                    c2b_mapping.append(map_entry)
                    ## this chain, however, can still span over the downstream
                    ## transcripts, so keep it as the iteration start
                    chain_anchor = chain_ind
                    ## go to the next transcript
                    continue
    return (c2b_mapping, skipped)


cdef int focalInd_(list transcript_entries, str focal_transcript):
    cdef int current_ind = 0, focal_ind
    cdef proceed = True
    cdef str transcript_record
    while current_ind < len(transcript_entries):
        transcript_record = transcript_entries[current_ind]
        if transcript_record == focal_transcript:
            focal_ind = current_ind
            return focal_ind
        current_ind += 1


def focalInd(
    transcript_entries: list,
    focal_transcript: str
) -> int:
    return focalInd_(transcript_entries, focal_transcript)


cdef dict parseExonFile_(str file):
    cdef str line
    cdef list line_data
    cdef dict out_dict = dict()
    with open(file, 'r') as h:
      for line in h.readlines():
        line = line.strip()
