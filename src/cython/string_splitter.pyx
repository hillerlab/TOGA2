import cython
cimport cython

cdef str EMPTY_STR = ''
cdef str NL = '\n'
cdef str PLUS = '+'
cdef SPACE = ' '
cdef TAB = '\t'

def split_at_(str string, str sep):
    cdef str line = EMPTY_STR
    cdef str symbol = EMPTY_STR
    for symbol in string:
        if symbol == sep:
            yield line
            line = EMPTY_STR
            continue
        line += symbol
    if line:
      yield line

cdef list split_at(str string, str sep):
    cdef str line = EMPTY_STR
    cdef str symbol = EMPTY_STR
    cdef list output = []
    for symbol in string:
        if symbol == sep:
            output.append(line)
            line = EMPTY_STR
            continue
        line += symbol
    if line:
        output.append(line)
    return output

cpdef transcriptwise_subchains(str chain, list names, list starts, list stops):
    cdef str line = EMPTY_STR
    cdef str symbol = EMPTY_STR
    cdef list block
    cdef str ref_chr, query_chr
    cdef int ref_chr_size, ref_start, ref_end, query_chr_size, query_start, query_end
    cdef bint ref_strand, query_strand
    cdef int size, dt, dq
    cdef int r_start = 0
    cdef int r_end = 0
    cdef int q_start = 0
    cdef int q_end = 0
    cdef int b_num = 1
    cdef int tr_num = len(names)
    cdef int tr_start, tr_end
    cdef int prev_tr_end = 0
    cdef int ref_first, ref_second
    cdef int t_chain_start, t_chain_stop
    cdef str gap_num = EMPTY_STR
    cdef bint up_block, down_block, discard_next_block, block_added, gap_added
    cdef bint prev_end_exceeded = False
    cdef int curr_tr = 0

    cdef abs_start = starts[0]
    cdef abs_stop = max(stops)
    cdef int[4] gap_coord_tuple# = (0,0,0,0)
    cdef int[4] coord_tuple# = (0,0,0,0)
    cdef bint gap_coords_exist = False
    cdef bint coords_exist = False ## TODO: Not sure if this flag is really needed

    cdef tuple chain_meta
    cdef dict tr2blocks = {}

    cdef str tr
    cdef int i

    cdef int counter = 0
    cdef bint first_block, last_block

    for tr in names:
        tr2blocks[tr] = {}

    # print(f'{names}\n{starts}\n{stops}')
    for symbol in chain:
        if symbol == NL:
            # print(line)
            if line == EMPTY_STR:
                break
            if line[:5] == 'chain':
                block = split_at(line, SPACE)
                ref_chr = block[2]
                ref_chr_size = int(block[3])
                ref_strand = block[4] == PLUS
                ref_start = int(block[5])
                ref_end = int(block[6])
                if not ref_strand:
                    t_chain_start: int = ref_chr_size - ref_end
                    t_chain_stop: int = ref_chr_size - ref_start
                    ref_start = ref_chr_size - ref_start
                    ref_end = ref_chr_size - ref_end
                    ## TODO: Think properly how to change iteration at this point
                    names = names[::-1]
                    starts = starts[::-1]
                    stops = stops[::-1]
                    # abs_start = starts[0]
                    # abs_stop = stops[-1]
                else:
                    t_chain_start, t_chain_stop = ref_start, ref_end
                query_chr = block[7]
                query_chr_size = int(block[8])
                query_strand = block[9] == PLUS
                query_start = int(block[10])
                query_end = int(block[11])
                if not query_strand:
                    query_start = query_chr_size - query_start
                    query_end = query_chr_size - query_end
                chain_meta = (
                    ref_chr, query_chr, t_chain_start, t_chain_stop, ref_strand, query_start, query_end, query_strand
                )
                # print('chain_meta=', chain_meta)
                line = EMPTY_STR
                continue
            block = split_at(line, TAB)
            counter += 1
            # print(block)
            size = int(block[0])
            if len(block) == 3:
                dt = int(block[1])
                dq = int(block[2])
            else:
                dt = 0
                dq = 0

            first_block = counter == 1
            last_block = len(block) != 3

            if not r_start:
                r_start = ref_start
            if not q_start:
                q_start = query_start

            if ref_strand:
                r_end = r_start + size
            else:
                r_end = r_start - size

            if query_strand:
                q_end = q_start + size
            else:
                q_end = q_start - size

            if r_end > r_start:
                ref_first = r_start
                ref_second = r_end
            else:
                ref_first = r_end
                ref_second = r_start
            if q_end > q_start:
                query_first = q_start
                query_second = q_end
            else:
                query_first = q_end
                query_second = q_start
            coord_tuple = (
                ref_first, ref_second, query_first, query_second
            )
            coords_exist = True
            ## reset the current start values
            r_start = r_end
            q_start = q_end
            ## first, check if the block intersects the transcript array
            up_block = coord_tuple[1] < abs_start
            down_block = coord_tuple[0] > abs_stop

            # print('coord_tuple=', coord_tuple, ', abs_start=', abs_start, ', abs_stop=', abs_stop, ', up_block=', up_block, ', down_block=', down_block)

            ## if iterator has not yet reached the first transcript's search space,
            ## there is no need to check for further intersections;
            ## however, in order to ensure that every transcript has a space to align,
            ## create a provisional gap object
            if (ref_strand and up_block) or (not ref_strand and down_block):
                if ref_strand:
                    r_end = r_start + dt
                else:
                    r_end = r_start - dt
                if query_strand:
                    q_end = q_start + dq
                else:
                    q_end = q_start - dq
                # if not gap_coord_tuple:
                gap_num = '0_0'
                gap_coord_tuple = (r_start, r_end, q_start, q_end)
                gap_coords_exist = True
                r_start = r_end
                q_start = q_end
                line = EMPTY_STR
                continue
            # print(f'Adding the block {b_num}, curr_tr={curr_tr}')
            block_added = False
            ## if iterator goes past the last transcript's locus, the process halts;
            ## however, for spanning chains the upstream gap should still be added
            if ((ref_strand and down_block) or (not ref_strand and up_block)):
                # print('The loop has passed the last transcript')
                discard_next_block = True
                for i in range(curr_tr, len(names)):
                    tr = names[i]
                    # if tr == 'ENST00000541819.GABRB3':
                    #     print(f'UBUB tr2blocks[tr]={tr2blocks[tr]}')
                    # print(f'tr={tr}, tr2blocks[tr]={tr2blocks[tr]}, tr_start={starts[i]}, tr_stop={stops[i]}, block={coord_tuple}, abs_start={abs_start}, abs_stop={abs_stop}')
                    # print(f'tr={tr}, gap_num in tr2blocks[tr]={gap_num in tr2blocks[tr]}, b_num not in tr2blocks[tr]={b_num not in tr2blocks[tr]}')
                    if gap_num in tr2blocks[tr] and b_num not in tr2blocks[tr]:
                        discard_next_block = False
                        break
                # print(f'b_num={b_num}, gap_num={gap_num} discard_next_block={discard_next_block}')
                if discard_next_block:
                    for i in range(curr_tr, len(names)):
                        tr = names[i]
                        if not tr2blocks[tr]:
                            # print(f'Adding gap {gap_num} for transcript {tr} before breaking the loop')
                            if not gap_coords_exist:
                                continue
                            # print(f'gap_num={gap_num}, gap_coord_tuple={gap_coord_tuple}, gap_coords_exist={gap_coords_exist}')
                            tr2blocks[tr][gap_num] = gap_coord_tuple
                    line = EMPTY_STR
                    # print(f'tr2blocks[tr]={tr2blocks[tr]}')
                    # print('BBBBBREAK')
                    break
            ## check if the block intersects any of the transcripts
            # print(f'Intersecting the block {b_num}')
            for i in range(curr_tr, len(names)):
                tr = names[i]
                tr_start = starts[i]
                tr_end = stops[i]
                prev_end_exceeded = tr_end >= prev_tr_end
                if prev_end_exceeded:
                    prev_tr_end = tr_end
                # print(f'{tr}, {tr_start}, {tr_end}, {tr_end-tr_start}, {b_num}, {coord_tuple}, {curr_tr}, {ref_strand}, {tr_num}')
                if ref_strand:
                    if tr_end < coord_tuple[0]:
                        ## the transcript lies upstream to the aligned block;
                        ## update the transcript pointer and proceed further
                        if prev_end_exceeded:
                            curr_tr += 1
                        # print(f'Updating counter: {tr}, {tr_end}, {coord_tuple}')
                        if curr_tr < tr_num:
                            abs_start = starts[curr_tr]
                        ## TODO: Uncomment the following print tomorrow;
                        ## according to it, excessive number of blocks is saved per transcript,
                        ## but why???
                        # print(f'{gap_num}, {tr2blocks[tr]}, {gap_num in tr2blocks[tr]}')
                        # if gap_num and gap_num in tr2blocks[tr]:
                        #     tr2blocks[tr][str(b_num)] = coord_tuple
                        if not tr2blocks[tr]:
                            # print(f'Adding emergency gap: {gap_coord_tuple}')
                            if not gap_coords_exist:
                                # gap_coord_tuple = (0, 0, 0, 0)
                                continue
                            tr2blocks[tr][gap_num] = gap_coord_tuple
                        # print(f'Incrementing current transcript number after {b_num}')
                        continue
                    if tr_start >= coord_tuple[1]:
                        # print(f'tr_start={tr_start}, block_end={coord_tuple[1]}')
                        ## somehow the transcript iterator moved downstream
                        ## to the current block; as long as transcripts are
                        ## sorted properly, there is no need to check the remaining
                        ## transcripts
                        ## if none of the previous blocks were added, then for
                        ## the next block, we can start iterating from the
                        ## current transcript
                        if not block_added:
                            # print(f'block_added={block_added}')
                            # print('Updating the transcript pointer after transcript {tr} and block {b_num}')
                            curr_tr = i
                        # if not tr2blocks[tr] and not first_block:
                        #     tr2blocks[tr][gap_num] = gap_coord_tuple
                        # print(f'Breaking at block {b_num} for transcript {tr}, positive strand')
                        break
                else:
                    ## logic is essentially the same except for now the iterator
                    ## goes upstream along the chromosome
                    if tr_start >= coord_tuple[1]:
                        curr_tr += 1
                        ## TODO: Shift abs_start/abs_stop?
                        if curr_tr < tr_num:
                            abs_stop = stops[curr_tr]
                        if gap_num and gap_num in tr2blocks[tr]:
                            tr2blocks[tr][str(b_num)] = coord_tuple
                        if not tr2blocks[tr]:
                            tr2blocks[tr][gap_num] = gap_coord_tuple
                        continue
                    if tr_end < coord_tuple[0]:
                        if not block_added:
                            curr_tr = i
                        # print(f'Breaking at block {b_num}, negative strand')
                        break
                # print(f'Transcript {tr} intersects block {b_num}')
                block_added = True
                if not tr2blocks[tr] and not first_block:
                    tr2blocks[tr][gap_num] = gap_coord_tuple
                # print(f'Adding block {coord_tuple} for transcript {tr}')
                tr2blocks[tr][str(b_num)] = coord_tuple
            # print(f'curr_tr={curr_tr}')
            if curr_tr >= len(names):
                line = EMPTY_STR
                break
            ## if this is the last block, there is no chain gap to follow
            if last_block:
                break
            ## otherwise, create a gap block
            if ref_strand:
                gap_num = str(b_num) + '_' + str(b_num+1)
            else:
                gap_num = str(b_num-1) + '_' + str(b_num)
            if ref_strand:
                r_end = r_start + dt
            else:
                r_end = r_start - dt
            if query_strand:
                q_end = q_start + dq
            else:
                q_end = q_start - dq
            gap_coord_tuple: Tuple[int] = (r_start, r_end, q_start, q_end)
            # print(f'Adding the gap {gap_num}, curr_tr={curr_tr}')
            ## and do the same trick
            gap_added = False
            for i in range(curr_tr, len(names)):
                tr = names[i]
                tr_start = starts[i]
                tr_end = stops[i]
                # print(f'tr={tr}, tr_start={tr_start}, tr_end={tr_end}, gap_coord_tuple={gap_coord_tuple}')
                if ref_strand:
                  ## do not add gaps lying downstream to the transcript
                    if tr_end < gap_coord_tuple[0]:
                        # print(f'Gap {gap_num} ({gap_coord_tuple}) starts downstream to tr {tr} ({tr_start}, {tr_end})')
                        # curr_tr += 1
                        continue
                        # pass
                    ## do not add gaps as a starting block unless we suspect a spanning chain
                    if not tr2blocks[tr] and tr_end > gap_coord_tuple[1]:
                        # print(f'Gap {gap_num} cannot be the first block for tr {tr}')
                        gap_added = True
                        continue
                        # break
                    ## do not gaps as last blocks unless we suspect a spanning chain
                    if tr_start < gap_coord_tuple[0] and tr_end < gap_coord_tuple[1]:
                        # print(f'Gap {gap_num} cannot be the last block for tr {tr}')
                        # continue
                        pass
                    ## if the gap lies upstream to the transcript, it does not
                    ## intersect neither it nor the following transcripts
                    if tr_start >= gap_coord_tuple[1]:
                        # print(f'Breaking at gap {gap_num}, positive strand')
                        if not gap_added:
                            curr_tr = i
                        break
                else:
                    if tr_start >= gap_coord_tuple[1]:
                        # curr_tr += 1
                        continue
                    if not tr2blocks[tr] and tr_start > gap_coord_tuple[0]: ## TODO: Check if the logic is inverted properly
                        # print(f'Gap {gap_num} cannot be the first block for tr {tr}')
                        gap_added = True
                        continue
                    if tr_end > gap_coord_tuple[1] and tr_start > gap_coord_tuple[0]:
                        # print(f'Gap {gap_num} cannot be the last block for tr {tr}')
                        continue
                    if tr_end < gap_coord_tuple[0]:
                        # print(f'Breaking at gap {gap_num}, negative strand')
                        if not gap_added:
                            curr_tr = i
                        break
                # print(f'Adding the gap {gap_num} for transcript {tr}')
                if not gap_added:
                    curr_tr = i
                gap_added = True
                # print(f'Adding gap {gap_num} ({gap_coord_tuple}) to tr {tr} ({tr_start}, {tr_end})')
                tr2blocks[tr][gap_num] = gap_coord_tuple
                if str(b_num) not in tr2blocks[tr]:
                    # print(f'{gap_coord_tuple}, {coord_tuple}')
                    # print(f'Adding the previous block {b_num} for transcript {tr} since the gap {gap_num} was also added')
                    tr2blocks[tr][str(b_num)] = coord_tuple
            r_start = r_end
            q_start = q_end
            # print(f'curr_tr={curr_tr}')
            if curr_tr >= len(names):
                # print(f'Breaking the outer loop at block {b_num}')
                if gap_num in tr2blocks[tr]:
                    tr2blocks[tr][str(b_num)] = coord_tuple
                line = EMPTY_STR
                break
            b_num += 1# if ref_strand else -1
            line = EMPTY_STR
            continue
        line += symbol
    # print(f'line={line}')
    # print(f'counter={counter}')
    # print('\n')
    # print(f'tr2blocks={tr2blocks}')
    # print('\n')
    return tr2blocks, chain_meta


# def split_at(string: str, sep: str):
#     return split_at_(string, sep)
