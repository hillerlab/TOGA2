import cython
cimport cython
from libc.stdio cimport FILE, fclose, fopen

cdef bytes CHAIN = b'chain'
cdef bytes NL = b'\n'
cdef bytes PLUS = b'+'
cdef bytes SPACE = b' '
cdef bytes TAB = b'\t'
cdef tuple FIELDS = (2, 3, 4, 5, 6, 12)
cdef int MIN_FIELD = 1
cdef int MAX_FIELD = 5

cdef int sort_coord_tuples(tuple coords):
    return coords[1]

cdef tuple get_chain_coords(bytes line):
    cdef int nf = 0
    cdef bytes field = b''
    cdef bytes sym = b''
    cdef str utf_field = ''
    cdef int size, start, stop
    cdef bint strand
    cdef str chrom, chain_id
    line += SPACE
    cdef int line_len = len(line)
    cdef int p = 6
    # while p < line_len:
        # print(line)
        # sym = line[p].to_bytes(1)
        # print(sym)
    for sym in line:
        if sym == SPACE:
            if nf in FIELDS:
                utf_field = field.decode('utf8')
                if nf == 2: ## reference chromosome name
                    chrom = utf_field
                if nf == 3: ## reference chromosome size
                    size = int(utf_field)
                if nf == 4: ## alignment strand in reference genome
                    strand = field == PLUS
                if nf == 5: ## strandwise start
                    if strand:
                        start = int(utf_field)
                    else:
                        stop = size - int(utf_field)
                if nf == 6: ## strandwise stop
                    if strand:
                        stop = int(utf_field)
                    else:
                        start = size - int(utf_field)
                if nf == 12: ## chain id
                    chain_id = utf_field
            field = b''
            utf_field = ''
            nf += 1
            continue
        field += sym
        p += 1
    return (chain_id, chrom, start, stop)



cpdef dict retrieve_chain_headers(str file):
    cdef bytes buff
    cdef bytes p
    cdef bytes bline = b''
    cdef tuple coords
    cdef str chain_id, chrom
    cdef int start, stop
    cdef dict output = {}
    # cdef FILE* f
    cdef bytes sym
    f = open(file, 'rb')
    # while fgetc(sym, 1, f)!=EOF:
    buff = f.read()
         # buff += sym
    f.close()
    for p in buff:
        if p == NL:
            if bline[:5] == CHAIN:
                coords = get_chain_coords(bline)
                chain_id = coords[0]
                chrom = coords[1]
                start = coords[2]
                stop = coords[3]
                if chrom not in output:
                    output[chrom] = []
                output[chrom].append((chain_id, start, stop))
            bline = b''
            continue
        bline += p
    for chrom in output:
        output[chrom] = sorted(output[chrom], key=sort_coord_tuples)
    return output


cdef tuple extract_bed_data(bytes line):
    cdef int nf = 1
    cdef str chrom, tr
    cdef int start, stop
    cdef bytes field = b''
    cdef str utf_field
    cdef bytes sym = b''
    for sym in line:
        if sym == TAB:
            utf_field = field.decode('utf8')
            if nf == 1:
                chrom = utf_field
            if nf == 4:
                tr = utf_field
            if nf == 7:
                start = int(utf_field)
            if nf == 8:
                stop = int(utf_field)
            nf += 1
            if nf > 8:
                break
            field = b''
            continue
        field += sym
    return (tr, chrom, start, stop)


cpdef dict get_bed_coords(str file):
    cdef bytes buff
    cdef bytes p
    cdef bytes bline = b''
    cdef tuple bed_data
    cdef dict output = {}
    # cdef FILE* f
    cdef bytes sym
    f = open(file, 'rb')
    # while fgetc(sym, 1, f)!=NULL:
    buff = f.read()
         # buff += sym
    f.close()
    for p in buff:
        if p == NL:
            bed_data = extract_bed_data(bline)
            tr = bed_data[0]
            chrom = bed_data[1]
            start = bed_data[2]
            stop = bed_data[3]
            if chrom not in output:
                output[chrom] = []
            output[chrom].append((tr, start, stop))
            bline = b''
            continue
        bline += p
    for chrom in output:
        output[chrom] = sorted(output[chrom], key=sort_coord_tuples)
    return output


cdef int intersection(int s1, int e1, int s2, int e2):
    cdef int min_end
    cdef int max_start
    if e1 > e2:
        min_end = e2
    else:
        min_end = e1
    if s1 > s2:
        max_start = s1
    else:
        max_start = s2
    return min_end - max_start
