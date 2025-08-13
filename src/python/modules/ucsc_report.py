"""
A collection of constants and executables related to UCSC BigBed report formatting
"""

from cesar_wrapper_constants import (
    FS_INDELS, COMPENSATION, DEL_MISS, DEL_EXON, MISS_EXON, SSM
)
from cesar_wrapper_executables import Mutation
from collections import defaultdict
from dataclasses import dataclass ## TODO: HIGHLY REDUNDANT
from shared import parts
from typing import Any, Dict, Iterable, List, Optional, Set, TextIO, Tuple 

__author__ = 'Yury V. Malovichko'
__credits__ = ('Bogdan Kirilenko', 'Björn Langer', 'Michael Hiller')
__year__ = '2024'

## define the constants
ABBR2STATUS: Dict[str, str] = {
    'FI': 'Fully intact', 
    'I': 'Intact', 
    'PI': 'Partially intact', 
    'M': 'Missing', 
    'PM': 'Partially missing', 
    'L': 'Lost', 
    'UL': 'Uncertain Loss', 
    'PG': 'Paralogous projection',
    'PP': 'Processed pseudogene'
}
MUT_SLOTS_TO_IGNORE: Tuple[str] = ('projection', 'chrom', 'start', 'stop')
OUT_OF_CHAIN_PLACEHOLDER: float = 0.0
BR: str = '<BR>'
BOX_START: str = '<TT>'
BOX_END: str = '</TT>'
BOLD_START: str = '<B>'
BOLD_END: str = '</B>'
HEAD_START: str = '<H{}>'
HEAD_END: str = '</H{}>'
SPACE: str = "&nbsp;"
REF_LINK_PLACEHOLDER: str = '<A HREF="javascript:;"></A>'
PLACE_HOLDER_EXON_MID: str = SPACE * 5
ALN_HEADER: str = 'Sequence alignment between reference and query exon:'
PROTEIN_HEADER: str = 'Predicted protein sequence ({} amino acids, {} nucleotides):<BR>'
PROTEIN_ALN_LINE: str = f'ref:{SPACE}{{}}{BR}{PLACE_HOLDER_EXON_MID}{{}}{BR}que:{SPACE}{{}}{BR}{BR}'
EXON_INTRO: str = 'Sequence alignment between reference and query exon:'
EXON_NUMBER: str = f'{HEAD_START}Exon number: {{}}{HEAD_END}'
EXON_REGION: str = f'{BOLD_START}Exon region:{BOLD_END} {{}}:{{}}-{{}}'
NUC_PERC_ID: str = f'{BOLD_START}Nucleotide identity:{BOLD_END} {{}}'
BLOSUM: str = f'{BOLD_START}BLOSUM:{BOLD_END} {{}}'
INTER_GAP: str = f'{BOLD_START}Intersects assembly gaps:{BOLD_END} {{}}'
ALN_CLASS: str = f'{BOLD_START}Exon alignment class:{BOLD_END} {{}}'
EXP_REGION: str = f'{BOLD_START}Detected within expected region ({{}}:{{}}-{{}}):{BOLD_END} {{}}'
EXON_ENTRY: str = f'{{}}{BR}{BR}{{}}{BR}{{}}{BR}'

## define the constants ## TODO: Check which of these are still needed
# Mutation classes
EX_DEL = "Deleted exon"
EX_MIS = "Missing exon"
START_MIS = "START_MISSED"
FS_DEL = "FS_DEL"
FS_INS = "FS_INS"
BIG_DEL = "BIG_DEL"
BIG_INS = "BIG_INS"
INTRON_GAIN: str = 'INTRON_GAIN'
INDELS: Tuple[str] = (FS_DEL, FS_INS, BIG_DEL, BIG_INS, INTRON_GAIN)
BIG_INDELS: Tuple[str] = (BIG_DEL, BIG_INS, INTRON_GAIN)
INS: Tuple = (FS_INS, BIG_INS)
DELS: Tuple[str] = (FS_DEL, BIG_DEL)
STOP: str = 'STOP'
# (ag)acceptor-EXON-donor(gt)
SSM_D: str = "SSMD"  # Donor, right, GT,GC
SSM_A: str = "SSMA"  # Acceptor, left, AG
COMP: str = 'COMPENSATION'

## define the template constants
TEMPLATE_PATH_1 = "svg_template.txt"
TEMPLATE_PATH_2 = "supply/svg_template.txt"
DEFAULT_TEMPLATE: str = (
    '<svg version="1.1" '
    'xmlns="http://www.w3.org/2000/svg" '
    'xmlns:xlink="http://www.w3.org/1999/xlink" '
    'width="{}" height="{}" viewBox = "0 0 {} {}" '
    'onclick="mouseclick(evt)">{}'
    '</svg>'
)

OPACITY = 100
HORIZONTAL: str = "horizontal"
MUT_LINE_FIELDS: int = 8
MUT_LINE_FIELDS_SP: int = MUT_LINE_FIELDS ## TODO: Redundant???

BLACK: str = "#121212"  # almost black
MISS_SEQ_COLOR: str = "#878787"  # grey
MASKED_MUT_COLOR: str = "#878787"  # grey
INACT_MUT_COLOR: str = "#cf232b"  # dark red
FP_EXON_DEL_COLOR: str = "#0c7bdc"  # blue-ish
BACKGROUND_COLOR: str = "#87bcbc"  # cyan-ish
STOP_CODON_COLOR: str = "#121212"  # almost black

## define size constants (all values given in pixels if not stated otherwise)
MAX_VERTICAL_SPACE: int = 100
HALF_EXON_HEIGHT: int = 15
HALF_UTR_HEIGHT: int = 5  ## legacy constant
UTR_WIDTH: int = 0  ## legacy constant
INTER_FRAGMENT_SPACE: float = 38.0
VERTICAL_OFFSET: float = 30.0
HORIZONTAL_OFFSET: float = 15.0
CHAIN_LABEL_OFFSET: float = 15.0
TRANSCRIPT_NAME_OFFSET: float = 10.0

EXON_BASE_SIZE: int = 0.8
INTRON_BASE_SIZE = 0.02
MAX_INTRON_SIZE: float = 2000.0  ## bases
MIN_INTRON_SIZE: float = 600.0  # bases
INTACT_INTRON_SIZE: int = MIN_INTRON_SIZE * INTRON_BASE_SIZE

GAP_WIDTH: int = 0  # pixels
HALF_GAP_HEIGHT: int = 0  # pixels

MIN_ARROW_SIZE: int = 10  # pixels
MAX_ARROW_SIZE: int = 20  # pixels
ARROW_SIZE: int = 1  # pixels per base

PRINT_EXON_SCHEME: bool = False
ROUND_EDGES: bool = False
INTRON_STYLE: str = 'stroke:#999; stroke-width:3;'
INTRON_OBJECT: str = f'<line x1="{{}}" y1="{{}}" x2="{{}}" y2="{{}}" style="{INTRON_STYLE};" />'
EXON_ANC_STYLE: str = "stroke-width:3;"
EXON_NON_ANC_STYLE: str = "stroke: black; stroke-width:3; stroke-dasharray: 5,5;"
EXON_STYLE_PLATE: str = 'fill:{};fill-opacity:1.00'
EXON_OBJECT: str = '<rect class="anc_exon" x="{}" y="{}" width="{}" height="{}" style="{}" />'
INSERT_STYLE: str = 'style="fill:{0}; stroke-opacity:1; fill-opacity:1"'
INSERT_TEMPLATE: str = '  <polygon points="{},{} {},{} {},{}" {}/>'
DEL_STYLE: str = 'style="stroke:{0}; stroke-width:{1}; stroke-opacity:1"'
MISSSEQ_STYLE: str = 'style="stroke:{0}; stroke-width:{1}; stroke-opacity:1"'
STOP_CODON_STYLE: str = 'style="stroke:{0};stroke-width:3;"'
COMP_INDEL_STYLE: str = 'style="fill:none;stroke-width:1;stroke:green;"'
COMPENSATION_ARC_COORDS: str = 'M {} {} C {} {} {} {} {} {}'
COMPENSATION_TEMPLATE: str = '  <path d="{}" {} />'
TEXT_HIGHLIGHT_STYLE: str = f'style="fill:{INACT_MUT_COLOR};"'
TEXT_STYLE: str = 'font-size:{}px;fill:{};'
TEXT_TEMPLATE: str = (
    '<text style="fill:{};" ><tspan x="{}" y="{}" style="{}" >{}</tspan></text>'
)
FONTFAMILY: str = "Courier New"
FONT_ASPECT_RATIO: float = 0.5 ## actually, it's 0.39
CHAIN_ID_FONTSIZE: int = 18
STOP_LABEL_FONTSIZE: int = 15
SS_LABEL_FONTSIZE: int = 15
MO_FONTSIZE: int = 15
EX_NUM_FONTSIZE: int = 18
TRANSCRIPT_NAME_FONTSIZE: int = 20
STACKING_THRESHOLD: int = 35  # pixels
J_STYLE: bool = False
REF_LINK_PLACEHOLDER: str = '<A HREF="javascript:;"></A>'


def _bold(string: str) -> str:
    """Formats input string as an HTML bold text"""
    return f'{BOLD_START}{string}{BOLD_END}'

def _td(string: str) -> str:
    """Formats input string as an HTML table cell"""
    return f'<td>{string}</td>'

def format_fasta_as_aln(
        seq1: str, 
        seq2: str, 
        w: Optional[int] = 80,
        protein: Optional[bool] = False, 
        bp_num: Optional[int] = None
    ) -> str:
    """Formats two strings as alignment pseudographics"""
    out_line: str = BOX_START
    header_line: str = ''
    if protein:
        aa_num: int = len(seq2.replace('-', ''))
        bp_num = bp_num if bp_num is not None else 'X'
        header_line = PROTEIN_HEADER.format(aa_num, bp_num)
    out_line += header_line
    chunked_aln: List[List[Tuple[str]]] = parts(list(zip(seq1, seq2)), w)
    for portion in chunked_aln:
        s1: str = ''.join((x[0].upper() for x in portion))
        middle: str = ''.join(
            (SPACE if x[0].upper() != x[1].upper() or x[0] == '-' else '|' for x in portion)
        )
        s2: str = ''.join((x[1].upper() for x in portion))
        out_line += PROTEIN_ALN_LINE.format(s1, middle, s2)
    out_line += BOX_END
    return out_line


def exon_aln_header(
    exon: int,
    chrom: str,
    start: int,
    end: int,
    exp_start: int,
    exp_end: int,
    found_in_exp: bool,
    nuc_id: float,
    blosum: float,
    gap_intersection: bool,
    aln_class: str
) -> str:
    """
    """
    number: str = EXON_NUMBER.format(5, exon, 5)
    region: str = EXON_REGION.format(chrom, start, end)
    nuc_id = round(nuc_id, 2)
    blosum = round(blosum, 2)
    id_line: str = f'{NUC_PERC_ID.format(nuc_id)} | {BLOSUM.format(blosum)}'
    gap_line: str = INTER_GAP.format('YES' if gap_intersection else 'NO')
    aln_class_line: str = ALN_CLASS.format(aln_class)
    exp_reg_line: str = EXP_REGION.format(chrom, exp_start, exp_end, 'YES' if found_in_exp else 'NO')
    header: str = BR.join([number, region, id_line, gap_line, aln_class_line, exp_reg_line])
    return header


def mutation_table(muts: List[Mutation]) -> str:
    """Given a list of Mutation objects, prepares an HTML data table"""
    out_line: str = ''
    for mut in muts:
        exon: str = _td(mut.exon)
        codon: str = _td(mut.ref_codon)
        mut_type: str = _td(
            mut.mutation_class if mut.description != '-' else '-'
        )
        descr: str = _td(
            mut.description if mut.description != '-' else mut.mutation_class
        )
        # masked: str = _td('YES' if mut.is_masked else 'NO')
        is_inactivating: str = _td('NO' if mut.is_masked else 'YES')
        mut_id: str = _td(mut.mutation_id)
        reason: str = _td(mut.masking_reason)
        out_line += f'<tr>{exon}{codon}{mut_type}{descr}{is_inactivating}{reason}{mut_id}</tr>'
    return out_line

def exon_aln_entry(seq1: str, seq2: str, header: str) -> str:
    """Generates header line for an exon alignment entry"""
    aln: str = format_fasta_as_aln(seq1, seq2)
    return EXON_ENTRY.format(header, _bold(EXON_INTRO), aln)


def get_chain_features(
    file: TextIO, str, proj_list: Iterable[str] = []
) -> Dict[str, Tuple[str]]:
    """
    Extracts feature values from the chain_results_df.tsv TOGA output file,
    returns the {projection: (feature_tuple)} dictionary
    """
    proj2features: Dict[str, Tuple[str]] = {}
    for line in TextIO:
        data: List[str] = line.rstrip().split('\t')
        if len(data) != 16:
            raise ValueError(
                'Classification feature file differs from the expected format'
            )
        if data[0] == 'transcript':
            continue
        trans: str = data[0]
        chain: str = data[2]
        proj: str= f'{trans}#{chain}'
        if proj_list and proj not in proj_list:
            continue
        synt_: str = data[3]
        gl_exo_: str = data[5]
        loc_exon_: str = data[8]
        exon_cover_: str = data[9]
        intron_cover_: str = data[10]
        exon_fract_: str = data[13]
        intron_fract_: str = data[14]
        flank_cov_: str = data[15]

        exon_cov: str = (
            str(float(exon_cover_) / float(exon_fract_))
            if float(exon_fract_) != 0
            else "0"
        )
        intron_cov: str = (
            str(float(intron_cover_) / float(intron_fract_))
            if float(intron_fract_) != 0
            else "0"
        )
        proj2features[proj] = (
            synt_, flank_cov_, gl_exo_, loc_exon_, exon_cov, intron_cov
        )
    return proj2features


@dataclass
class IntronDash:
    __slots__ = ['x1', 'x2', 'y']
    x1: float
    x2: float
    y: float

    def line(self) -> str:
        """Returns inton SVG representation"""
        return INTRON_OBJECT.format(self.x1, self.y, self.x2, self.y)


@dataclass
class ExonBox:
    __slots__ = ['num', 'color', 'x', 'y', 'width', 'height']
    num: int
    color: str
    x: float
    y: float
    width: float
    height: float

    def line(self) -> str:
        """Returns exon SVG representation"""
        style: str = EXON_STYLE_PLATE.format(self.color)
        return EXON_OBJECT.format(self.x, self.y, self.width, self.height, style)


@dataclass
class TextStack:
    __slots__ = ['x', 'y', 'label', 'size', 'color']
    x: float
    y: float
    label: str
    size: float
    color: str

    def line(self) -> str:
        """Returns SVG string representation of the text field"""
        style: str = TEXT_STYLE.format(self.size, self.color)
        return TEXT_TEMPLATE.format(self.color, self.x, self.y, style, self.label)


@dataclass
class Insertion:
    __slots__ = ['length', 'x', 'y', 'masked']
    length: int
    x: float
    y: float
    masked: bool

    def line(self) -> str:
        """Returns a superscript inverted triangle marking a point mutation"""
        height: float = min(max(self.length * ARROW_SIZE, MIN_ARROW_SIZE), MAX_ARROW_SIZE)
        width: float = height / 2
        ## get the coordinates of the right left angle
        p2x_: float = self.x + width / 2
        p2y_: float = self.y - height
        ## and for upper left angle
        p3x_: float = self.x - width / 2
        p3y_: float = self.y - height
        color: str = MASKED_MUT_COLOR if self.masked else INACT_MUT_COLOR
        style: str = INSERT_STYLE.format(color)
        return INSERT_TEMPLATE.format(self.x, self.y, p2x_, p2y_, p3x_, p3y_, style)


@dataclass
class BlockMutation:
    __slots__ = ['x', 'y', 'width', 'height', 'color']
    x: float
    y: float
    width: float
    height: float
    color: bool

    def line(self) -> str:
        """Returns SVG string representation of the text field"""
        style: str = EXON_STYLE_PLATE.format(self.color)
        return EXON_OBJECT.format(self.x, self.y, self.width, self.height, style)


class ProjectionPlotter:
    # __slots__ = [
    #     'ref_tr', 'projection', 'mutations', 'fragm2ex',
    #     'max_width', 'max_height', 'exon2muts', 'compensation2fs',
    #     'mut_masked', 'mut_label', 'missing_exons', 'deleted_exons',
    #     'last_exon', 'tr_label'
    # ]
    def __init__(
        self,
        tr: str,
        ref_tr: Dict[int, int],
        # projection: Iterable[ProjectionFeatures],
        mutations: Iterable[Mutation],
        exon2chain: Dict[int, str]
    ) -> None:
        self.tr_label: str = tr
        self.ref_tr: Dict[int, int] = ref_tr
        # self.projection: Iterable[ProjectionFeatures] = projection
        self.mutations: Iterable[Mutation] = mutations
        self.fragm2ex: Dict[int, List[int]] = defaultdict(list)
        for x,y in exon2chain.items():
            self.fragm2ex[y].append(x)
        self.max_width: int = 0
        self.max_height: int = 0
        self.exon2muts: Dict[int, List[str]] = defaultdict(list)
        self.exon2comps: Dict[int, List[str]] = defaultdict(list)
        self.mut2position: Dict[str, int] = {}
        self.mut_masked: Dict[str, bool] = {}
        self.mut_label: Dict[str, str] = {}
        self.mut2comp: Dict[str, str] = {}
        self.comp2mut: Dict[str, Set[str]] = defaultdict(list)
        self.mut2exon: Dict[str, int] = {}
        self.missing_exons: Set[int] = set()
        self.critical_deletions: Set[int] = set()
        self.safe_deletions: Set[int] = set()
        self.last_exon: int = max(ref_tr)
        # chain_postfix: str = ','.join(map(str, sorted(self.fragm2ex, key=lambda x: min(self.fragm2ex[x]))))
        # self.tr_label: str = f'{self.ref_tr.transcript}#{chain_postfix}'

        self.exon_boxes: Dict[int, ExonBox] = {}
        self.regular_intron_dashes: Dict[int, IntronDash] = {}
        self.terminal_intron_dashes: Dict[int, IntronDash] = {}
        self.mut_objects: Dict[str, Any] = {}

        self._run()

    def _run(self) -> None:
        """
        """
        self._compute_plot_width()
        self._compute_plot_height()
        self._prepare_mutations()
        self._generate_filebuffer()
        # svg = self.plot_svg()
        # print(svg)

    def _compute_plot_width(self) -> None: ## DONE!
        """
        Computes future picture width. In the simplest case, calculates the width
        needed to accommodate for the whole transcript; if projection comprises of
        multiple fragments, return the maximum width needed across all fragments.
        """
        max_width: int = 0
        ## for fragmented projections, an additional left-side offset is needed
        ## to accommodate for chain labels
        add_ids: bool = len(self.fragm2ex) > 1
        longest_chain_label: int = max(len(str(x)) for x in self.fragm2ex)
        ## for each fragment, calculate the width of the respective subplot
        for fragm, exons in self.fragm2ex.items():
            gene_width: int = 2 * HORIZONTAL_OFFSET + (
                CHAIN_ID_FONTSIZE * longest_chain_label * FONT_ASPECT_RATIO + 5.0
            ) * add_ids
            for exon in exons:
                ## add the exon width
                exon_size: int = self.ref_tr[exon] * EXON_BASE_SIZE
                gene_width += exon_size
                ## if this is the fragment's first exon and it's not exon 1,
                ## add half the length of the previous intron
                if exon == min(exons) and exon != 1:
                    # this_exon_start: int = self.ref_tr.exons[exon].__getattribute__(
                    #     'start' if self.strand else 'stop'
                    # )
                    # prev_exon_stop: int = self.ref_tr.exons[exon-1].__getattribute__(
                    #     'stop' if self.strand else 'start'
                    # )
                    # prev_intron_len: int = (1 if self.strand else -1) * (this_exon_start - prev_exon_stop)
                    intron_size: float = self._optimal_intron(exon)
                    gene_width += intron_size
                ## otherwise, add the next intron's length
                ## if this is the fragment's last exon and it's not the terminal exon,
                ## add half the length; otherwise, add the full intron
                ## for the C-terminal exon, add nothing
                if exon != self.last_exon:
                    # this_exon_stop: int = self.ref_tr.exons[exon].__getattribute__(
                    #     'stop' if self.strand else 'start'
                    # )
                    # next_exon_start: int = self.ref_tr.exons[exon+1].__getattribute__(
                    #     'start' if self.strand else 'stop'
                    # )
                    # next_intron_len: int = (1 if self.strand else -1) * (next_exon_start - this_exon_stop)
                    intron_size: float = self._optimal_intron(exon)
                    ## if this is the last exon of the segment yet not the terminal
                    ## exon overall, add half the next intron; otherwise, add
                    ## the full intron
                    gene_width += intron_size
            max_width = max(max_width, gene_width)
        text_label_width: int = len(self.tr_label) * TRANSCRIPT_NAME_FONTSIZE
        self.max_width = max(max_width, text_label_width)

    def _compute_plot_height(self) -> None: ## TODO: MODIFY AS NEEDED
        ## TODO: Revise once done with generate_filebuffer() update
        """
        Compute figure height based on the number of fragments in the projection
        """
        frag_num: int = len(self.fragm2ex)
        ## sum the cumulative heights of all fragments
        picture_height: int =  frag_num * HALF_EXON_HEIGHT * 2
        ## and add the space between fragments
        picture_height += (frag_num - 1) * INTER_FRAGMENT_SPACE
        ## and a little bit of offset from above and below
        picture_height += 2 * VERTICAL_OFFSET
        ## and some more for the exon numeration line
        picture_height += EX_NUM_FONTSIZE * frag_num
        ## and, finally, add the space for the transcript name
        picture_height += TRANSCRIPT_NAME_FONTSIZE + TRANSCRIPT_NAME_OFFSET

        self.max_height = picture_height

    def _prepare_mutations(self) -> None:
        """
        Retrieves mutations which will be added to the HTML report plot.
        These contain:
        * frameshifting indels;
        * nonsense mutations;
        * missing and deleted exons;
        * frameshift compensations (as a decorator over plotted compensations)
        """
        for mut in self.mutations:
            ex: int = int(mut.exon)
            missing: bool = ex in self.missing_exons
            deleted: bool = ex in self.safe_deletions or ex in self.critical_deletions
            if missing or deleted:
                continue
            mut_class: str = mut.mutation_class
            mut_meta: str = mut.description
            mut_id: str = mut.mutation_id
            masked: bool = mut.is_masked #mut.is_masked == 'MASKED' #
            ## if it's a frameshift or a stop codon, the following data are collected:
            ## 1) location in the reference (TODO: currently alignment is being tracked)
            ## 2) masking status
            ## 3) mutation label
            # print(f'{mut_class=}, {mut_meta=}, {DEL_MISS=}, {mut_class in DEL_MISS=}')
            if mut_class in FS_INDELS or mut_class == STOP:
                rel_position: int = int(mut.ref_codon) * 3 - 2 ## TODO: Recalculate codons into reference codon numbers
                self.mut2position[mut_id] = rel_position
                self.mut_masked[mut_id] = masked
                if mut_class in FS_INDELS:
                    self.mut_label[mut_id] = (
                        mut_meta if int(mut_meta) < 0 else f'+{mut_meta}'
                    )
                else:
                    self.mut_label[mut_id] = mut_meta.split('->')[1]
            ## or, if it's a compensation, add it to the compensation slot
            ## and keep the IDs of respective frameshifts
            elif mut_class == COMPENSATION:
                fs_start, fs_stop = map(int, mut_meta.split('_')[1].split('-'))
                for fs in range(fs_start, fs_stop + 1):
                    self.mut2comp[f'FS_{fs}'] = mut_id
                    self.comp2mut[mut_id].append(f'FS_{fs}')
                self.exon2comps[ex].append(mut_id)
                self.mut2exon[mut_id] = ex
                continue
            ## missing and deleted exons are inferred from respective mutations
            elif mut_class in DEL_MISS:#== '-' and mut_meta in DEL_MISS:
                if mut_class == MISS_EXON:
                    self.missing_exons.add(ex)
                elif mut_class == DEL_EXON:
                    if masked:
                        self.safe_deletions.add(ex)
                    else:
                        self.critical_deletions.add(ex)
                else:
                    raise Exception('Missing/Deleted exon meta is corrupted for %i' % ex)
                for m_ in self.exon2muts[ex]:
                    if m_ in self.mut2exon:
                        del self.mut2exon[m_]
                    if m_ in self.mut_label:
                        del self.mut_label[m_]
                    if m_ in self.mut_masked:
                        del self.mut_masked[m_]
                self.exon2muts[ex].clear()
            ## FAFO how to handle splice shift mutations
            elif mut_class in SSM:
                sites: List[str] = mut_meta.split('->')
                query_site: str = '--' if len(sites) < 2 else sites[1]
                self.mut_label[mut_id] = query_site
                self.mut_masked[mut_id] = masked
            elif mut_class in BIG_INDELS:
                continue
            ## other mutations are not currently added to the plot
            else:
                continue
            self.exon2muts[ex].append(mut_id)
            self.mut2exon[mut_id] = ex

    def _generate_filebuffer(self) -> str:
        """
        Generates the file buffer, i.e. lines standing for actual graphic
        representation in the resulting SVG file
        """
        self.buffer_lines: List[str] = [] ## TODO: Should the buffer be kep as a string by default?
        fragms: List[int] = sorted(
            self.fragm2ex.keys(), key=lambda x: min(self.fragm2ex[x])
        )
        add_labels: bool = len(self.fragm2ex) > 1
        if add_labels:
            longest_chain_label: int = max(len(str(x)) for x in self.fragm2ex)
            chain_id_offset: float = (
                CHAIN_ID_FONTSIZE * longest_chain_label * FONT_ASPECT_RATIO + 5.0
            )
        else:
            chain_id_offset: float = 0.0
        for i, fragm in enumerate(fragms, start=1):
            exons: Set[int] = self.fragm2ex[fragm]
            first_ex: int = min(exons)
            last_ex: int = max(exons)
            ## estimate the x coordinate
            ex_y: int = (
                VERTICAL_OFFSET + HALF_EXON_HEIGHT * 2 * (i - 1) +
                INTER_FRAGMENT_SPACE * (i - 1)
            )
            in_y: int = ex_y + HALF_EXON_HEIGHT
            ## set the x start coordinate
            fragm_x: float = HORIZONTAL_OFFSET
            # print(f'{self.fragm2ex=}, {add_labels=}')
            # print(f'{self.safe_deletions=}, {self.critical_deletions=}')
            if add_labels:
                chain_label_x: float = fragm_x #+ CHAIN_LABEL_OFFSET / 2
                chain_label_y: float = in_y + CHAIN_ID_FONTSIZE // 2
                chain_id_label: TextStack = TextStack(
                    chain_label_x, chain_label_y, fragm, CHAIN_ID_FONTSIZE, BLACK
                )
                self.buffer_lines.append(chain_id_label.line())
                fragm_x += chain_id_offset + CHAIN_LABEL_OFFSET
            ## if the fragment starts with an exon other than 1,
            ## add half the previous intron at the beginning of the subplot
            if first_ex > 1:
                # this_exon_start: int = self.ref_tr.exons[first_ex].__getattribute__(
                #     'start' if self.strand else 'stop'
                # )
                # prev_exon_stop: int = self.ref_tr.exons[first_ex-1].__getattribute__(
                #     'stop' if self.strand else 'start'
                # )
                # prev_intron_len: int = (1 if self.strand else -1) * (
                #     this_exon_start - prev_exon_stop
                # )
                intron_size: float = self._optimal_intron(first_ex - 1)
                intron_start: int = fragm_x
                intron_stop: int = fragm_x + intron_size // 2
                _intron: IntronDash = IntronDash(intron_start, intron_stop, in_y)
                self.terminal_intron_dashes[first_ex] = _intron
                self.buffer_lines.append(_intron.line())
                fragm_x += intron_size // 2
            ## then, plot all the exons in this fragment
            for ex in sorted(exons):
                exon_width: float = self.ref_tr[ex] * EXON_BASE_SIZE
                exon_height: int = HALF_EXON_HEIGHT * 2
                ex_color: str = self._exon_color(ex)
                _exon: ExonBox = ExonBox(ex, ex_color, fragm_x, ex_y, exon_width, exon_height)
                self.exon_boxes[ex] = _exon
                exon_line: str = _exon.line()
                self.buffer_lines.append(exon_line)
                ex_mid: float = fragm_x + exon_width / 2
                ex_num_y: float = ex_y + (HALF_EXON_HEIGHT * 2 + EX_NUM_FONTSIZE)
                exon_label: TextStack = TextStack(
                    ex_mid, ex_num_y, ex, EX_NUM_FONTSIZE, BLACK
                )
                self.buffer_lines.append(exon_label.line())
                ## for each exon, plot all the respective mutations
                for mut in self.exon2muts[ex]:
                    if mut not in self.mut_label:
                        continue
                    label: str = self.mut_label[mut]
                    masked: bool = self.mut_masked[mut]
                    if mut.split('_')[0] in SSM:
                        if mut.split('_')[0] == SSM_A:
                            ss_x: float = fragm_x - len(label) * SS_LABEL_FONTSIZE * FONT_ASPECT_RATIO
                            ss_y: float = ex_y + HALF_EXON_HEIGHT + SS_LABEL_FONTSIZE
                        else:
                            ss_x: float = fragm_x + exon_width
                            ss_y: float = ex_y + HALF_EXON_HEIGHT
                        ss_col: str = MASKED_MUT_COLOR if masked else INACT_MUT_COLOR
                        ss_text: TextStack = TextStack(
                            ss_x, ss_y, label, SS_LABEL_FONTSIZE, ss_col
                        )
                        self.buffer_lines.append(ss_text.line())
                    elif 'FS' in mut:
                        label: int = int(label)
                        pos: int = self.mut2position[mut]
                        indel_x: float = fragm_x + (pos - self.len_till_exon(ex)) * EXON_BASE_SIZE
                        if label > 0:
                            ## create a short insertion object
                            ins_color: str = MASKED_MUT_COLOR if masked else INACT_MUT_COLOR
                            ins: Insertion = Insertion(abs(label), indel_x, ex_y, ins_color)
                            self.mut_objects[mut] = ins
                            self.buffer_lines.append(ins.line())
                            ins_label_text: str = f'+{label}'
                            ins_label: TextStack = TextStack(
                                indel_x - abs(label) * EXON_BASE_SIZE,
                                ex_y - HALF_UTR_HEIGHT - min(max(abs(label), MIN_ARROW_SIZE), MAX_ARROW_SIZE),
                                ins_label_text,
                                MO_FONTSIZE,
                                ins_color
                            )
                            self.buffer_lines.append(ins_label.line())
                        else:
                            del_width: float = abs(label) * EXON_BASE_SIZE
                            del_height: float = exon_height
                            del_color: str = MASKED_MUT_COLOR if masked else INACT_MUT_COLOR
                            del_: BlockMutation = BlockMutation(
                                indel_x, ex_y, del_width, del_height, del_color
                            )
                            self.mut_objects[mut] = del_
                            self.buffer_lines.append(del_.line())
                            del_label_text: str = str(label)
                            del_label: TextStack = TextStack(
                                indel_x - abs(label) * EXON_BASE_SIZE,
                                ex_y - HALF_UTR_HEIGHT,
                                del_label_text,
                                MO_FONTSIZE,
                                del_color
                            )
                            self.buffer_lines.append(del_label.line())
                    elif 'STOP' in mut and 'LOSS' not in mut:
                        pos: int = self.mut2position[mut]
                        stop_x: float = fragm_x + (pos - self.len_till_exon(ex)) * EXON_BASE_SIZE
                        stop_width: float = 3 * EXON_BASE_SIZE
                        stop_height: float = exon_height
                        stop_color: str = MASKED_MUT_COLOR if masked else STOP_CODON_COLOR
                        stop_mut: BlockMutation = BlockMutation(
                            stop_x, ex_y, stop_width, stop_height, stop_color
                        )
                        self.mut_objects[mut] = stop_mut
                        self.buffer_lines.append(stop_mut.line())
                        stop_label_text: str = str(label)
                        stop_label: TextStack = TextStack(
                            stop_x - 3 * EXON_BASE_SIZE,
                            ex_y - HALF_UTR_HEIGHT,
                            stop_label_text,
                            STOP_LABEL_FONTSIZE,
                            stop_color
                        )
                        self.buffer_lines.append(stop_label.line())
                ## for each exon except for the last one, add the trailing intron;
                ## if this is the last exon for the fragment but not for the
                ## projection, add half the exon instead
                fragm_x += exon_width
                if ex != self.last_exon:
                    # this_exon_stop: int = self.ref_tr.exons[ex].__getattribute__(
                    #     'stop' if self.strand else 'start'
                    # )
                    # next_exon_start: int = self.ref_tr.exons[ex+1].__getattribute__(
                    #     'start' if self.strand else 'stop'
                    # )
                    # next_intron_len: int = (1 if self.strand else -1) * (next_exon_start - this_exon_stop)
                    intron_size: float = self._optimal_intron(ex)
                    if ex == last_ex:
                        intron_size //= 2
                    intron_start: int = fragm_x
                    intron_stop: int = fragm_x + intron_size
                    _intron = IntronDash(intron_start, intron_stop, in_y)
                    self.regular_intron_dashes[ex] = _intron
                    intron_line: str = _intron.line()
                    self.buffer_lines.append(intron_line)
                    fragm_x += intron_size
        for ex in sorted(self.ref_tr):
            for comp in self.exon2comps[ex]:
                self.buffer_lines.extend(self._compensation_lines(comp))
                    # for comp_line in self._compensation_lines(comp):
                    #     self.buffer_lines.append(comp_line)

        tr_label_x: float = HORIZONTAL_OFFSET
        tr_label_y: float = self.max_height - VERTICAL_OFFSET
        transcript_label: TextStack = TextStack(
            tr_label_x, tr_label_y, self.tr_label, TRANSCRIPT_NAME_FONTSIZE, BLACK
        )
        self.buffer_lines.append(transcript_label.line())
        # print(f'{self.max_height=}, {self.max_width=}')
        # print('\n'.join(self.buffer_lines))


    def plot_svg(self) -> str:
        """
        Returns a string representation of the SVG plot
        """
        buffer: str = '\n'.join(self.buffer_lines)
        plot: str = DEFAULT_TEMPLATE.format(
            self.max_width, self.max_height,
            self.max_width, self.max_height,
            buffer
        )
        return plot

    def _optimal_intron(self, exon: int) -> float:
        """Returns the optimal plotted intron length"""
        affected_intron_size: float = 3.5 * SS_LABEL_FONTSIZE * FONT_ASPECT_RATIO
        for mut in self.exon2muts[exon]:
            if mut.split('_')[0] == SSM_D:
                return affected_intron_size
        if exon + 1 in self.exon2muts:
            for mut in self.exon2muts[exon + 1]:
                if mut.split('_')[0] == SSM_A:
                    return affected_intron_size
        return INTACT_INTRON_SIZE
        # return max(MIN_INTRON_SIZE, min(MAX_INTRON_SIZE, intron_length)) * INTRON_BASE_SIZE

    def _exon_color(self, exon: int) -> str:
        """
        Returns exon box color
        """
        color: str = ''
        if exon in self.missing_exons:
            color = MISS_SEQ_COLOR
        elif exon in self.critical_deletions:
            color = INACT_MUT_COLOR
        elif exon in self.safe_deletions:
            color = FP_EXON_DEL_COLOR
        else:
            color = BACKGROUND_COLOR
        return color

    def len_till_exon(self, exon: int) -> int:
        """Returns the number of the first codon for the given exon"""
        return sum(y for x,y in self.ref_tr.items() if x < exon)

    def _compensation_lines(self, id: str) -> Iterable[str]:
        """
        For a given Compensation entry, return the SVG line depicting an arc
        connecting the compensating frameshif mutations
        """
        compensated_muts: List[str] = [x for x in self.comp2mut[id] if x in self.mut2exon]
        for i, first_mut in enumerate(compensated_muts[:-1]):
            second_mut: str = compensated_muts[i+1]
            arc_start_x: float = self.mut_objects[first_mut].x
            arc_end_x: float = self.mut_objects[second_mut].x
            arc_start_exon: int = self.mut2exon[first_mut]
            arc_start_y: float = self.exon_boxes[arc_start_exon].y
            arc_end_exon: int = self.mut2exon[second_mut]
            arc_end_y: float = self.exon_boxes[arc_end_exon].y
            arc_center_y: float = None
            arc_center_x: float = (arc_start_x + arc_end_x) // 2
            if arc_start_y == arc_end_y:
                arc_center_y = arc_start_y + HALF_UTR_HEIGHT
            else:
                arc_center_y = arc_start_y + 2 * HALF_EXON_HEIGHT + INTER_FRAGMENT_SPACE // 2
            loc: str = COMPENSATION_ARC_COORDS.format(
                arc_start_x, arc_start_y, arc_center_x, arc_center_y,
                arc_center_x, arc_center_y, arc_end_x, arc_end_y
            )
            comp_line: str = COMPENSATION_TEMPLATE.format(loc, COMP_INDEL_STYLE)
            yield comp_line
