"""
Constants for CESAR wrapper module
"""

from typing import Dict, Set, Tuple

import os

## default values for profile files
DEF_BLOSUM_FILE: Tuple[str] = ('supply', 'BLOSUM62.txt')
HL_COMMON_ACCEPTOR: Tuple[str] = ('extra', 'tables', 'human', 'acc_profile.txt')
HL_COMMON_DONOR: Tuple[str] = ('extra', 'tables', 'human', 'do_profile.txt')
HL_FIRST_ACCEPTOR: Tuple[str] = ('extra', 'tables', 'human', 'firstCodon_profile.txt')
HL_LAST_DONOR: Tuple[str] = ('extra', 'tables', 'human', 'lastCodon_profile.txt')
HL_EQ_ACCEPTOR: Tuple[str] = ('supply', 'eq_acc_profile.txt')
HL_EQ_DONOR: Tuple[str] = ('supply', 'eq_donor_profile.txt')

CESAR_PROFILE_DIR: Tuple[str] = ('supply', 'CESAR2.0', 'profiles')
HG38_CANON_U2_ACCEPTOR: Tuple[str] = (
    *CESAR_PROFILE_DIR, 'human', 'canon_U2_acceptor.tsv'
)
HG38_CANON_U2_DONOR: Tuple[str] = (
    *CESAR_PROFILE_DIR, 'human', 'canon_U2_donor.tsv'
)
HG38_NON_CANON_U2_ACCEPTOR: Tuple[str] = (
    *CESAR_PROFILE_DIR, 'human', 'nonCanon_U2_acceptor.tsv'
)
HG38_NON_CANON_U2_DONOR: Tuple[str] = (
    *CESAR_PROFILE_DIR, 'human', 'nonCanon_U2_donor.tsv'
)
HG38_CANON_U12_ACCEPTOR: Tuple[str] = (
    *CESAR_PROFILE_DIR, 'human', 'canon_U12_acceptor.tsv'
)
HG38_CANON_U12_DONOR: Tuple[str] = (
    *CESAR_PROFILE_DIR, 'human', 'canon_U12_donor.tsv'
)
HG38_NON_CANON_U12_ACCEPTOR: Tuple[str] = (
    *CESAR_PROFILE_DIR, 'human', 'canon_U12_acceptor.tsv'
)
HG38_NON_CANON_U12_DONOR: Tuple[str] = (
    *CESAR_PROFILE_DIR, 'human', 'canon_U12_acceptor.tsv'
)
FIRST_ACCEPTOR: Tuple[str] = (*CESAR_PROFILE_DIR, 'firstCodon_profile.tsv')
LAST_DONOR: Tuple[str] = (*CESAR_PROFILE_DIR, 'lastCodon_profile.tsv')

## warning messages for unaligned exon portions
MEM_UNALIGNED_WARNING: str = (
    'Group {} containing exons {} exceeds either memory ({}) or search space '
    '({}) limits and will not be aligned.'
)
SHORT_SPACE_UNALGNED_WARNING: str = (
    'Group {} containing exons {} has a very short query sequence ({})) '
    'and will not be aligned'
)
LARGE_EXON_UNALIGNED_WARNING: str = (
    'Group {} containing exons {} has one or several of its exon ({}) '
    'longer than 2kb missing from the chain and will not be aligned'
)

HL_CESAR_PATH: str = os.path.join(
    os.path.sep,
    'projects',
    'hillerlab',
    'genome',
    'src',
    'TOGA_pub',
    'CESAR2.0',
    'cesar'
)
HL_BW2W_PATH: str = os.path.join(
    os.path.sep,
    'projects',
    'hillerlab',
    'genome',
    'bin',
    'x86_64',
    'bigWigToWig'
)

## Shell oneliners
SPLICEAI_PROCESS_SCRIPT: str = (
    'awk \'{{if($1 ~ /coordinates/){{split($2,chr,"="); split($3,start,"="); chrom=chr[2]; '
    'pos=start[2];next}}; if($1 >= 0.001){{print chrom"\t"pos"\t"$0}}; pos++}}\''
)

## accepted U2 splice sites
ACCEPTOR_SITE: Tuple[str] = ("ag",)
DONOR_SITE: Tuple[str] = ("gt", "gc")
ACCEPTOR_SITE_U12: str = 'ag'
DONOR_SITE_U12: str = 'gt'
U2: str = 'U2'
U12: str = 'U12'

## special codon shortcuts
XXX_CODON: str = "XXX"
GAP_CODON: str = "---"
NNN_CODON: str = "NNN"
STOPS: Set[str] = {"TAG", "TGA", "TAA"}
START: str = 'ATG'

AA_CODE: Dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "---": "-",
    "NNN": "X",
}

## Preprocessing constants
MAX_CHAIN_GAP_SIZE: int = 1000000
MIN_REF_LEN_PERC: float = 0.5
SINGLE_EXON_MIN_REF_LEN_PERC: float = 0.3
MIN_PROJ_OVERLAP_THRESHOLD: float = 0.5
EXTRA_FLANK: float = 0.1
SEARCH_OFFSET: int = 12 ## offset to add to the potential gained intron sequence

## length and size constants
FLANK_SPACE: int = 50
MAX_UNCOV_LEN: int = 50
MIN_ASMBL_GAP_SIZE: int = 10
MIN_EXON_LENGTH: int = 7
MIN_INTRON_LENGTH: int = 10
ORTHOLOGY_THRESHOLD: float = 0.5
SS_SIZE: int = 2
MAX_DEV_FROM_SPLICEAI: int = 600 ## TODO: To be discussed with Michael
MAX_CHAIN_INTRON_LEN: int = 500000

## Exon identity thresholds for intact exons
MIN_ID_THRESHOLD: float = 45.0
MIN_BLOSUM_THRESHOLD: float = 20.0

## default idnel penalties taken from the parse_cesar_out.py
INS_PEN: int = -1
DEL_PEN: int = -1

## quality class-wise exon identity threholds
HQ_PID: int = 75
HQ_BLOSUM: int = 65

AB_INCL_PID: int = 25
AB_INCL_BLOSUM: int = 25

C_INCL_PID: int = 65
C_INCL_BLOSUM: int = 50

A_T_PID: int = 65
A_T_BLOSUM: int = 40

LO_T_PID: int = 45
LO_T_BLOSUM: int = 25

## mutation class shortcuts
LEFT_SPLICE_CORR: Tuple[str] = ('ag',)  # acceptor
LEFT_SPLICE_CORR_U12: str = 'ag'  # acceptor
RIGHT_SPLICE_CORR: Tuple[str] = ('gt', 'gc')
RIGHT_SPLICE_CORR_U12: str = 'gt'
MISS_EXON: str = "Missing exon"
DEL_EXON: str = "Deleted exon"
DEL_MISS: Tuple[str] = (MISS_EXON, DEL_EXON)
COMPENSATION: str = "COMPENSATION"
SSM_A: str = "SSMA"
SSM_D: str = "SSMD"
SSM: Tuple[str, str] = (SSM_A, SSM_D)
START_MISSING: str = "START_MISSING"
ATG: str = "ATG"
FS_DEL: str = "FS_DEL"
FS_INS: str = "FS_INS"
FS_INDELS: Tuple[str, str] = (FS_DEL, FS_INS)
BIG_DEL: str = "BIG_DEL"
BIG_INS: str = "BIG_INS"
BIG_INDEL: Tuple[str] = (BIG_DEL, BIG_INS)
STOP: str = 'STOP'
STOP_MISSING: str = 'STOP_MISSING'
INTRON_GAIN: str = 'INTRON_GAIN'
DEFAULT_STOP_MISSING: str = 'Missing stop masked'
ALT_MASKING_REASON: str = 'Alternative exon splitting'
INTRON_GAIN_MASK_REASON: str = 'Intron gain masked'
INTRON_DEL_REASON: str = 'Intron deletion'
COMPENSATION_REASON: str = 'Compensated'
ALT_FRAME_REASON: str = 'Alternative frame found'
EX_DEL_REASON: str = 'Exon is deleted'
EX_MISS_REASON: str = 'Exon is missing'
U12_REASON: str = 'U12 intron'
NON_CANON_U2_REASON: str = 'Non-canonical U2 intron'
OBSOLETE_COMPENSATION: str = 'Treated as inactivating'
SAFE_SPLICE_SITE_REASONS: Tuple[str, str] = (
    U12_REASON, NON_CANON_U2_REASON, INTRON_DEL_REASON
)
## reasons under which masked mutations are not unmasked 
SAFE_UNMASKABLE_REASONS: Tuple[str] = (
    U12_REASON, NON_CANON_U2_REASON, INTRON_DEL_REASON, ALT_FRAME_REASON,
    EX_DEL_REASON, EX_MISS_REASON
)
## mutation types which do not get unmasked
SAFE_UNMASKABLE_TYPES: Tuple[str] = (
    BIG_INS, BIG_DEL, COMPENSATION, MISS_EXON, START_MISSING, STOP_MISSING, INTRON_GAIN
)


## projection classification constants
## minimal present codon faction to automatically classify projection as intact
STRICT_FACTION_INTACT_THRESHOLD: float = 0.6
## minimal intact codon threshold among non-missing codons for a projection
## not to be automatically classified as lost
INTACT_CODON_LOSS_THRESHOLD: float = 0.35
## minimal absent/unaligned exon portion to classify a missing exon as partially missing
OUT_OF_CHAIN_MISSING_THRESHOLD: float = 0.5
## minimal non-deleted codon fraction for a projection
## not to be automatically classified as lost
NON_DEL_LOSS_THRESHOLD: float = 0.2
## minimal intact non-missing threshold portion for a mutation-unaffected projection
## to be not classified as an uncertain loss
MIN_INTACT_UL_FRACTION: float = 0.49
## maximal missing reference sequence portion for a projection with misssing exons
## to be classified as partially intact
MAX_MISSING_PM_THRESHOLD: float = 0.5
## maximal length of 'retained' intron to be considered a precise intron deletion
## coupled with insertion and incorrectly processed by CESAR
MAX_RETAINED_INTRON_LEN: int = 27
## maximum indel size, in base
BIG_INDEL_SIZE: int = 50
## Maximum size of an exon safe for deletion no matter what (deprecated?)
SAFE_EXON_DEL_SIZE: int = 40  # actually 39
## Maximum size of the first/last exon to be classified as Missing instead of Deleted
TERMINAL_EXON_DEL_SIZE: int = 20

## ORTHOLOGY LABELS
ORTHOLOG: str = 'ORTHOLOG'
PARALOG: str = 'PARALOG'
PROC_PSEUDOGENE: str = 'PROCESSED_PSEUDOGENE'

## projection loss classes
FI: str = 'FI'
I: str = 'I'
PI: str = 'PI'
M: str = 'M'
PM: str = 'PM'
L: str = 'L'
UL: str = 'UL'
PG: str = 'PG'
PP: str = 'PP'
N: str = 'N'
LOSS_STATUSES: Tuple[str] = (
    FI, I, PI, M, PM, L, UL, PG, PP, N
)

## numeric orthology resolution constants
MIN_COV_FOR_ORTH: float = 0.5
MAX_QLEN_FOR_ORTH: float = 0.2
MIN_INTRON_COV_FOR_ORTH: float = 0.4

## track coloring constants
DARK_BLUE: str = '0,0,100'
BLUE: str = "0,0,200"
LIGHT_BLUE: str = "0,200,255"
LIGHT_RED: str = "255,50,50"
SALMON: str = "255,160,120"
GREY: str = "130,130,130"
BROWN: str = "159,129,112"
BLACK: str = "10,10,10"
PINK: str = '250,50,200'
CLASS_TO_COL: Dict[str, str] = {
    PP: PINK,
    PG: BROWN,
    PM: GREY,
    L: LIGHT_RED,
    M: GREY,
    UL: SALMON,
    PI: LIGHT_BLUE,
    I: BLUE,
    FI: DARK_BLUE
}
CLASS_TO_NAME: Dict[str, str] = {
    'FI': 'Fully Intact',
    'I': 'Intact',
    'PI': 'Partially Intact',
    'UL': 'Uncertain Loss',
    'M': 'Missing',
    'L': 'Lost',
    'PG': 'Paralogous Projection',
    'PP': 'Processed Pseudogene'
}
NUM_TO_CLASS: Dict[int, str] = {
    -1: N, 0: PP, 1: PG, 2: PM, 3: M, 4: L, 5: UL, 6: PI, 7: I, 8: FI
}
CLASS_TO_NUM: Dict[str, int] = {v: k for k, v in NUM_TO_CLASS.items()}