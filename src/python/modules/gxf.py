"""
A collection of executables for annotation retrieval from GTF/GFF3 (GXF) files
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Union


class Gxf:
    GENE: str = "gene"
    TRANSCRIPT: str = "transcript"
    EXON: str = "exon"
    CDS: str = "CDS"
    EXON_TYPES: Tuple[str, str] = (EXON, CDS)
    ID: str = "ID"
    GTF: str = "gtf"
    GFF3: str = "gff"
    GFX: Tuple[str, str] = (GTF, GFF3)
    BED: str = "bed"
    ANNOT_FORMAT: str = (BED, GTF, GFF3)
    BIOTYPES: Tuple[str, ...] = ("biotype", "gene_biotype", "transcript_biotype",)
    GENE_ID: str = "gene_id"
    TR_ID: str = "transcript_id"
    TR_PREFICES: Tuple[str, ...] = ("transcript:", "rna-", "gene-",)
    EXON_PARENTS: Tuple[str, ...] = ("transcript_id", "Parent", "ID",)
    TR_PARENTS: Tuple[str, ...] = ("gene_name", "gene", "gene_id", "Name", "ID",)
    CODING_TAGS: Tuple[str, str] = ("mRNA", "protein_coding", "transcript",)
    VDJ_TAGS: Tuple[str, ...] = ("V_segment", "C_region", "V_gene_segment", "C_gene_segment",)

    @classmethod
    def get_exon_parent(cls, attrs: Dict[str, str]) -> str:
        for parent_slot in cls.EXON_PARENTS:
            if parent_slot in attrs and attrs[parent_slot]:
                return attrs[parent_slot]
        return ""

    @classmethod
    def get_biotype(cls, attrs: Dict[str, str]) -> Union[str, None]:
        for bio_slot in cls.BIOTYPES:
            if bio_slot in attrs and attrs[bio_slot]:
                return attrs[bio_slot]
        return

    @classmethod
    def remove_tr_prefices(cls, name: str) -> str:
        for prefix in cls.TR_PREFICES:
            if name.startswith(prefix):
                name = name[len(prefix):]
        return name

    @classmethod
    def get_transcript_parent(cls, attrs: Dict[str, str]):
        """
        Prioritize gene name from attributes

        Args:
            attrs: dictionary representation of a Gxf item's attribute field

        Returns:
            str: gene name, "Unknown" if all attributes are missing or empty 
        """
        for key in cls.TR_PARENTS:
            if key in attrs and attrs[key]:
                return attrs[key]
        return "Unknown"


@dataclass(slots=True)
class GxfTranscript:
    chrom: str
    start: int
    end: int
    strand: bool
    name: str
    parent: str
    is_coding: bool
    exons: List[Tuple[int, int]] = field(default_factory=list)
    cds_start: Union[int, None] = None
    cds_end: Union[int, None] = None


def parse_gxf_attrs(attrs: str, gff3: bool = False) -> Dict[str, str]:
    """
    Parses the GXF-formatted attribute field.

    Args:
        *attrs*: a GXF-formatted attribute field (column 9 in the GXF file).
        Individual attributes are expected to be separated by a semicolon (";") 
        or a semicolon with a whitespace ("; "). 
        Attribute names and values are expected to be separated with a format-specific symbol 
        (see below).

        *gff3*: a boolean argument indicating whether the input attribute field comes from 
        a GFF3 file; otherwise the string is treated as GTF-formatted. For GTF files, attribute 
        names and values are expected to be separated by whitespace (" "), while GFF3 files are 
        expected to use equal sign ("=") as key-value separator. Set to False by default.

    Returns:
        dictionary {attribute_name: attribute_value}
    """
    sep: str = "=" if gff3 else " "
    pairs: List[List[str]] = [
        x.strip().split(sep, 1)
        for x in attrs.strip().split(";")
        if x.strip()
    ]
    return {x[0].strip(): x[1].strip().strip('"') for x in pairs if len(x) == 2}