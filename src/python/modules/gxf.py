"""
A collection of executables for annotation retrieval from GTF/GFF3 (GXF) files
"""

from dataclasses import dataclass, field


class Gxf:
    GENE: str = "gene"
    TRANSCRIPT: str = "transcript"
    EXON: str = "exon"
    CDS: str = "CDS"
    STOP_CODON: str = "stop_codon"
    CHILD_TYPES: tuple[str, str] = (CDS, EXON, STOP_CODON)
    ID: str = "ID"
    PARENT: str = "Parent"
    GTF: str = "gtf"
    GFF3: str = "gff"
    GFX: tuple[str, str] = (GTF, GFF3)
    BED: str = "bed"
    ANNOT_FORMAT: str = (BED, GTF, GFF3)
    BIOTYPES: tuple[str, ...] = (
        "biotype", "gene_biotype", "transcript_biotype", "transcript_type", "gbkey"
    )
    GENE_ID: str = "gene_id"
    TR_ID: str = "transcript_id"
    TR_PREFICES: tuple[str, ...] = ("transcript:", "rna-", "gene-",)
    EXON_PARENTS: tuple[str, ...] = ("transcript_id", PARENT, ID,)
    TR_PARENTS: tuple[str, ...] = ("gene_id", "gene_name", "gene", "Name", ID,)
    CODING_TAGS: tuple[str, str] = (
        "mRNA", 
        "protein_coding", 
        "transcript",
        "nonsense_mediated_decay",
        "protein_coding_LoF",
    )
    VDJ_TAGS: tuple[str, ...] = (
        "V_segment", 
        "C_region", 
        "V_gene_segment", 
        "C_gene_segment", 
        "IG_V_gene", 
        "TR_C_gene",
        "IG_C_gene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_LV_gene",
        "IG_V_gene",
        "TR_C_gene",
        "TR_J_gene",
        "TR_V_gene",
        "TR_D_gene",
    )

    @classmethod
    def remove_tr_prefices(cls, name: str) -> str:
        for prefix in cls.TR_PREFICES:
            name = name.removeprefix(prefix)
        return name

    @classmethod
    def get_exon_parent(cls, attrs: dict[str, str]) -> str:
        for parent_slot in cls.EXON_PARENTS:
            if attrs.get(parent_slot):
                return cls.remove_tr_prefices(attrs[parent_slot])
        return ""

    @classmethod
    def get_biotype(cls, attrs: dict[str, str]) -> str | None:
        for bio_slot in cls.BIOTYPES:
            if attrs.get(bio_slot):
                return attrs[bio_slot]
        return

    @classmethod
    def get_transcript_parent(cls, attrs: dict[str, str]):
        """
        Prioritize gene name from attributes

        Args:
            attrs: dictionary representation of a Gxf item's attribute field

        Returns:
            str: gene name, "Unknown" if all attributes are missing or empty 
        """
        for key in cls.TR_PARENTS:
            if attrs.get(key):
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
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds_start: int | None = None
    cds_end: int | None = None
    has_cds: bool = False

    def add_exon(self, start: int, end: int) -> None:
        """
        Adds a new exon to the exon collection, merging the overlapping intervals

        Arguments:
            start: int, new exon's start coordinate
            end: int, new exon's end coordinate

        Returns:
            None

        Note:
            Does not account for nested exons
        """ 
        merged: list[tuple[int, int]] = []
        for s, e in self.exons:
            if s <= end and start <= e:
                start, end = min(start, s), max(end, e)
            else:
                merged.append((s, e))
        merged.append((start, end))
        merged.sort()
        self.exons = merged

    def add_cds(self, start: int, end: int, stop: bool = False) -> None:
        """
        Adds coding sequence exon/frame/codon to the transcript object to ``self.exons`` 
        and corrects the coding sequence boundaries in the transcript

        Arguments:
            start: int, CDS item start
            end: int, CDS item end
            stop: bool, indicator whether the added item is a stop codon. 
                Individual stop codon does not count towards the presence of the coding sequence 

        Returns:
            None
        """
        self.cds_start = start if self.cds_start is None else min(self.cds_start, start)
        self.cds_end = end if self.cds_end is None else max(self.cds_end, end)
        self.has_cds |= not stop
        self.add_exon(start, end)

    def add_stop(self, start: int, end: int) -> None:
        """
        A shortcut for stop codon addition with ``add_cds()``

        Arguments:
            start: int, stop codon start
            end: int, stop codon end

        Returns:
            None
        """
        self.add_cds(start, end, stop=True)



def parse_gxf_attrs(attrs: str, gff3: bool = False) -> dict[str, str]:
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
    pairs: list[list[str]] = [
        x.strip().split(sep, 1)
        for x in attrs.strip().split(";")
        if x.strip()
    ]
    return {x[0].strip(): x[1].strip().strip('"') for x in pairs if len(x) == 2}