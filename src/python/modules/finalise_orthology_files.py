#!/usr/bin/env python3

"""
Assign the names based on referene gene IDs to query genes, finalising the orthology step
"""

import os
from collections import defaultdict
from typing import Dict, List, Optional, Set, TextIO, Tuple, Union

import click

from .constants import Headers
from .shared import (
    CONTEXT_SETTINGS,
    CommandLineManager,
    base_proj_name,
    get_proj2trans,
    parse_single_column,
)

MANY2MANY: str = "many2many"
NONE: str = "None"
ONE2MANY: str = "one2many"
T_GENE: str = "t_gene"
ZEE: int = 26
BACKTICK: int = 96
MAX_LEGAL_COPY_NUM: int = 300


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    "orthology_results",
    type=click.File("r", lazy=True),
    metavar="ORTHOLOGY_CLASSIFICATION",
)
@click.argument(
    "query_gene_file", type=click.File("r", lazy=True), metavar="QUERY_GENE_TABLE"
)
@click.argument(
    "output_dir",
    type=click.Path(exists=False),
    metavar="OUTPUT_DIR",
    default="orthology_output",
)
@click.option(
    "--query_gene_bed_file",
    "-qb",
    type=click.File("r", lazy=True),
    metavar="QUERY_GENE_BED_FILE",
    help="A path to query gene BED file",
)
@click.option(
    "--loss_summary",
    "-l",
    type=click.File("r", lazy=True),
    metavar="LOSS_SUMMARY_FILE",
    help=(
        "A path to loss summary file. "
        "Projection loss statuses are used to discriminate between "
        "lost and missing orthologous loci"
    ),
)
@click.option(
    "--discarded_projections",
    "-d",
    type=click.File("r", lazy=True),
    metavar="DISCARDED_PROJECTIONS_LIST",
    default=None,
    show_default=True,
    help="A single-column file containing names of discarded projections",
)
@click.option(
    "--paralogs",
    "-p",
    type=click.File("r", lazy=True),
    metavar="PROCESSED_PSEUDOGENES_LIST",
    default=None,
    show_default=True,
    help=(
        "A single-column file containing names of paralogous projections. "
        'Genes comprising of these projections get the "paralog_" prefix'
    ),
)
@click.option(
    "--processed_pseudogenes",
    "-pp",
    type=click.File("r", lazy=True),
    metavar="PROCESSED_PSEUDOGENES_LIST",
    default=None,
    show_default=True,
    help=(
        "A single-column file containing names of processed pseudogene/retrogene projections. "
        'Genes comprising of these projections get the "retro_" prefix'
    ),
)
@click.option(
    "--reference_isoforms",
    "-i",
    type=click.File("r", lazy=True),
    metavar="REF_ISOFORMS_FILE",
    default=None,
    show_default=True,
    help=(
        "A two-column tab-separate table containing gene-to-isoforms mapping for the reference. "
        "This mapping is used for proper missing/lost loci naming; if not provided, "
        "these non-functional loci will be assigned technical names with no orthology indication"
    ),
)
@click.option(
    "--log_file",
    "-lf",
    type=click.Path(exists=False),
    metavar="LOG_FILE",
    default=None,
    show_default=True,
    help="A path to write execution log to",
)
@click.option(
    "--log_name",
    "-ln",
    type=str,
    metavar="STR",
    default=None,
    show_default=True,
    help="Logger name to use; relevant only upon main class import",
)
@click.option(
    "--verbose",
    "-v",
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help="Controls execution verbosity",
)
class QueryGeneNamer(CommandLineManager):
    """
    Renames the query genes from the reg_{x} notation after the reference genes they correspond to:\n
        * Genes with a sole orthologous gene in the query gain its name after it;\n
        * Genes with up to three orthologs in the reference get a composite name;\n
        * Genes with more than three orthologs are named after
        the highest scoring orthology with plus symbol at the end\n
        Arguments are:\n
        * ORTHOLOGY_CLASSIFICATION is the final orthology file created by TOGA2
        (orthology_classification.tsv by default);\n
        * QUERY_GENE_TABLE is a query gene-to-transcript (projection) mapping file created by TOGA2
        (query_genes.tsv by default);\n
        * OUTPUT_DIR is an output directory to store the results in
    """

    __slots__: Tuple[str] = (
        "log_file",
        "orthology_file",
        "query_gene_table",
        "query_gene_bed",
        "gene2new_name",
        "tr2new_gene_name",
        "ref_gene2tr",
        "ref_tr2gene",
        "discarded",
        "paralogs",
        "ppgenes",
        "loss_status",
    )

    def __init__(
        self,
        orthology_results: click.File,
        query_gene_file: click.File,
        output_dir: click.Path,
        query_gene_bed_file: Optional[click.File],
        loss_summary: Optional[click.File],
        discarded_projections: Optional[click.File],
        paralogs: Optional[click.File],
        processed_pseudogenes: Optional[click.File],
        reference_isoforms: Optional[click.File],
        log_file: Optional[click.Path],
        log_name: Optional[str],
        verbose: Optional[bool],
    ) -> None:
        self.v: bool = verbose
        self.log_file: Union[str, None] = log_file
        self.set_logging(name=log_name, toga_module="query_gene_naming")

        self._to_log("Creating output directory for orthology wrap-up")
        self._mkdir(output_dir)
        self.orthology_file: str = os.path.join(
            output_dir, "orthology_classification.tsv"
        )
        self.query_gene_table: str = os.path.join(output_dir, "query_genes.tsv")
        self.query_gene_bed: str = os.path.join(output_dir, "query_genes.bed")

        self.ref_gene2tr: Dict[str, List[str]] = defaultdict(list)
        self.ref_tr2gene: Dict[str, str] = {}
        if reference_isoforms is not None:
            self._to_log("Parsing the reference isoforms file")
            self.parse_ref_isoforms(reference_isoforms)

        self.discarded: Set[str] = set()
        if discarded_projections is not None:
            self._to_log("Parsing the discarded projections file")
            self.discarded = parse_single_column(discarded_projections)

        self.ppgenes: Set[str] = set()
        if processed_pseudogenes is not None:
            self._to_log("Parsing the processed pseudogene list file")
            self.ppgenes = parse_single_column(processed_pseudogenes)
        self.paralogs: Set[str] = set()
        if paralogs is not None:
            self._to_log("Parsing the paralog list file")
            self.paralogs = parse_single_column(paralogs)

        self.loss_status: Dict[str, str] = {}
        if loss_summary is not None:
            self._to_log("Extracting loss status for query gene naming")
            self.parse_loss_summary(loss_summary)

        self._to_log(
            "Inferring orthology-based query gene names "
            "and adding them to the orthology classification file"
        )
        self.gene2new_name: Dict[str, str] = {}
        self.tr2new_gene_name: Dict[str, str] = {}
        self.modify_orthology_classification(orthology_results)
        self._to_log("Renaming query genes in the gene-to-transcript mapping file")
        self.modify_gene_mapping(query_gene_file)
        if query_gene_bed_file:
            self._to_log("Renaming query genes in the gene BED file")
            self.modify_gene_bed(query_gene_bed_file)

    def parse_ref_isoforms(self, file: TextIO) -> None:
        """Parses the reference isoforms file"""
        for i, line in enumerate(file, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 2:
                self._die(
                    (
                        "Improper formatting at isoforms file line %i; "
                        "expected 2 fields, got %i"
                    )
                    % (i, len(data))
                )
            gene, tr = data
            self.ref_tr2gene[tr] = gene

    def parse_loss_summary(self, file: TextIO) -> None:
        """Extracts loss statuses for query projections"""
        for i, line in enumerate(file, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) != 3:
                self._die(
                    (
                        "Improper formatting at loss summary file line %i; "
                        "expected 3 fields, got %i"
                    )
                    % (i, len(data))
                )
            if data[0] != "PROJECTION":
                continue
            name: str = data[1]
            status: str = data[2]
            self.loss_status[name] = status

    def modify_orthology_classification(self, file: TextIO) -> None:
        """
        Parses the original orthology file, inferres query gene names,
        and writes the updated orthology file
        """
        one2many_naming: Dict[str, int] = defaultdict(int)
        many2many_naming: Dict[str, int] = defaultdict(int)
        query_gene2orth_data: Dict[str, List[str]] = defaultdict(list)
        for i, line in enumerate(file, start=1):
            data: List[str] = line.rstrip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 5:
                self._die(
                    "Line %i in the orthology classification file has less than 5 fields"
                    % i
                )
            if data[0] == T_GENE:
                continue
            query_gene: str = data[2]
            query_gene2orth_data[query_gene].append(data)
        with open(self.orthology_file, "w") as h:
            h.write(Headers.ORTHOLOGY_TABLE_HEADER)
            for gene, lines in query_gene2orth_data.items():
                if gene == NONE:
                    new_query_name: str = NONE
                else:
                    ref_gene_names: List[str] = []
                    lines = sorted(
                        lines, key=lambda x: int(x[3].split("#")[-1].split(",")[0])
                    )
                    for line in lines:
                        if line[0] not in ref_gene_names:
                            ref_gene_names.append(line[0])
                    status: str = lines[0][4]
                    upd_ref_gene_names: List[str] = []
                    for ref_gene in ref_gene_names:
                        if status == ONE2MANY:
                            one2many_naming[ref_gene] += 1
                            suff: str = self._one2many_suffix(one2many_naming[ref_gene], ref_gene)
                            # ref_gene += f"_{chr(96 + one2many_naming[ref_gene])}"
                            ref_gene += f"_{suff}"
                        elif status == MANY2MANY:
                            many2many_naming[ref_gene] += 1
                            ref_gene += f"_{many2many_naming[ref_gene]}"
                        upd_ref_gene_names.append(ref_gene)
                    if len(upd_ref_gene_names) == 1:
                        new_query_name: str = upd_ref_gene_names.pop()
                    elif len(upd_ref_gene_names) <= 3:
                        new_query_name: str = ",".join(upd_ref_gene_names)
                    else:
                        new_query_name: str = upd_ref_gene_names[0] + "+"
                    self.gene2new_name[gene] = new_query_name
                for line in lines:
                    line[2] = new_query_name
                    h.write("\t".join(line) + "\n")
                    self.tr2new_gene_name[line[3]] = new_query_name

    def modify_gene_mapping(self, file: TextIO) -> None:
        """
        Modifies the query gene-to-transcript mapping file, substituting reg_{x}
        symbols with the corresponding orthology-based names
        """
        gene2tr: Dict[str, List[str]] = defaultdict(list)
        ref_gene_counter: Dict[str, int] = defaultdict(int)
        with open(self.query_gene_table, "w") as h:
            h.write(Headers.QUERY_GENE_HEADER)
            for i, line in enumerate(file, start=1):
                data: List[str] = line.rstrip().split("\t")
                if not data or not data[0]:
                    continue
                if len(data) < 2:
                    self._die(
                        "Line %i in the query gene mapping file has less than two columns"
                        % i
                    )
                proj: str = data[1]
                if proj in self.discarded or base_proj_name(proj) in self.discarded:
                    self._to_log(
                        "Removing projection %s from the final isoform file" % proj,
                        "warning",
                    )
                    continue
                old_name: str = data[0]
                ## gene has no orthologs in the orthology classification file
                if old_name not in self.gene2new_name:
                    ## another workaround for patching purposes
                    if proj in self.tr2new_gene_name:
                        gene_name: str = self.tr2new_gene_name[proj]
                    else:
                        gene2tr[old_name].append(data[1])
                        continue
                else:
                    gene_name: str = self.gene2new_name[old_name]
                data[0] = gene_name
                h.write("\t".join(data) + "\n")
            ## rename the remaining loci
            for gene, projs in gene2tr.items():
                # projs = [base_proj_name(x) for x in projs]
                ## check if any of the projections must be removed from the final file
                projs = [
                    x for x in projs if base_proj_name(x) not in self.discarded or x in self.discarded
                ]
                if not projs:
                    self._to_log(
                        "Dropping gene %s which no longer has any valid projections"
                        % gene
                    )
                    continue
                trs: List[str] = [get_proj2trans(x)[0] for x in projs]
                genes: Set[str] = {
                    self.ref_tr2gene.get(base_proj_name(x), base_proj_name(x))
                    for x in trs
                }
                upd_genes: List[str] = []
                for _gene in genes:
                    ref_gene_counter[_gene] += 1
                    if ref_gene_counter[_gene] > 1:
                        _gene = f"{_gene}_{ref_gene_counter[_gene]}"
                    upd_genes.append(_gene)
                if len(upd_genes) == 1:
                    new_query_name = upd_genes.pop()
                elif len(upd_genes) <= 3:
                    new_query_name = ",".join(upd_genes)
                else:
                    new_query_name = upd_genes.pop() + "+"
                ## check if this is a processed pseudogene/retrogene locus
                ppnum: int = sum(base_proj_name(x) in self.ppgenes for x in projs)
                if ppnum == len(projs):
                    new_query_name = f"retro_{new_query_name}"
                else:
                    projs = [x for x in projs if base_proj_name(x) not in self.ppgenes]
                    if not projs:
                        self._to_log(
                            "Dropping gene %s which no longer has any valid projections"
                            % gene
                        )
                        continue
                    paralog_num: int = sum(
                        base_proj_name(x) in self.paralogs for x in projs
                    )
                    if paralog_num == len(projs):
                        new_query_name = f"paralog_{new_query_name}"
                    else:
                        projs = [
                            x for x in projs if base_proj_name(x) not in self.paralogs
                        ]
                        if not projs:
                            continue
                        ## a workaround for post-hoc patching purposes (patch for v2.0.6)
                        has_known_name: bool = False
                        known_names: Set[str] = {
                            self.tr2new_gene_name.get(base_proj_name(x), "")
                            for x in projs
                        }
                        if len(known_names) == 1:
                            candidate_name: str = known_names.pop()
                            if candidate_name and "reg_" not in candidate_name:
                                new_query_name = candidate_name
                                has_known_name = True
                        if not has_known_name:
                            statuses: List[str] = [
                                self.loss_status.get(base_proj_name(x), "L")
                                for x in projs
                            ]
                            if all(x == "M" for x in statuses):
                                new_query_name = f"missing_{new_query_name}"
                            else:
                                new_query_name = f"lost_{new_query_name}"
                self.gene2new_name[gene] = new_query_name
                for proj in projs:
                    h.write(new_query_name + "\t" + proj + "\n")

    def modify_gene_bed(self, file: TextIO) -> None:
        """
        Modifies the query gene BED file, substituting reg_{x}
        symbols with the corresponding orthology-based names
        """
        ## safeguard against duplicates
        encountered_coords: Set[str] = set()
        with open(self.query_gene_bed, "w") as h:
            for i, line in enumerate(file, start=1):
                data: List[str] = line.rstrip().split("\t")
                if not data or not data[0]:
                    continue
                if len(data) < 4:
                    self._die(
                        "Line %i in the query gene BED file has less than four columns"
                        % i
                    )
                ## generally the absence of the reg_ name might result only from the gene 
                ## having all of its projections discarded for one reason or another
                old_name: str = data[3]
                key: Tuple[str, str, str] = tuple(data[:3])
                if key in encountered_coords:
                    self._to_log(
                        'Duplicate gene coordinates encountered: %s' % (f'{key[0]}:{key[1]}-{key[2]}')
                    )
                encountered_coords.add(key)
                if old_name not in self.gene2new_name:
                    self._to_log(
                        "Gene %s at line %i in the gene table has no proven orthologs in the reference"
                        % (old_name, i),
                        "warning",
                    )
                    # gene_name: str = old_name
                    ## at some point, these genes were still recorded in the output
                    ## but now we assume them to be dead on arrival
                    continue
                else:
                    gene_name: str = self.gene2new_name[old_name]
                data[3] = gene_name
                h.write("\t".join(data) + "\n")


    def _one2many_suffix(self, num: int, name: str) -> str:
        """
        Returns the one-letter suffix for one2many genes

        Args:
            num: Current gene number in the co-ortholog family
            name: Reference gene name

        Returns:
            Letter a-z if number is equal to or lower than 26, a binary combination of {a-z}{a-z} otherwise

        Fails:
            If a number of copies exceeds the maximal allowed copy number
        """
        if num <= ZEE:
            return chr(BACKTICK + num)
        if num >= MAX_LEGAL_COPY_NUM:
            self._die(
                "Maximum one2many copy number (%i) exceeded for gene %s (%i genes found)" % 
                (MAX_LEGAL_COPY_NUM, name, num)
            )
        first_letter: int = chr(BACKTICK + num // ZEE)
        last_letter: int = chr(BACKTICK + num % ZEE + int(num // ZEE > 1))
        return first_letter + last_letter


if __name__ == "__main__":
    QueryGeneNamer()
