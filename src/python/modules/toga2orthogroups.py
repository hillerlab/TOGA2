#!/usr/bin/env python3
"""
toga2orthogroups — Build orthogroups from TOGA2 pairwise orthology annotations.

Constructs orthogroups by modelling the (reference_gene - query_gene) orthology
relationships as a graph and extracting connected components with Union-Find.
Optionally merges with PANTHER family assignments.

Usable as:
    - CLI:  python toga2orthogroups.py -t DIR -s FILE -b FILE -i FILE -o DIR
    - Module:  from src.python.modules.toga2orthogroups import build_orthogroups, generate_count_table
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass, field
from math import inf, sqrt


__author__ = "Michele Albertini"
__email__ = "michelealbertini30@gmail.com"
__github__ = "https://github.com/michelealbertini30"


# Module-level logger
log = logging.getLogger(__name__)

# Some TOGA2 orthology_classification.tsv files contain very long fields.
# Raise the csv field size limit to the largest value the platform allows.
_csv_limit = sys.maxsize
while True:
    try:
        csv.field_size_limit(_csv_limit)
        break
    except OverflowError:
        _csv_limit //= 2

# QC: FamZ
_FAMILY_COUNT_Z = 3.0       # |z| threshold for flagging a species in a single family
_MAX_ZERO_FRAC = 0.5        # skip families where >= this fraction of species have 0 copies

# QC: OrthoZ
_MIN_SPAN_SPECIES = lambda n: max(3, min(n // 3, 10))   # min species-with-genes for a scoreable family

# ---------------------------------------------------------------------------
# Union-Find
# ---------------------------------------------------------------------------
# Union-Find tracks which elements belong to the same connected component.
#
# Core architecture:
#   - Every reference gene starts as its own singleton component
#   - For each query species, an inverted index is built: query gene -> ref genes
#   - If a query gene shares multiple ref genes, those get unioned in the same component
#   - This is repeated for all species, modifying the same components
#   - After all species are processed, connected components are extracted (orthogroups)

class UnionFind:
    """Weighted Union-Find with path compression"""

    __slots__ = ("_parent", "_rank")

    def __init__(self) -> None:
        # _parent[x] -> points towards the root of component x.
        self._parent: dict[str, str] = {}

        # _rank[x] -> height of the subtree rooted at x.
        self._rank: dict[str, int] = {}

    def find(self, x: str) -> str:
        """Find representative with path compression.

        If x has not been seen before, register it as a new singleton.
        Otherwise, walk up the parent chain to find the root, then
        compress the path so every node on the way points directly to the root.
        """
        # First encounter: initialise x as its own parent with rank 0.
        if x not in self._parent:
            self._parent[x] = x
            self._rank[x] = 0
            return x

        # Iterative root-finding.
        root = x
        while self._parent[root] != root:
            root = self._parent[root]

        # Path compression.
        while self._parent[x] != root:
            self._parent[x], x = root, self._parent[x]

        return root

    def union(self, a: str, b: str) -> None:
        """Merge the components containing a and b (union by rank).

        Otherwise, attach the root with the lower rank under the root with
        the higher rank. If ranks are equal, pick arbitrarily and increment
        the winner's rank.
        """
        ra, rb = self.find(a), self.find(b)

        # Already in the same component
        if ra == rb:
            return

        # Ensure ra is always the higher-rank root
        if self._rank[ra] < self._rank[rb]:
            ra, rb = rb, ra

        # Attach rb's tree under ra.
        self._parent[rb] = ra

        # Increment rank when the two trees have equal rank
        if self._rank[ra] == self._rank[rb]:
            self._rank[ra] += 1

    def components(self) -> dict[str, list[str]]:
        """Return {root: [members]} for every element currently tracked.

        Each key is the canonical root of one connected component.
        The associated list contains all members of that component,
        including the root itself.
        """
        groups: dict[str, list[str]] = defaultdict(list)
        for x in self._parent:
            # Add each member to its root's list
            # self.find(x) returns the root
            groups[self.find(x)].append(x)
        return dict(groups)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class ReferenceGeneSet:
    """Reference genes and their transcripts, filtered by a user-supplied blacklist.

    Transcripts on blacklisted chromosomes/scaffolds are excluded before
    building orthogroups to avoid copy-number biases.
    """
    genes: set[str]
    transcripts: set[str]
    gene_to_transcripts: dict[str, list[str]] = field(default_factory=dict)


@dataclass
class SpeciesOrthologs:
    """Per-species mapping of reference genes to query genes.

    Populated from two TOGA2 output files for each query species:
      - loss_summary.tsv -> which query transcripts have acceptable loss status.
      - orthology_classification.tsv → ref_gene + query_gene pairs
    """
    species: str
    # Set of query gene IDs with at least one intact transcript.
    ref_to_query: dict[str, set[str]]


@dataclass
class Orthogroups:
    """Connected-component orthogroups with per-species copy-number info.

    The families dict encodes the full membership needed to build count table.
    """
    # family_id -> set of "species|query_gene" strings, one per orthologous gene
    families: dict[str, set[str]]
    # Reference gene to family membership
    ref_gene_to_family: dict[str, str]
    # ref_gene -> set of "species|query_gene" orthologs (family-level QC)
    ref_gene_orthologs: dict[str, set[str]] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Main Functions
# ---------------------------------------------------------------------------

def _parse_species_list(value: str) -> list[str]:
    """Return an ordered list of species names.

    Accepts either a comma-separated string (e.g. 'human,mouse') or a path to a
    file with one name per line (tab-separated lines use the first column only).
    """
    p = Path(value)
    if p.is_file():
        result = []
        for line in p.read_text().splitlines():
            sp = line.strip().split("\t")[0]
            if sp:
                result.append(sp)
        return result
    return [sp.strip() for sp in value.split(",") if sp.strip()]


def _parse_blacklist(value: str | None) -> set[str]:
    """Return a set of chromosome/scaffold names to exclude.

    Accepts either a comma-separated string (e.g. 'chrX,chrY') or a path to a
    file with one name per line.
    """
    if not value:
        return set()
    p = Path(value)
    if p.is_file():
        return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}
    return {name.strip() for name in value.split(",") if name.strip()}


def load_reference_genes(
    isoforms_path: str | Path,
    transcripts_bed_path: str | Path,
    blacklist: set[str] | None = None,
) -> ReferenceGeneSet:
    """Load TOGA isoforms, filter out blacklisted chromosomes, return gene set.

    Parameters
    ----------
    isoforms_path : path to TOGA isoforms TSV (gene \\t transcript)
    transcripts_bed_path : path to transcript BED (chr in col1, name in col4)
    blacklist : chromosome/scaffold names to exclude (default: none)
    """
    _blacklist = blacklist or set()

    # --- Collect non-blacklisted transcript IDs from the BED file ---
    kept_transcripts: set[str] = set()
    n_excluded = 0

    with open(transcripts_bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split("\t", 5)
            chrom, name = parts[0], parts[3]
            if chrom in _blacklist:
                n_excluded += 1
            else:
                kept_transcripts.add(name)

    log.info(
        "Transcripts BED: %d kept, %d excluded (%s)",
        len(kept_transcripts), n_excluded, ", ".join(sorted(_blacklist)) or "none",
    )

    # --- Read isoforms TSV, keep only entries with at least one non-blacklisted transcript ---
    genes: set[str] = set()
    gene_to_tx: dict[str, list[str]] = defaultdict(list)
    n_iso_total = 0
    n_iso_kept = 0

    with open(isoforms_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            n_iso_total += 1
            gene, transcript = row[0], row[1]
            if transcript in kept_transcripts:
                n_iso_kept += 1
                genes.add(gene)
                gene_to_tx[gene].append(transcript)

    log.info(
        "TOGA isoforms: %d total, %d kept, %d genes",
        n_iso_total, n_iso_kept, len(genes),
    )

    return ReferenceGeneSet(
        genes=genes,
        transcripts=kept_transcripts,
        gene_to_transcripts=dict(gene_to_tx),
    )


def load_species_orthologs(
    species: str,
    toga_dir: str | Path,
    reference_genes: set[str],
    include_ul: bool = False,
) -> SpeciesOrthologs:
    """Load orthology data filtered by loss status.

    Only keeps relationships where the query transcript loss status is
    FI/I/PI, and additionally UL if --include-ul is provided.
    """
    toga_dir = Path(toga_dir)
    species_dir = toga_dir / species
    # The set of loss status codes we consider as gene presence (intact).
    intact_statuses = {"I", "FI", "PI"}
    if include_ul:
        intact_statuses.add("UL")

    # --- Parse loss_summary.tsv to find intact query transcripts ---
    intact_transcripts: set[str] = set()
    loss_path = species_dir / "loss_summary.tsv"

    with open(loss_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 3 and row[2] in intact_statuses:
                intact_transcripts.add(row[1])

    log.debug(
        "  %s loss_summary: %d intact transcripts",
        species, len(intact_transcripts),
    )

    # --- Parse orthology_classification.tsv ---
    ref_to_query: dict[str, set[str]] = defaultdict(set)
    ortho_path = species_dir / "orthology_classification.tsv"
    n_rows = 0
    n_kept = 0   # not None, intact, autosomal

    with open(ortho_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 5:
                continue
            n_rows += 1
            t_gene, q_gene, q_transcript = row[0], row[2], row[3]

            if q_gene == "None":
                continue

            if q_transcript not in intact_transcripts:
                continue

            if t_gene not in reference_genes:
                continue

            n_kept += 1
            ref_to_query[t_gene].add(q_gene)

    log.debug(
        "  %s orthology: %d rows, %d kept, %d ref genes with orthologs",
        species, n_rows, n_kept, len(ref_to_query),
    )

    return SpeciesOrthologs(species=species, ref_to_query=dict(ref_to_query))


def build_orthogroups(
    ref_genes: ReferenceGeneSet,
    species_list: list[str],
    toga_dir: str | Path,
    panther_path: str | Path | None = None,
    isoforms_path: str | Path | None = None,
    include_ul: bool = False,
) -> Orthogroups:
    """Build gene families using Union-Find over reference genes.

    For each query species, reference genes that share any query ortholog are
    unioned into the same component. Cross-species families emerge naturally
    as a single Union-Find structure.Optionally, PANTHER assignments add
    extra ref→ref edges before components are extracted.

    Returns an Orthogroups object with family memberships and a mapping from
    reference genes to family IDs.
    """

    uf = UnionFind()

    # Pre-register reference gene so that genes with no orthologs in any
    # species still appear as singleton families in the final output.
    for g in ref_genes.genes:
        uf.find(g)

    # For each reference gene, the set of "species|query_gene"
    ref_gene_orthologs: dict[str, set[str]] = defaultdict(set)

    # --- Process each species ---
    for spe in species_list:
        log.debug("Processing species: %s", spe)

        # Load the filtered ref->query mappings for this species.
        sp_ortho = load_species_orthologs(spe, toga_dir, ref_genes.genes, include_ul=include_ul)

        # Build the inverted index: query_gene -> ref_genes
        query_to_refs: dict[str, list[str]] = defaultdict(list)

        for ref_gene, q_genes in sp_ortho.ref_to_query.items():
            for qg in q_genes:
                query_to_refs[qg].append(ref_gene)
                ref_gene_orthologs[ref_gene].add(f"{spe}|{qg}")

        # Union all reference genes that share any query gene.
        for qg, ref_list in query_to_refs.items():
            if len(ref_list) > 1:
                anchor = ref_list[0]
                for rg in ref_list[1:]:
                    uf.union(anchor, rg)

        log.debug("  %s: %d query genes", spe, len(query_to_refs))

    # --- Optional PANTHER merge ---
    gene_to_panther: dict[str, str] = {}
    if panther_path is not None and isoforms_path is not None:
        log.info("Loading PANTHER families...")
        panther_families, gene_to_panther = _load_panther_families(panther_path, isoforms_path)
        for gene_list in panther_families.values():
            ref_in_family = [g for g in gene_list if g in ref_genes.genes]
            if len(ref_in_family) < 2:
                continue
            anchor = ref_in_family[0]
            for g in ref_in_family[1:]:
                uf.union(anchor, g)

    # --- Extract connected components and build family objects ---
    components = uf.components()

    families: dict[str, set[str]] = {}
    ref_gene_to_family: dict[str, str] = {}

    for _root, members in components.items():
        # Family ID: use PANTHER ID when available (PANTHER mode)
        if gene_to_panther:
            panther_ids = {gene_to_panther[m] for m in members if m in gene_to_panther}
            family_id = min(panther_ids) if panther_ids else min(members)
        else:
            family_id = min(members)

        # Aggregate all "species|query_gene" orthologs from every reference
        family_orthologs: set[str] = set()
        for rg in members:
            family_orthologs.update(ref_gene_orthologs.get(rg, set()))

        families[family_id] = family_orthologs

        # Record the family ID for every member reference gene.
        for rg in members:
            ref_gene_to_family[rg] = family_id

    log.info(
        "Orthogroups: %d families from %d reference genes",
        len(families), len(ref_gene_to_family),
    )

    return Orthogroups(
        families=families,
        ref_gene_to_family=ref_gene_to_family,
        ref_gene_orthologs=dict(ref_gene_orthologs),
    )


def _load_panther_families(
    panther_path: str | Path,
    isoforms_path: str | Path,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Load PANTHER DB and map gene symbols to reference gene IDs.

    Returns:
        family_to_genes : {panther_family_id: [ref_gene_ids]}
        gene_to_panther : {ref_gene_id: panther_family_id}

    The PANTHER flat file uses gene symbols (e.g. "BRCA2") rather than ENSG
    IDs. The isoforms TSV encodes the gene symbol in the transcript ID field
    after the last "#" character (e.g. "ENST00000380152.7#BRCA2").
    We exploit this to build a symbol -> ENSG mapping and then look up each
    reference gene in the PANTHER table.
    """
    # --- Build gene-symbol -> PANTHER family ID mapping ---
    symbol_to_family: dict[str, str] = {}

    with open(panther_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 5:
                continue
            symbol = row[2]
            # Trim to family-level ID (9 chars = "PTHRxxxxx")
            family_id = row[3][:9] if len(row[3]) >= 9 else row[3]
            symbol_to_family[symbol] = family_id

    # --- Build reference gene ID -> gene symbol mapping via isoforms ---
    gene_to_symbol: dict[str, str] = {}

    with open(isoforms_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            gene = row[0]
            symbol = row[1].rsplit("#", 1)[-1]
            if gene not in gene_to_symbol:
                gene_to_symbol[gene] = symbol

    # --- Group reference genes by their PANTHER family ---
    family_to_genes: dict[str, list[str]] = defaultdict(list)
    gene_to_panther: dict[str, str] = {}
    matched = 0  # ref genes matched to PANTHER

    for gene, symbol in gene_to_symbol.items():
        if symbol in symbol_to_family:
            pid = symbol_to_family[symbol]
            family_to_genes[pid].append(gene)
            gene_to_panther[gene] = pid
            matched += 1

    log.info(
        "PANTHER: %d symbols matched, %d families",
        matched, len(family_to_genes),
    )
    return dict(family_to_genes), gene_to_panther


def generate_count_table(
    orthogroups: Orthogroups,
    species_list: list[str],
) -> tuple[list[str], list[list]]:
    """Build the copy-number count table.

    CAFE5 expects a tab-separated file where:
      - Header: [Desc, Family ID, species1, species2, ...]
      - Rows: NA, Family ID, and copy number of that gene family in each species.
    """
    header = ["Desc", "Family ID"] + species_list
    rows = []

    # Iterate over families in sorted order so the output is deterministic
    for family_id in sorted(orthogroups.families):
        members = orthogroups.families[family_id]

        # Count unique query genes per species.
        counts: dict[str, int] = {sp: 0 for sp in species_list}
        for member in members:
            if "|" in member:
                sp = member.split("|", 1)[0]
                if sp in counts:
                    counts[sp] += 1

        rows.append(["NA", family_id] + [counts[sp] for sp in species_list])

    return header, rows


def apply_size_filter(rows: list[list], threshold: int) -> tuple[list[list], int]:
    """Drop families where any species has >= threshold gene copies.

    Row format expected: [Desc, family_id, count1, count2, ...]
    Returns (filtered_rows, n_removed).
    """
    filtered = [row for row in rows if all(c < threshold for c in row[2:])]
    return filtered, len(rows) - len(filtered)


def write_count_table(
    header: list[str],
    rows: list[list],
    output_path: str | Path,
) -> None:
    """Write copy-number count table to disk."""

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    log.info("Wrote count table: %s (%d families)", output_path, len(rows))


def write_orthogroup_membership(
    orthogroups: Orthogroups,
    output_path: str | Path,
) -> None:
    """Write orthogroup membership file.

    The first field is the family ID, remaining fields are the "species|query_gene"
    """
    with open(output_path, "w") as fh:
        for family_id in sorted(orthogroups.families):
            members = sorted(orthogroups.families[family_id])
            fh.write(family_id + "\t" + "\t".join(members) + "\n")
    log.info("Wrote membership: %s (%d families)", output_path, len(orthogroups.families))


# ---------------------------------------------------------------------------
# Species QC — inflation diagnostics
# ---------------------------------------------------------------------------

def species_qc_diagnostics(
    orthogroups: Orthogroups,
    species_list: list[str],
    ortho_z: float = 3.0,
    fam_z: float = 3.0,
    output=None,
) -> dict[str, dict]:
    """Detect species with low orthology resolution via two complementary signals.

    OrthoZ — spanning-rate signal:
        SumZ = Sum of per family sum(n_refs_spanned - 1) / total queries per species z-scores.
        OrthoZ = Z-score of SumZ across species.

        Family filtering:
        - Skipped if fewer than max(3, min(n//3, 10)) species have genes, or var==0.

    FamZ — copy-number signal:
        Per family, z-score raw copy counts across all species.
        Count families where a species is a count outlier (|z| > _FAMILY_COUNT_Z).
        Z-score those counts across species → FamZ.

        Family filtering:
        - Skipped if fewer than 50% of species have genes or var==0

    Flag a species if OrthoZ > ortho_z OR FamZ > fam_z.

    Returns
    -------
    dict mapping species -> {flagged_count, count_scoreable,
                             fam_z, ortho_z, flagged}
    """
    if output is None:
        output = sys.stderr

    n_total = len(species_list)

    min_fam_species = _MIN_SPAN_SPECIES(n_total)

    sp_set = set(species_list)

    # Invert ref_gene_to_family to iterate families by their ref genes.
    family_to_refs: dict[str, list[str]] = defaultdict(list)
    for rg, fid in orthogroups.ref_gene_to_family.items():
        family_to_refs[fid].append(rg)

    # --- Spanning-rate → OrthoZ ---
    # Capture species with high spanning depth over ref genes -> orthogroup inflation
    sum_z: dict[str, float] = defaultdict(float)
    n_scoreable_span = 0

    for _, ref_genes_fam in family_to_refs.items():
        # Skip singleton families
        if len(ref_genes_fam) < 2:
            continue

        # Build inverted index: query_gene -> {ref_genes}.
        query_to_refs: dict[str, set[str]] = defaultdict(set)
        for rg in ref_genes_fam:
            for sq in orthogroups.ref_gene_orthologs.get(rg, set()):
                query_to_refs[sq].add(rg)

        # sp_weighted accumulates (n_refs - 1) for every spanning query gene
        sp_total: dict[str, int] = defaultdict(int)
        sp_weighted: dict[str, float] = defaultdict(float)
        for sq, refs in query_to_refs.items():
            sp = sq.split("|", 1)[0]
            if sp not in sp_set:
                continue
            sp_total[sp] += 1
            if len(refs) > 1:
                sp_weighted[sp] += len(refs) - 1

        # Require at least min_fam_species species with genes for a stable z-score.
        present = [sp for sp in species_list if sp_total.get(sp, 0) > 0]
        if len(present) < min_fam_species:
            continue

        rates = {
            sp: sp_weighted.get(sp, 0.0) / sp_total[sp]
            for sp in present
        }
        vals = list(rates.values())
        mean_r = sum(vals) / len(vals)
        var_r = sum((v - mean_r) ** 2 for v in vals) / (len(vals) - 1)

        # Skip families where all species have the same spanning rate
        if var_r == 0.0:
            continue

        n_scoreable_span += 1
        std_r = sqrt(var_r)

        # Accumulate positive z-scores (inflation)
        for sp in present:
            z = (rates[sp] - mean_r) / std_r
            sum_z[sp] += max(0.0, z)

    # Cross-species z-score of SumZ values → OrthoZ.
    sum_z_vals = [sum_z.get(sp, 0.0) for sp in species_list]
    grand_mean_span = sum(sum_z_vals) / n_total if n_total > 0 else 0.0
    var_span = sum((v - grand_mean_span) ** 2 for v in sum_z_vals) / (n_total - 1) if n_total > 1 else 0.0
    span_std = sqrt(var_span) if var_span > 0 else inf

    # --- Copy-number outlier count → FamZ ---
    # Capture assembly artefacts (M, L) or inflated one-to-many
    n_count_flagged: dict[str, int] = defaultdict(int)
    n_count_scoreable = 0
    n_families = len(orthogroups.families)
    n_one2one = 0

    for members in orthogroups.families.values():
        # Count query genes per species in this family.
        sp_counts: dict[str, int] = {sp: 0 for sp in species_list}
        for m in members:
            if "|" in m:
                sp = m.split("|", 1)[0]
                if sp in sp_counts:
                    sp_counts[sp] += 1

        # One-to-one check: every species has exactly 1 copy.
        if all(sp_counts[sp] == 1 for sp in species_list):
            n_one2one += 1

        vals_c = [sp_counts[sp] for sp in species_list]

        # Skip zero-inflated families.
        n_zero = sum(1 for v in vals_c if v == 0)
        if n_zero / n_total >= _MAX_ZERO_FRAC:
            continue

        # Skip families with no copy-number variation.
        mean_c = sum(vals_c) / n_total
        var_c = sum((v - mean_c) ** 2 for v in vals_c) / (n_total - 1) if n_total > 1 else 0.0
        if var_c == 0.0:
            continue

        n_count_scoreable += 1
        std_c = sqrt(var_c)

        # Flag species deviating from _FAMILY_COUNT_Z
        for sp in species_list:
            if abs((sp_counts[sp] - mean_c) / std_c) > _FAMILY_COUNT_Z:
                n_count_flagged[sp] += 1

    # Cross-species z-score of n_count_flagged → FamZ.
    count_flagged_vals = [n_count_flagged.get(sp, 0) for sp in species_list]
    grand_mean_cf = sum(count_flagged_vals) / n_total if n_total > 0 else 0.0
    var_cf = sum((v - grand_mean_cf) ** 2 for v in count_flagged_vals) / (n_total - 1) if n_total > 1 else 0.0
    fam_std = sqrt(var_cf) if var_cf > 0 else inf

    # --- Flag species based on two metric results ---
    results: dict[str, dict] = {}
    flagged_species: list[str] = []

    for sp in species_list:
        # OrthoZ
        sz = sum_z.get(sp, 0.0)
        ortho_z_val = (sz - grand_mean_span) / span_std if span_std != inf else 0.0

        # FamZ
        cf = n_count_flagged.get(sp, 0)
        fam_z_val = (cf - grand_mean_cf) / fam_std if fam_std != inf else 0.0

        # Flag on either signal exceeding its respective threshold.
        flagged = ortho_z_val > ortho_z or fam_z_val > fam_z
        results[sp] = {
            "flagged_count": cf,
            "count_scoreable": n_count_scoreable,
            "fam_z": round(fam_z_val, 3),
            "ortho_z": round(ortho_z_val, 3),
            "flagged": flagged,
        }
        if flagged:
            flagged_species.append(sp)

    # --- Formatted report to stderr ---
    output.write("\n=== Species QC Diagnostics ===\n")
    output.write(f"  ortho-z threshold: {ortho_z}\n")
    output.write(f"  fam-z threshold:  {fam_z}\n")
    output.write(f"  one-to-one orthologs: {n_one2one}\n")
    output.write(f"  total families: {n_families}\n")
    output.write(f"  scored families (OrthoZ): {n_scoreable_span}\n")
    output.write(
        f"  {'Species':<25} {'FlagFam':>12} {'FamZ':>8} {'OrthoZ':>8} {'Flag':>6}\n"
    )
    output.write(f"  {'-'*25} {'-'*12} {'-'*8} {'-'*8} {'-'*6}\n")

    for sp in species_list:
        r = results[sp]
        flagfam_str = f"{r['flagged_count']}/{n_count_scoreable}"
        flag_str = " ***" if r["flagged"] else ""
        output.write(
            f"  {sp:<25} {flagfam_str:>12} {r['fam_z']:>8.3f} {r['ortho_z']:>8.3f} {flag_str:>6}\n"
        )

    span_flagged = [sp for sp in flagged_species if results[sp]["ortho_z"] > ortho_z]
    fam_flagged  = [sp for sp in flagged_species if results[sp]["fam_z"] > fam_z and results[sp]["ortho_z"] <= ortho_z]

    if span_flagged:
        output.write(f"\n  WARNING: {len(span_flagged)} species flagged as potential orthogroup inflators (ortho_z > {ortho_z}):\n")
        for sp in span_flagged:
            r = results[sp]
            output.write(f"    - {sp}: ortho_z={r['ortho_z']:.2f}\n")

    if fam_flagged:
        output.write(f"\n  WARNING: {len(fam_flagged)} species flagged for abnormal copy-number estimates (fam_z > {fam_z}):\n")
        for sp in fam_flagged:
            r = results[sp]
            output.write(f"    - {sp}: fam_z={r['fam_z']:.2f} ({r['flagged_count']}/{n_count_scoreable} families)\n")

    if flagged_species:
        output.write("\n  Consider re-running without these species.\n")
    else:
        output.write("\n  No species flagged.\n")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    app = argparse.ArgumentParser(
        description=(
            "Build orthogroups from TOGA2 pairwise orthology annotations.\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s [-h] -t DIR -s FILE -b FILE -i FILE -o DIR [options]",
        epilog=(
            "example:\n"
            "  toga2orthogroups.py \\\n"
            "    -t TOGA2 \\\n"
            "    -s species.lst \\\n"
            "    -b toga2.transcripts.bed \\\n"
            "    -i toga2.isoforms.tsv \\\n"
            "    -o orthogroups \\\n"
        ),
    )

    req = app.add_argument_group("required")
    req.add_argument(
        "-t", "--toga-dir",
        required=True,
        metavar="DIR",
        help="Directory with per-species TOGA2 output subdirs",
    )
    req.add_argument(
        "-s", "--species-list",
        required=True,
        metavar="LIST|FILE",
        help="Comma-separated species names or path to a file with one name per line",
    )
    req.add_argument(
        "-b", "--transcripts-bed",
        required=True,
        metavar="FILE",
        help="Reference transcript BED file",
    )
    req.add_argument(
        "-i", "--isoforms",
        required=True,
        metavar="FILE",
        help="TOGA2 isoforms file",
    )
    req.add_argument(
        "-o", "--out-dir",
        required=True,
        metavar="DIR",
        help="Output directory",
    )

    opt = app.add_argument_group("optional")
    opt.add_argument(
        "-f", "--force",
        action="store_true",
        help="Overwrite output files if they already exist",
    )
    opt.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print per-species processing stats",
    )
    opt.add_argument(
        "-ul", "--include-ul",
        action="store_true",
        help="Include UL (Uncertain Loss) transcripts",
    )
    opt.add_argument(
        "-sf", "--size-filter",
        type=int,
        metavar="INT",
        default=None,
        help="Remove gene families where any species has >= INT gene copies",
    )
    opt.add_argument(
        "-bl", "--blacklist",
        metavar="LIST|FILE",
        default=None,
        help=(
            "Comma-separated chr/scaffold names or path to a file with one name per line"
        ),
    )
    opt.add_argument(
        "--panther",
        metavar="FILE",
        default=None,
        help="PANTHER database flat file; if provided, enables PANTHER-guided family merging",
    )
    opt.add_argument(
        "--one-to-one",
        action="store_true",
        help="Per-reference gene matrix with query gene name if single-copy",
    )

    qc = app.add_argument_group("QC")
    qc.add_argument(
        "--ortho-z",
        type=float,
        metavar="FLOAT",
        default=3.0,
        help="OrthoZ threshold for orthogroup inflation detection",
    )
    qc.add_argument(
        "--fam-z",
        type=float,
        metavar="FLOAT",
        default=3.0,
        help="FamZ threshold for family copy-number outlier detection",
    )
    qc.add_argument(
        "--no-qc",
        action="store_true",
        help="Skip species QC diagnostics",
    )

    return app.parse_args(argv)


def _log_elapsed(t0: float) -> None:
    elapsed = time.perf_counter() - t0
    if elapsed < 60:
        log.info("\nTime elapsed: %.1fs", elapsed)
    else:
        m, s = divmod(int(elapsed), 60)
        log.info("\nTime elapsed: %dm %ds", m, s)


def run(
    toga_dir: str,
    species_list_path: str,
    transcripts_bed: str,
    isoforms: str,
    out_dir: str,
    force: bool = False,
    include_ul: bool = False,
    blacklist: str | None = None,
    panther: str | None = None,
    ortho_z: float = 3.0,
    fam_z: float = 3.0,
    no_qc: bool = False,
    one_to_one: bool = False,
    size_filter: int | None = None,
) -> None:
    """Core runner — called by both the argparse main() and the Click subcommand.

    Running: Filesystem I/O, orthogroup building, and QC.
    Logging is configured by the caller before invoking this function.
    """
    t0 = time.perf_counter()
    log.info("\n=== Running toga2orthogroups ===")

    # --- Resolve output directory and expected file paths ---
    _out_dir = Path(out_dir)
    _out_dir.mkdir(parents=True, exist_ok=True)

    out_matrix = _out_dir / "one2one_matrix.tsv"
    out_tsv    = _out_dir / "orthogroups_matrix.tsv"
    out_map    = _out_dir / "orthogroups_map.tsv"

    if not force:
        if one_to_one:
            if out_matrix.exists():
                log.error("Output file already exists: %s\nUse -f / --force to overwrite.", out_matrix)
                sys.exit(1)
        else:
            existing = [f for f in (out_tsv, out_map) if f.exists()]
            if existing:
                log.error(
                    "Output file(s) already exist: %s\nUse -f / --force to overwrite.",
                    ", ".join(str(f) for f in existing),
                )
                sys.exit(1)

    # --- Load species list ---
    species_list = _parse_species_list(species_list_path)
    log.info("Species: %d loaded", len(species_list))

    # --- Load reference gene set ---
    ref_genes = load_reference_genes(isoforms, transcripts_bed, blacklist=_parse_blacklist(blacklist))

    # --- One-to-one mode: per-reference-gene matrix, bypasses Union-Find ---
    if one_to_one:
        if panther:
            log.warning("--panther is ignored in --one-to-one mode")
        if size_filter is not None:
            log.warning("--size-filter is ignored in --one-to-one mode")

        # Collect raw per-(ref_gene, species) individual query gene sets.
        per_gene: dict[str, dict[str, set[str]]] = {rg: {} for rg in ref_genes.genes}
        for spe in species_list:
            log.debug("Processing species: %s", spe)
            sp_ortho = load_species_orthologs(spe, toga_dir, ref_genes.genes, include_ul=include_ul)
            for rg, q_genes in sp_ortho.ref_to_query.items():
                # q_gene strings may be comma-separated lists (many2many rows);
                # split to collect distinct individual query gene loci.
                per_gene[rg][spe] = {g for qg in q_genes for g in qg.split(",")}

        header = ["Ref_ID"] + species_list
        rows = []
        for rg in sorted(per_gene):
            sp_genes = per_gene[rg]
            n_ones = sum(1 for sp in species_list if len(sp_genes.get(sp, set())) == 1)
            if n_ones < 2:
                continue
            row: list = [rg]
            for sp in species_list:
                genes = sp_genes.get(sp, set())
                row.append(next(iter(genes)) if len(genes) == 1 else "-")
            rows.append(row)

        with open(out_matrix, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(header)
            writer.writerows(rows)
        log.info("Wrote one2one matrix: %s (%d genes)", out_matrix, len(rows))
        _log_elapsed(t0)
        return

    # --- Build orthogroups ---
    if panther:
        log.info("Running in PANTHER mode...")
        orthogroups = build_orthogroups(
            ref_genes, species_list, toga_dir,
            panther_path=panther,
            isoforms_path=isoforms,
            include_ul=include_ul,
        )
    else:
        orthogroups = build_orthogroups(
            ref_genes, species_list, toga_dir,
            include_ul=include_ul,
        )

    header, rows = generate_count_table(orthogroups, species_list)

    if size_filter is not None:
        rows, n_removed = apply_size_filter(rows, size_filter)
        log.info("Size filter (threshold=%d): %d families removed, %d kept", size_filter, n_removed, len(rows))
        kept_ids = {row[1] for row in rows}
        orthogroups.families = {fid: m for fid, m in orthogroups.families.items() if fid in kept_ids}

    write_count_table(header, rows, out_tsv)
    write_orthogroup_membership(orthogroups, out_map)

    if not no_qc:
        species_qc_diagnostics(
            orthogroups, species_list,
            ortho_z=ortho_z,
            fam_z=fam_z,
        )

    _log_elapsed(t0)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    # --- Configure logging ---
    # With --verbose: DEBUG log level with per-species stats.
    logging.basicConfig(format="%(message)s", stream=sys.stderr)
    logging.getLogger().setLevel(logging.DEBUG if args.verbose else logging.INFO)

    run(
        toga_dir=args.toga_dir,
        species_list_path=args.species_list,
        transcripts_bed=args.transcripts_bed,
        isoforms=args.isoforms,
        out_dir=args.out_dir,
        force=args.force,
        include_ul=args.include_ul,
        blacklist=args.blacklist,
        panther=args.panther,
        ortho_z=args.ortho_z,
        fam_z=args.fam_z,
        no_qc=args.no_qc,
        one_to_one=args.one_to_one,
        size_filter=args.size_filter,
    )


if __name__ == "__main__":
    main()
