#!/usr/bin/env python3
"""
Assembles projections from fragments scattered across different contigs in the query assembly
"""

# import os
# import sys

# LOCATION: str = os.path.dirname(os.path.abspath(__file__))
# PARENT: str = os.sep.join(LOCATION.split(os.sep)[:-1])
# sys.path.extend([LOCATION, PARENT])

import sys
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from os import PathLike
from typing import Any, TextIO

import click

from .constants import Headers
from .shared import (
    CONTEXT_SETTINGS,
    CommandLineManager,
    intersection,
    read_tab,
)  ## TODO: most of legacy T1 code in this import line is no longer needed; remove in the stable version

# artificial 0-scored points
SOURCE = "SOURCE"
SINK = "SINK"
SCORE_THRESHOLD = 0.5
EXON_COV_THRESHOLD = 1.33
MAX_OVERLAP = 250  # TODO: check whether 250 is a plausible value
CHAIN: str = "chain"
TRANSCRIPT: str = "transcript"

__author__ = "Ekaterina Osipova & Bogdan M. Kirilenko"
__credits__ = "Yury V. Malovichko"
__year__ = "2024"
__all__ = ()


def all_unique(collection: Iterable[Any]) -> bool:
    """
    Check if all members in a collection are unique.
    Based on solution from https://stackoverflow.com/questions/5278122

    Args:
        collection: an iterable collection of arbitrary elements

    Returns:
        bool: True if all elements are unique, False otherwise
    """
    ## TODO: Move to src.python.modules.shared before the release
    seen: set[Any] = set()
    return not any(x in seen or seen.add(x) for x in collection)


@dataclass
class ChainData:
    __slots__ = ("chrom", "end", "id", "score", "start")
    id: str
    chrom: str
    start: int
    end: int
    score: float


class Vertex:
    __slots__ = ("children", "chrom", "end", "name", "score", "start")
    """Graph vertex."""

    def __init__(self, name: str, chrom: str, start: int, end: int, score: float) -> None:
        self.name: str = name
        self.chrom: str = chrom
        self.start: int = start
        self.end: int = end
        self.score: float = score
        self.children: list[str] = []

    def add_child(self, v):
        if v not in self.children:
            self.children.append(v)


class Graph:
    """Build a directed graph using adjacency list."""

    def __init__(self):
        self.vertices = {}

    def add_vertex(self, vertex):
        """Add vertex if it's not in the vertices list."""
        if isinstance(vertex, Vertex) and vertex.name not in self.vertices:
            self.vertices[vertex.name] = vertex
            return True
        else:
            return False

    def add_edge(self, parent, child):
        """Find vertices with parent and child names."""
        if parent in self.vertices and child in self.vertices:
            self.vertices[parent].add_child(child)
            return True
        else:
            return False

    def topological_sort_util(self, v, visited, stack):
        """Perform Depth First Search

        Mark the current node as visited.
        """
        visited[v] = True
        # check all children of this vertex if they're visited
        for i in self.vertices[v].children:
            if visited[i] is False:
                self.topological_sort_util(i, visited, stack)
        # add current vertex to stack
        stack.insert(0, v)

    def topological_sort(self) -> list[int]:
        """Perform topological sort.

        Use recursive function topological_sort_util().
        Mark all the vertices as not visited.
        """
        visited = {v: False for v in self.vertices}
        # initiate stack to store sorted vertices
        stack = []
        for vertex in self.vertices:
            if visited[vertex] is False:
                self.topological_sort_util(vertex, visited, stack)
        # return sorted list of vertices
        return stack

    def add_source_sink_graph(self) -> None:
        """Add artificial Source and Sink vertices to the chain graph.

        Assign them zero length and zero score.
        """
        source_end = min([self.vertices[vertex].start for vertex in self.vertices])
        source_start = source_end
        sink_start = max([self.vertices[vertex].end for vertex in self.vertices])
        sink_end = sink_start
        self.add_vertex(Vertex(SOURCE, "", source_start, source_end, 0))
        self.add_vertex(Vertex(SINK, "", sink_start, sink_end, 0))

        ## add edges from source to each vertex
        for vertex in self.vertices:
            if vertex != SOURCE:
                self.add_edge(SOURCE, vertex)

        ## add edges from each vertex to sink
        for vertex in self.vertices:
            if vertex != SINK:
                self.add_edge(vertex, SINK)

    def find_shortest_path(
        self, 
        sorted_vertices: list[int], 
        diff_contigs_only: bool
    ) -> tuple[float, list[int]]:
        """Find the shortest path in directed acyclic graph.

        Initiate dictionary with the shortest paths to each node:
        {vertex: (value, path itself)}.
        """
        shortest_paths = {}
        for sorted_vertex in sorted_vertices:
            shortest_paths[sorted_vertex] = (0, [SOURCE])

        ## check each child of the current vertex
        ## and update shortest path to this vertex in the dictionary
        for sorted_vertex in sorted_vertices:
            children = self.vertices[sorted_vertex].children
            for child in children:
                current_score: float = shortest_paths[child][0]
                path_score: float = shortest_paths[sorted_vertex][0]
                path: list[str] = shortest_paths[sorted_vertex][1]
                child_score = self.vertices[child].score
                upd_score = path_score + child_score
                if diff_contigs_only:
                    all_chroms: list[str] = [self.vertices[x].chrom for x in path]
                    all_chroms.append(self.vertices[sorted_vertex].chrom)
                    all_chroms.append(self.vertices[child].chrom)
                    all_chroms = [x for x in all_chroms if x]
                    all_chroms_unique: bool = all_unique(all_chroms)
                    # print(f"{all_chroms=}, {all_chroms_unique=}")
                else:
                    all_chroms_unique: bool = True
                ## update the current path to the child
                ## if its score is lower than current recorded one
                if upd_score < current_score and all_chroms_unique:
                    new_path = list(shortest_paths[sorted_vertex][1])
                    if sorted_vertex not in new_path:
                        new_path.append(sorted_vertex)
                    shortest_paths[child] = (upd_score, new_path)
        return shortest_paths[SINK]

    def __repr__(self) -> str:
        lines = []
        for elem in self.vertices:
            line = f"{elem}\t{self.vertices[elem].children}\n"
            lines.append(line)
        return "".join(lines)


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument("chain_file", type=click.Path(exists=True), metavar="CHAIN_FILE")
@click.argument("orthology_scores", type=click.Path(exists=True), metavar="ORTH_SCORES")
@click.argument("bed_file", type=click.Path(exists=True), metavar="REF_BED")
@click.option(
    "--output",
    "-o",
    type=click.File("w", lazy=True),
    metavar="FILE",
    default=sys.stdout,
    show_default=False,
    help="A path to an output file [default: stdout]",
)
@click.option(
    "--fragmented_only",
    "-f",
    is_flag=True,
    default=False,
    show_default=True,
    help=("If set, only fragmented (multi-chain) projections are listed in the output"),
)
@click.option(
    "--orthology_threshold",
    "-ot",
    type=float,
    metavar="FLOAT",
    default=SCORE_THRESHOLD,
    show_default=True,
    help="Probability threshold for considering projections as orthologous",
)
@click.option(
    "--diff_contigs_only",
    is_flag=True,
    default=False,
    show_default=True,
    help=(
        "If set, stitching fragments located on the same query contig "
        "(scaffold, chromosome, etc.) is deprecated"
    ),
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
    "--verbose", "-v", is_flag=True, default=False, help="Control logging verbosity"
)
@click.option(
    "--debug",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Increase logging verbosity by logging debugging messages. 
Automatically enables \"--verbose\" flag""",
)
class FragmentStitcher(CommandLineManager):
    """ """

    __slots__ = (
        "bed_file",
        "chain_coordinates",
        "chain_file",
        "diff_contigs_only",
        "fragmented_only",
        "orth_chains",
        "orth_threshold",
        "output",
        "score_file",
        "tr2chains2score",
        "tr2exons",
        "tr2fragments",
    )

    def __init__(
        self,
        chain_file: str | PathLike,
        orthology_scores: TextIO,
        bed_file: TextIO,
        output: TextIO,
        fragmented_only: bool,
        orthology_threshold: float,
        diff_contigs_only: bool,
        log_name: str | None,
        verbose: bool,
        debug: bool,
    ) -> None:
        self.v: bool = verbose
        self.debug: bool = debug
        self.set_logging(name=log_name, toga_module="stitch_fragments")
        self.chain_file: str | PathLike = chain_file
        self.score_file: TextIO = orthology_scores
        self.bed_file: TextIO = bed_file
        self.output: TextIO = output
        self.orth_threshold = orthology_threshold
        self.fragmented_only: bool = fragmented_only
        self.diff_contigs_only: bool = diff_contigs_only
        ## all orthologous chain storage
        self.orth_chains: set[str] = set()
        ## chain coordinates
        self.chain_coordinates: dict[str, tuple[str, int, int]] = {}
        ## transcript-to-chain mapping
        self.tr2exons: dict[str, list[tuple[int, ...]]] = {}
        ## input transcript:chains mapping: orthologous chain lists
        self.tr2chains2score: dict[str, dict[str, float]] = defaultdict(dict)
        ## output transcript:chains mapping:
        self.tr2fragments: dict[str, list[str]] = {}
        self.output: TextIO = output

        self.run()

    def run(self) -> None:
        """Main method"""
        self.read_scores()
        self.read_bed()
        self.read_chain_file()
        self.stitch_fragments()

    def read_scores(self) -> None:
        for i, data in enumerate(read_tab(self.score_file), start=1):
            if len(data) < 3:
                self._die(
                    "Improper score file formatting at line %i: expected at least three fields, got %i"
                    % (i, len(data))
                )
            if data[0] == TRANSCRIPT:
                continue
            tr, chain, score = data
            try:
                score: float = float(score)
            except ValueError:
                self._die(
                    "Improper score file formatting at line %i: orthology score %s is not a valid decimal"
                    % (i, score)
                )
            if score < self.orth_threshold:
                continue
            self.tr2chains2score[tr][chain] = score
            self.orth_chains.add(chain)

    def read_chain_file(self) -> None:
        ## TODO: A good target for Rust/Maturin implementation
        ## the following values are required: chromosome, start, end, strand
        for i, data in enumerate(read_tab(self.chain_file, sep=" "), start=1):
            if data[0] != CHAIN:
                continue
            ## the fields needed are: start and end in **reference**
            ## and chromosome in **query**
            if len(data) != 13:
                self._die(
                    (
                        "Wrong chain header formatting in chain file %s at line %i; "
                        "expected 13 fields, got %i"
                    )
                    % (self.chain_file, i, len(data))
                )
            chain_id: str = data[12]
            ## do not bother with chains that are not orthologous to any known ref locus
            if chain_id not in self.orth_chains:
                continue
            if data[4] not in ("+", "-"):
                self._die(
                    (
                        "Wrong chain header formatting in chain file %s at line %i; "
                        "expected either '+' or '-' for strand symbol, got %s"
                    )
                    % (self.chain_file, i, data[4])
                )
            strand: bool = data[4] == "+"
            if strand: 
                try:
                    ref_start: int = int(data[5])
                    ref_end: int = int(data[6])
                except ValueError:
                    self._die(
                        (
                            "Wrong chain header formatting in chain file %s at line %i; "
                            "reference coordinates are not numeric values"
                        )
                        % (self.chain_file, i)
                    )
            else:
                try:
                    ref_length: int = int(data[3])
                    ref_start: int = int(ref_length) - int(data[6])
                    ref_end: int = int(ref_length) - int(data[5])
                except ValueError:
                    self._die(
                        (
                            "Wrong chain header formatting in chain file %s at line %i; "
                            "reference coordinates are not numeric values"
                        )
                        % (self.chain_file, i)
                    )
            ## then, get the **query** chromosome
            query_chrom: str = data[7]
            ## record the retrieved data in self.chain_coordinates
            self.chain_coordinates[chain_id] = (query_chrom, ref_start, ref_end)

    def read_bed(self) -> dict[str, list[tuple[int, ...]]]:
        for i, data in enumerate(read_tab(self.bed_file), start=1):
            start: int = int(data[1])
            name = data[3]
            if name not in self.tr2chains2score:
                continue
            cds_start: int = int(data[6])
            cds_end: int = int(data[7])
            strand: bool = data[5] == "+"
            sizes = [int(x) for x in data[10].split(",") if x != ""]
            starts = [int(x) for x in data[11].split(",") if x != ""]
            exon_num: int = len(sizes)
            exon_coords: list[tuple[int, ...]] = []
            for i in range(exon_num):
                exon_start: int = start + starts[i]
                if exon_start >= cds_end:
                    break
                exon_end: int = exon_start + sizes[i]
                if exon_end <= cds_start:
                    continue
                exon_start = max(exon_start, cds_start)
                exon_end = min(exon_end, cds_end)
                num: int = i + 1 if strand else exon_num - i
                exon_coords.append((num, exon_start, exon_end))
            self.tr2exons[name] = exon_coords

    def stitch_fragments(self) -> None:
        """
        Main fragment stitching method:
        For each transcript, assess exon-to-chain intersection for all coding exons, construct
        an ordered chain graph and sort it
        """
        header_written: bool = False
        for tr, chains2scores in self.tr2chains2score.items():
            exons: list[tuple[int, int]] | None = self.tr2exons.get(tr)
            if exons is None:
                self._die("No exon records found for reference transcript %s" % tr)
            chain2exons: dict[str, tuple[int, ...]] = self._exon_cov_by_chain(
                exons, chains2scores
            )
            if len(chain2exons) <= 1:
                continue
            ## if any chain covers the entirety of transcript,
            ## there is no need to assemble fragmented copies
            # if any(all(y) for x, y in chain2exons.items()):
            if any(all(x in y for x in range(1, len(exons) + 1)) for y in chain2exons.values()):
                continue
            ## next, check exon-to-chain intersection for all exon ("exon coverage")
            avg_exon_cov: float = self._get_exon_coverage(chain2exons, len(exons))
            # print(f"AFTER: {tr=}, {chains2scores=}, {avg_exon_cov=}")
            ## do not stitch fragments with high average exon-to-chain intersection rate
            if avg_exon_cov > EXON_COV_THRESHOLD:
                continue
            ## otherwise, proceed with chain graph construction
            chain_graph: Graph = self._build_graph(chains2scores)
            ## add dummy source and sink nodes
            chain_graph.add_source_sink_graph()
            ## perform topological sort
            sorted_vertices: list[int] = chain_graph.topological_sort()
            ## find the shortest path in the graph;
            ## this be the fragment projection configuration
            _, shortest_path = chain_graph.find_shortest_path(sorted_vertices, self.diff_contigs_only)
            # if tr == "XM_011532850.3#RGPD1":
            #     print(f"{tr=}, {chain2exons=}, {avg_exon_cov=}, {sorted_vertices=}, {shortest_path=}")
            #     for chain in chain2exons:
            #         print(f"{chain=}, {self.chain_coordinates[chain]=}")
            ## remove SINK from the path
            # print(f"{shortest_path=}")
            shortest_path = shortest_path[1:]
            ## do not report single-chain transcripts
            if self.fragmented_only and len(shortest_path) < 2:
                continue
            self.tr2fragments[tr] = shortest_path
            del chain_graph
            if not header_written:
                self.output.write(Headers.FRAGM_PROJ_HEADER)
                header_written = True
            chain_str: str = ",".join(map(str, sorted(map(int, shortest_path))))
            self.output.write(f"{tr}\t{chain_str}\n")

    def _exon_cov_by_chain(
        self, exons: list[tuple[int, ...]], chains2scores: list[str]
    ) -> dict[str, tuple[int, ...]]:
        """ """
        chains: list[str] = sorted(
            chains2scores.keys(), key=lambda x: self.chain_coordinates[x][1]
        )
        output: dict[str, tuple[int, ...]] = {}
        # curr_exon: int = 0
        for chain in chains:
            covered_exons: list[str] = []
            start: int = self.chain_coordinates[chain][1]
            end: int = self.chain_coordinates[chain][2]
            for exon, exon_start, exon_end in exons:
                if intersection(start, end, exon_start, exon_end) > 0:
                    # if exon_start >= start and exon_end <= end:
                    covered_exons.append(exon)
            output[chain] = tuple(covered_exons)
        return output

    def _get_exon_coverage(
        self, chain2exons: dict[str, tuple[int, ...]], exon_num: int
    ) -> float:
        """
        Calculates the average number of exon occurrences within chain boundaries

        Args:
            chain2exons: {chain: {exon_num: bool}} dictionary of exon occurrence within each chain
            exon_num: number of exons in the focal transcript

        Returns:
            float average number of exon occurrences
        """
        exon_cov: int = sum([len(x) for x in chain2exons.values()])
        average_cov = exon_cov / exon_num
        return average_cov

    def _build_graph(self, chains2scores: dict[str, float]) -> Graph:
        """
        Constructs a Graph object with chains overlapping the transcript as vertices

        Args:
            chains2scores: {chain_id: orthology_score} dictionary

        Returns:
            A Graph object for the focal transcript
        """
        chain_graph: Graph = Graph()
        # add all vertices to the chain graph
        for chain_id, score in chains2scores.items():
            ## safeguard check for coordinate presence in the collection
            if chain_id not in self.chain_coordinates:
                self._die("Chain %s has no recorded coordinates" % chain_id)
            chrom, start, end = self.chain_coordinates.get(chain_id)
            v: Vertex = Vertex(chain_id, chrom, start, end, -1 * score)
            chain_graph.add_vertex(v)

        # add edges to the chain graph
        for i in chain_graph.vertices:
            for j in chain_graph.vertices:
                if i == j:
                    ## ignore self-edges
                    continue
                i_vertex: Vertex = chain_graph.vertices[i]
                j_vertex: Vertex = chain_graph.vertices[j]
                ## do not connect fragments located on the same query chromosome
                if i_vertex.chrom == j_vertex.chrom and self.diff_contigs_only:
                    continue
                ## allow some overlap between chains
                if (
                    intersection(
                        i_vertex.start, i_vertex.end, j_vertex.start, j_vertex.end
                    )
                    > MAX_OVERLAP
                ):
                    continue
                chain_graph.add_edge(i_vertex.name, j_vertex.name)
        return chain_graph


if __name__ == "__main__":
    FragmentStitcher()
