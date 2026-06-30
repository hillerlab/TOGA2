#!/usr/bin/env python3

__author__ = "Tim Stadager"
__credits__ = ("Leon Hilgers", "Yury V. Malovichko")
__year__ = "2026"

from collections import defaultdict
from .shared import CommandLineManager, hex_dir_name, read_tab
from typing import Dict, List, Optional, Set, Tuple, Union

import os

CHROM_SIZES: str ="chrom.sizes"
ORTH_FILE: str = "orthology_classification.tsv"
TOGA2: str = "TOGA2"
ONE2ONE: str = "one2one"
QUERY_GENE_BED: str = "query_genes.bed"
REF_GENES: str = "{ref}.toga.genes.bed"
MODULE_NAME: str = "toga2agora"
GENES: str = "genes"
ORTH_GROUPS: str = "orthology_groups"
GENE_FILE: str = "genes.{species}.list"
ANC_ORTH_FILE: str = "orthologyGroups.{node}.list"
GENE_FILE_LINE: str = "{chrom}\t{start}\t{end}\t{strand}\t{gene}\n"


class TreeNode:

    __slots__ = ("name", "ancestor", "children", "chromosomes")

    def __init__(
        self, 
        name: str, 
        ancestor: Optional[str] = None,
        children: Optional[List[str]] = None,
        chromosomes: Optional[List[str]] = None,
    ) -> None:
        self.name: str = name
        self.ancestor: Union[str, None] = ancestor
        self.children: List[str] = [] if children is None else children
        self.chromosomes: Dict[str, str] = {} if chromosomes is None else chromosomes

    def update_ancestor(self, ancestor: str) -> None:
        self.ancestor = ancestor

    def update_chromosomes(self, chromosomes: Dict[str, str]) -> None:
        self.chromosomes = chromosomes


class Leaf(TreeNode):

    __slots__ = ("genes",)

    def __init__(self, *args, genes: Optional[Dict[str, List[Tuple[str, ...]]]] = None,) -> None:
        super().__init__(*args)
        self.genes: Dict[str, List[str, bool]] = genes if genes is not None else {}

    def update_genes(self, genes: Dict[str, List[Tuple[str, ...]]]) -> None:
        self.genes = genes

    def all_genes(self) -> Set[str]:
        out_set: Set[str] = set()
        for gene_list in self.genes.values():
            out_set.update((x[2] for x in gene_list))
        return out_set


class InternalNode(TreeNode):

    __slots__ = ("genes",)

    def __init__(self, *args, genes: Optional[Set[str]] = None) -> None:
        super().__init__(*args)
        self.genes: Set[str] = genes if genes is not None else set()

    def update_genes(self, genes: Set[str]) -> None:
        if not self.genes:
            self.genes = genes
        else:
            self.genes = self.genes.intersection(genes)

    def all_genes(self) -> Set[str]:
        return self.genes


class Toga2Agora(CommandLineManager):
    """
    \b
    MMP""MM""YMM   .g8""8q.    .g8\"""bgd      db          `7MMF'`7MMF'
    P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM
         MM     dM'      `MM dM'       `     ,V^MM.          MM    MM
         MM     MM        MM MM             ,M  `MM          MM    MM
         MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM
         MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM
       .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

    \b
    toga2agora - Prepare AGORA input from multiple TOGA2 pairwise annotation
    WARNING: This module is currently under development. Alternative input configuration options might be added in the future

    \b
    Mandatory arguments are:
    *INPUT_DIR* is the path to the master directory of the following structure:
        * reference
            * TOGA2
                * *ANNOTATION*
                    * ${REF_NAME}.toga.genes.bed
                * vs_query1
                * vs_query2
                * ...
                * vs_queryN
        * query1
            * chrom.sizes
        * query2
            * chrom.sizes
        * ...
        * queryN
            * chrom.sizes
    *TREE* is a path to the phylogenetic tree in Newick format. All leaves are expected to be named after the reference/query directories in *INPUT_DIR*
    *REF_NAME* is reference species names; the same name is expected to be used in both the tree and the *INPUT_DIR* structure.

    \b
    Output directory has the following structure:
    * `genes` - a directory with ordered gene names lists in AGORA-compatible format (chromosome-start-end-strand-gene);
    * `orthology_groups` - a directory with ancestral nodes' gene lists;
    * `AGORA_tree.nwk` - a labeled version of the input tree; each internal node is named after its immediate descendants
    """

    __slots__ = (
        "input_dir",
        "tree",
        "reference",
        "min_chrom_size",
        "output",
        "ref_path",
        "ref_annot_path",
        "ref_toga_path",
        "gene_dir",
        "orth_dir",
        "v",
        "query2chrom_size",
        "agora_tree_path",
        "species",
        "root_name",
        "tree_nodes",
    )

    def __init__(
        self,
        input_directory: os.PathLike,
        tree: os.PathLike,
        reference: str,
        annotation: Optional[os.PathLike],
        min_chrom_size: Optional[int] = 0,
        output: Optional[os.PathLike] = None,
        verbose: Optional[bool] = False,
        debug: Optional[bool] = False,
    ) -> None:
        self.debug: bool = debug
        self.v: bool = verbose | debug
        self.set_logging(name=MODULE_NAME, toga_module=MODULE_NAME)
        self.input_dir: os.PathLike = input_directory
        self.reference: str = reference
        self.ref_path: str = os.path.join(self.input_dir, self.reference)
        if not os.path.exists(self.ref_path):
            self._die(
                "Reference directory %s does not exist; exiting" % self.ref_path
            )
        self.ref_toga_path: str = os.path.join(self.ref_path, TOGA2)
        if not os.path.exists(self.ref_toga_path):
            self._die(
                "Reference TOGA2 directory %s does not exist" % self.ref_toga_path
            )
        self.ref_annot_path: str = os.path.join(self.ref_toga_path, annotation)
        if not os.path.exists(self.ref_annot_path):
            self._die(
                "Reference annotation directory %s does not exist" % self.ref_annot_path
            )
        
        self.tree: os.PathLike = tree
        self.output: os.PathLike = output if output is not None else hex_dir_name(MODULE_NAME)
        self.min_chrom_size: int = min_chrom_size
        self.agora_tree_path: str = os.path.join(self.output, "AGORA_tree.nwk")
        self.gene_dir: str = os.path.join(self.output, GENES)
        self.orth_dir: str = os.path.join(self.output, ORTH_GROUPS)

        self.query2chrom_size: Dict[str, os.PathLike] =  {}
        self.tree_nodes: Dict[str, TreeNode] = {}
        self.species: Set[str] = set()
        self.root_name: Union[str, None] = None

        self.run()

    def run(self) -> None:
        """Create output directories, process the tree and all species, propagate ancestral gene sets, and write output."""
        ## create the input directory
        self._to_log("Create output directory")
        self._mkdir(self.output)
        self._mkdir(self.gene_dir)
        self._mkdir(self.orth_dir)
        ## process the input phylogenetic tree; extract node names, name internal nodes if necessary
        self.process_tree()
        ## process query data
        for query in self.species:
            if query == self.reference:
                self._to_log("Processing data for reference species %s" % self.reference)
                self.process_reference()
                continue
            self._to_log("Processing data for query %s" % query)
            self.process_query(query)
        ## create ortholog lists for ancestors
        self.propagate_ancestral_orthologs()
        self.write_output()
        self._to_log("toga2agora finished successfully")

    def process_tree(self) -> None:
        """
        Processes the input phylogenetic tree as follows:
        * Checks if internal nodes in the tree are labeled:
        *     if internal nodes are labeled, saves their names;
        *     otherwise, assigns the internal nodes and saves them;
        * Saves children nodes for each internal node;
        * Writes/copies the updated input tree to the output directory

        Args:
            None

        Returns:
            None
        """
        try:
            from ete3 import Tree
        except ImportError:
            self._die("ETE3 is required for tree processing: pip install ete3")
            return

        # format=0: flexible Newick that stores internal labels as node.name
        t = Tree(str(self.tree), format=0)

        ## traverse the tree starting from the leaf nodes;
        ## save the node names and their immediate ancestors, name the ancestral node if needed
        ancestral_nodes: int = 0
        for node in t.traverse("postorder"):
            ## TODO: Input parsing, gene list propagation, and output writing can be safely handled here
            if node.is_leaf():
                node_object: TreeNode = Leaf(node.name)
                self.tree_nodes[node.name] = node_object
                self.species.add(node.name)
            else:
                if not (node.name and node.name.strip()):
                    node.name = "_".join(sorted(node.get_leaf_names()))
                name = node.name
                node_object: TreeNode = InternalNode(node.name)
                self._debug(
                    "Internal node %s has the following children nodes: %s" 
                    % (name, ",".join([x.name for x in node.get_children()]))
                )
                for child in node.get_children():
                    if child.name not in self.tree_nodes:
                        self._die("Child node %s has been skipped for some reason" % child)
                    self.tree_nodes[child.name].update_ancestor(name)
                    node_object.children.append(child.name)
                self.tree_nodes[node.name] = node_object
                ancestral_nodes += 1
                if node.is_root():
                    self.root_name = name

        # Write the labelled tree to AGORA_tree.nwk
        t.write(format=0, outfile=self.agora_tree_path)

        self._to_log(
            "Tree processed: %d species (reference + %d queries), %d ancestral nodes"
            % (len(self.species), len(self.species) - 1, ancestral_nodes)
        )

    def process_reference(self) -> None:
        """Parse chromosome sizes and the reference gene BED file; populate genes for the reference TreeNode."""
        ref_path: str = os.path.join(self.input_dir, self.reference)
        if not os.path.exists(ref_path):
            self._die("Reference directory %s does not exist" % ref_path)
        chrom_sizes: str = os.path.join(ref_path, CHROM_SIZES)
        if not os.path.exists(chrom_sizes):
            self._die(
                "Chromosome file %s for the reference genome %s does not exist" 
                % (chrom_sizes, self.reference)
            )
        mapping, chrom_order = self._parse_chromosomes(chrom_sizes)
        self.tree_nodes[self.reference].update_chromosomes(chrom_order)
        if not os.path.exists(self.ref_annot_path):
            self._die(
                "Annotation directory %s for the reference assembly %s does not exist"
                % (self.ref_annot_path, self.reference)
            )
        ref_gene_file: str = os.path.join(self.ref_annot_path, REF_GENES.format(ref=self.reference))
        if not os.path.exists(ref_gene_file):
            self._die(
                "Reference gene BED file %s does not exist" % ref_gene_file
            )
        chrom2items: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)
        for data in read_tab(ref_gene_file):
            chrom: str = data[0]
            new_name: Union[str, None] = mapping.get(chrom)
            if new_name is None:
                continue
            start: int = int(data[1])
            end: int = int(data[2])
            name: str = data[3]
            strand: str = data[5] + "1"
            chrom2items[new_name].append((start, end, name, strand))
        for chrom in chrom2items:
            chrom2items[chrom].sort(key=lambda x: (x[0], x[1]))
        self.tree_nodes[self.reference].update_genes(chrom2items)
        self._to_log(
            "Reference species %s has %i genes" 
            % (self.reference, len(self.tree_nodes[self.reference].all_genes()))
        )

    def process_query(self, query: str) -> None:
        """Parse chromosome sizes, extract 1:1 orthologs from the TOGA2 output, and populate genes for the query TreeNode."""
        ## check that the chrom.sizes file exists
        query_path: str = os.path.join(self.input_dir, query)
        if not os.path.exists(query_path):
            self._die(
                "Directory %s does not exist" % query_path 
            )
        chrom_sizes: str = os.path.join(query_path, CHROM_SIZES)
        if not os.path.exists(chrom_sizes):
            self._die("Chromosome size file %s does not exist" % chrom_sizes)
        mapping, chrom_order = self._parse_chromosomes(chrom_sizes)
        self.tree_nodes[query].update_chromosomes(chrom_order)
        ## now, check whether the TOGA2 output directory exists and contains a orthology classification file
        toga_path: str = os.path.join(self.ref_toga_path, "vs_" + query)
        if not os.path.exists(toga_path):
            self._die(
                "TOGA2 output directory for query %s (%s) does not exist" % (
                    query, toga_path
                )
            )
        orth_file: str = os.path.join(toga_path, ORTH_FILE)
        if not os.path.exists(orth_file):
            self._die(
                "Orthology classification file %s does not exist for query %s" % (
                    orth_file, query
                )
            )
        ## extract the 1:1 orthologs
        query_genes: Set[str] = set()
        for data in read_tab(orth_file):
            if data[4] != ONE2ONE:
                continue
            query_genes.add(data[2])
        ## check if the query bed file exists
        chrom2items: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)
        query_gene_file: str = os.path.join(toga_path, QUERY_GENE_BED)
        if not os.path.exists(query_gene_file):
            self._die(
                "Query gene file %s does not exist" % query_gene_file
            )
        for data in read_tab(query_gene_file):
            name: str = data[3]
            if name not in query_genes:
                continue
            chrom: str = data[0]
            new_name: Union[str, None] = mapping.get(chrom)
            if new_name is None:
                continue
            start: int = int(data[1])
            end: int = int(data[2])
            strand: str = data[5] + "1"
            chrom2items[new_name].append((start, end, name, strand))
        ## add the gene data to the respective TreeNode object
        for chrom in chrom2items:
            chrom2items[chrom] = sorted(chrom2items[chrom], key=lambda x: (x[0], x[1]))
        self.tree_nodes[query].update_genes(chrom2items)
        self._to_log(
            "Query %s has %i genes"
            % (query, len(self.tree_nodes[query].all_genes()))
        )

    def propagate_ancestral_orthologs(self) -> None:
        """
        Infers shared gene sets for common ancestors in the tree.
        The tree is traversed in postorder (leaves-to-root); for each ancestral node,
        its gene set is defined as intersection of all its children's gene sets.

        Args:
            None

        Returns:
            None
        """
        stack: List[str] = []
        root_index: int = 0
        root_name: Union[str, None] = self.root_name
        while root_name is not None or stack:
            ## traverse towards the leaves until a leaf node is encountered
            if root_name is not None:
                stack.append((root_name, root_index))
                root_index = 0
                if len(self.tree_nodes[root_name].children):
                    root_name = self.tree_nodes[root_name].children[0]
                else:
                    root_name = None
                continue
            ## current node is a leaf; pop it
            curr_node_name, curr_node_idx = stack.pop()
            ## update the parent's gene list
            parent: Union[str, None] = self.tree_nodes[curr_node_name].ancestor
            if parent is None:
                break
            self._debug("Current node: %s, parent node: %s" % (curr_node_name, parent))
            self.tree_nodes[parent].update_genes(self.tree_nodes[curr_node_name].all_genes())
            ## recursively remove last children of the respective parent from the stack
            while stack and curr_node_idx == len(self.tree_nodes[parent].children) - 1:
                curr_node_name, curr_node_idx = stack.pop()
                ## and again, update the gene set for the last ancestor
                parent: Union[str, None] = self.tree_nodes[curr_node_name].ancestor
                if parent is None:
                    break
                self._debug("Current node: %s, parent node: %s" % (curr_node_name, parent))
                self.tree_nodes[parent].update_genes(self.tree_nodes[curr_node_name].all_genes())

            ## proceed to the firs unvisited child of the stack's last element
            ## leaves and nodes with all children visited are guaranteed to have been 
            ## removed from the stack beforehand
            if stack:
                root_name = self.tree_nodes[stack[-1][0]].children[curr_node_idx + 1]
                root_index = curr_node_idx + 1
    
        if self.debug:
            for node in self.tree_nodes:
                self._debug(
                        "Node %s has %i genes in the gene set" 
                        % (node, len(self.tree_nodes[node].all_genes()))
                    )

    def write_output(self) -> None:
        """
        Writes ordered gene coordinates for leaf nodes at genes/ and 
        ancestral gene lists for internal nodes at orthology_groups.

        Args:
            None
        
        Returns:
            None
        """
        for name, node in self.tree_nodes.items():
            if name in self.species:
                output_file: str = os.path.join(self.gene_dir, GENE_FILE.format(species=name))
                with open(output_file, "w") as h:
                    for chrom in node.chromosomes:
                        for gene in node.genes[chrom]:
                            out_line: str = GENE_FILE_LINE.format(
                                chrom=chrom,
                                start=gene[0],
                                end=gene[1],
                                strand=gene[3],
                                gene=gene[2],
                            )
                            h.write(out_line)
            else:
                output_file: str = os.path.join(self.orth_dir, ANC_ORTH_FILE.format(node=name))
                with open(output_file, "w") as h:
                    for gene in node.all_genes():
                        h.write(gene + "\n")

    def _parse_chromosomes(self, file: str) -> Tuple[Dict[str, str], List[str]]:
        """
        Parses the chromosome file, filters and sorts chromosomes by their size, 
        and assigns them a pretty name (Chr${num}, numbers sorted in descending ortder)

        Args:
            * file: a path to the chrom.sizes file

        Return:
            dictionary {old_name: new_name}
        """
        chroms_with_sizes: List[Tuple[str, int]] = []
        chrom_counter: int = 1
        for data in read_tab(file):
            if not data[1].isdigit():
                continue
            size: str = int(data[1])
            if size < self.min_chrom_size:
                continue
            chroms_with_sizes.append((data[0], size))
        mapping: Dict[str, str] = {}
        chrom_order: List[str] = []
        for chrom, _ in sorted(chroms_with_sizes, key=lambda x: -x[1]):
            new_name: str = f"Chr{chrom_counter}"
            mapping[chrom] = new_name
            chrom_order.append(new_name)
            chrom_counter += 1
        return mapping, chrom_order
