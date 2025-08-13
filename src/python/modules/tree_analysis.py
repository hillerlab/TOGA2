from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from collections import Counter
# from io import StringIO
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import os
import sys

__author__ = 'Amy Stephen'
__credits__ = ('Yury V. Malovichko')
__year__ = '2024'

ENSEMBL_PATH = '/beegfs/home/astephen/resolveTOGA/files/human_mouse_homologs.tsv'
# ENSEMBL_PATH = '/beegfs/home/astephen/resolveTOGA/files/human_bosTau9_homologs.tsv'
# ENSEMBL_PATH = '/beegfs/home/astephen/resolveTOGA/files/human_panTro6_homologs.tsv'

"""
D = 0
P = 1
S = 2
Dq = 3
Dr = 4
R = 5
Q = 6
"""
CAT_DICT_NUM = {0 : 'D', 1 : 'P', 2: 'S', 3 : 'Dq', 4 : 'Dr', 5 : 'R', 6 : 'Q'}
CAT_DICT = {'D' : 0, 'P' : 1, 'S' : 2, 'Dq' : 3, 'Dr' : 4, 'R' : 5, 'Q' : 6}
CAT_DDICT = [
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 2, 1, 1, 1, 1], ## YM: So ancestor for S-S is also S?
            [0, 0, 1, 3, 0, 0, 3],
            [0, 0, 1, 0, 4, 4, 0],
            [0, 0, 1, 0, 4, 4, 2],
            [0, 0, 1, 3, 0, 2, 3]]

def label_clade(clade):
    """labels the nodes with D, P, S, Dq, Dr, R, Q"""

    if not clade.is_terminal():
        # gets category of children nodes
        left = clade.clades[0]
        right = clade.clades[1]
        l = left.name
        r = right.name
        if left.is_terminal():
            l = 'Q'
            if left.name[:3] == '#R#' or left.name[:3] == '_R_':
                l = 'R'
        if right.is_terminal():
            r = 'Q'
            if right.name[:3] == '#R#' or right.name[:3] == '_R_':
                r = 'R'
        # renames current node based on category of children
        cat_num = CAT_DDICT[CAT_DICT[l]][CAT_DICT[r]]
        clade.name = CAT_DICT_NUM[cat_num]

def make_cat_tree(tree):
    # traverse tree by reverse level-order
    for clade in list(tree.find_clades(order="level"))[::-1]:
        label_clade(clade)


def special_case(tree, bootthresh: float, iqtree: bool = False):
    """
    special case: when we could resolve into 2 one2ones but the bootstrap score is low - have to do an extra check
    some bootstrap scores are implicit so this checks in 2 one2ones case bootstrap > threshold
    """
    count = 0
    remove = []
    attr: str = 'confidence' if iqtree else 'comment'
    # count number of S
    for clade in list(tree.find_clades()):
        print(f'{clade=}, {clade.root=}')
        print(f'{clade.name=}, {clade.is_terminal()=}, {clade.comment=}, {clade.confidence=}, {clade.__getattribute__(attr)=}, {iqtree=}, {bootthresh=}')
        if clade.name == 'S':
            count = count + 1
    # only checks for specific case
    if count == 3:
        # if any(clade.comment is not None and int(clade.comment)<int(bootthresh) for clade in tree.find_clades()):
        if any(
            clade.__getattribute__(attr) is not None and float(clade.__getattribute__(attr)) < float(bootthresh) for clade in tree.find_clades()
        ):
            for leaf in list(tree.get_terminals()):
                remove.append(leaf.name)
    return remove

def can_resolve(
    tree: Phylo.BaseTree.Tree, bootthresh: float, iqtree: bool = False
) -> List[Tuple[str, str]]:
    """
    returns list of tuples of isoforms that can be resolved as one2one
    algo: For each one2one candidate leaf couple (parent is S ie. speciation node)
            If P is not in the path from leaf to root AND
            If all the bootstrap scores on path from leaf to root > boot threshold
            Then the couple can be resolved
    """
    resols = []
    attr: str = 'confidence' if iqtree else 'comment'
    for leaf in list(tree.get_terminals()):
        node_path = tree.get_path(leaf)
        parent = node_path[-1]
        if len(node_path) > 1:
            parent = node_path[-2]
        # print(f'{node_path=}')
        # for i, clade in enumerate(node_path):
        #     print(f'{i=}, {clade=}, {clade.name=}, {clade.comment=}')
        # of the leaves that are a result of speciation: do...
        if str(parent) == 'S':
            # check if 'P'roblem node on path back to root
            resol_P = not (any(clade.name == 'P' for clade in node_path))
            # check bootstrap values for path back to root except for root's bootstrap (hence the path[1:])
            resol_bootstrap = not any(
                clade.__getattribute__(attr) is not None and float(clade.__getattribute__(attr)) < float(bootthresh) for clade in node_path
            )
            # if not: append leaf name to resols
            if (resol_P and resol_bootstrap):
                resols.append(leaf.name)

    # special case where could be 2 one2ones but check if bootstrap at root high enough
    if len(tree.get_terminals()) == 4:
        remove = special_case(tree, bootthresh, iqtree)
        # print(f'{remove=}')
        resols = [item for item in resols if item not in remove]

    # more digestible format
    res = []
    ## YM: And a final sanity check: Remove the resolved pairs involving duplicated projections
    ## TODO: Discuss with Michael
    case_num: Dict[str, int] = Counter([get_tr(x) for i, x in enumerate(resols) if x % 2])
    for i in range(0, len(resols), 2):
        query_tr: str = get_tr(resols[i])
        if case_num[query_tr] > 1:
            continue
        res.append((resols[i], resols[i+1]))
    return res



def get_gene(name):
    return name.split('#')[-2]


def get_tr(proj: str) -> str:
    return '#'.join(proj.split('#')[:-1])


def check_ensembl(one_ones, ensembl_path, out_path):
    if one_ones == []:
        return
    with open(out_path, 'w') as out:
        for one2one in one_ones:
            with open(ensembl_path, 'r') as file:
                gene1 = get_gene(one2one[0])
                gene2 = get_gene(one2one[1])
                string = '>"' + str(gene1) + ', ' + str(gene2) + '\n'
                out.write(string)
                for line in file:
                    entries = line.split()
                    if (gene1 in entries) or (gene2 in entries):
                        out.write(line.strip())
                        out.write('\n')

def do_dir(d, bootthresh):
    for file_name in os.listdir(d):
        if 'nwk' in file_name:
            file_path = os.path.join(d, file_name)
            tree = Phylo.read(file_path, "newick")
            # make the categorized tree
            make_cat_tree(tree)
            msa_x = file_name.split('.')[-2]
            out_svg = msa_x + '_cat.svg'

            # make the visual
            visualize_newick(tree, out_svg)

            res = can_resolve(tree, bootthresh)
            out_name = msa_x + '_' + bootthresh + '_ensemblout.txt'
            # compare against ensembl
            check_ensembl(res, ENSEMBL_PATH, out_name)

def process(
    file_path: str,
    bootthresh: int,
    out_svg: Optional[Union[str, None]] = None,
    out_ens: Optional[Union[str, None]] = None,
    ensembl_path: Optional[Union[str, None]] = None
):
    tree = Phylo.read(file_path, "newick")
    # make the categorized tree
    make_cat_tree(tree)

    # visualize_newick(tree, out_svg)

    res = can_resolve(tree, bootthresh)
    return res
    # compare against ensembl
    # check_ensembl(res, ensembl_path, out_ens)

def main():
    # quick_test()

    d = sys.argv[1]
    bootthresh = sys.argv[2]
    do_dir(d, bootthresh)

# def quick_test():
#     make_cat_tree(tree)
#     # Phylo.write(tree,"tmp.nwk", format='newick')
#     can_resolve(tree)

def visualize_newick(newick_file, out_file):
    tree = newick_file
    # Start with a reasonably large size
    fig_width, fig_height = 14, 14
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)

    # Check if the content is fully visible
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    if not ax.get_autoscale_on() and bbox.width < fig_width and bbox.height < fig_height:
        # Content is fully visible, save the figure and return
        plt.savefig(out_file, format='svg')
        plt.close(fig)
        return

    # If content is not fully visible, increase the figure size
    while True:
        fig_width *= 1.1
        fig_height *= 1.1
        fig.set_size_inches(fig_width, fig_height, forward=True)
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        if not ax.get_autoscale_on() and bbox.width >= fig_width and bbox.height >= fig_height:
            break

    # Save the figure
    plt.savefig(out_file, format='svg')
    plt.close(fig)


def midpoint_root(unrooted_path, rooted_path):
    tree = Phylo.read(unrooted_path, "newick")
    tree.root_at_midpoint()
    Phylo.write(tree, rooted_path, "newick")

# if __name__ == '__main__':
#     main()
