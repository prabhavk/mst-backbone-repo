# mst-backbone-sem-gm

mst-backbone(sem-gm) comprises two algorithms:

1) mst-backbone: A minimum spanning tree framework for constraining search through phylogenetic tree-space.

2) sem-gm: A structural EM method for performing tree search under the general Markov model (GM).

Here is the link to the preprint describing the performance of mst-backbone(sem-gm)

https://www.biorxiv.org/content/10.1101/2020.06.30.180315v1

mst-backbone(sem-gm) outputs phylogenetic trees that are rooted under the GM model. 
Results on virus datasets (HIV and H3N2) suggests that it is possible to infer a realistic unrooted topology using mst-backbone(sem-gm) but the GM model may be overparameterized for inferring the location of the root.

Usage: 

mst-backbone-sem-gm   multiple_sequence_alignment_file  tree_depth

Paramters:

A multiple sequence alignment file either in FASTA format or PHYLIP format

tree_depth: the size of the vertex sets that induce a subtree and a connected graph in the MST; The sequences of the selected vertices (taxa) are used to infer local phylogenetic trees that are subsequently combined to construct a global phylogenetic tree. Default value of 10 seems adeqate and was tested on multiple empirical and simulated datasets.
