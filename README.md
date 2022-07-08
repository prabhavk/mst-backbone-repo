# About mst-backbone

mst-backbone is a phylogeny inference program that uses minimum spanning trees to perform constrained search through tree-space. 

Here is the link to the preprint describing the performance of mst-backbone 

https://www.biorxiv.org/content/10.1101/2020.06.30.180315v1

# Installation


```console
git clone https://github.com/prabhavk/mst-backbone-repo.git 
cd mst-backbone-repo/
make all

```

mst-backbone constructs phylogenetic trees by constraining the search through tree-space. Tree search is performed under the general Markov model using structural EM. 

 
UNRESTselector.py performs model selection by fitting non-stationary non-reversible non-homogeneous CT-HMM to phylogenetic trees. 

# Using mst-backbone

```console
mst-backbone --seq alignment.fas --size subtree_size --out output_prefix
```
Input file:

alignment.fas: a multiple sequence alignment file in FASTA format

Parameters:

subtree_size: constraint on the minimum number of non-zero weighted edges in the subgraphs induced by vertex sets $V_S$ and $V_O$ 

output_prefix: used for naming output files. Default prefix is mstbackbone_output

Output files:

output_prefix.mst: the maximum parsimony spanning tree used for constraining search through tree-space 

output_prefix.newick: the output tree in NEWICK format

output_prefix.edges: the output tree in EDGE LIST format

output_prefix.mstbackbone_log: contains all the messages that is printed to terminal

# Using UNRESTselector.py

The model selection procedure described in the paper can be reproduced as follows.

```console
python3 UNRESTselector.py --seq alignment.fas --out output_prefix --input_tree tree.newick
```

Input files:

alignment.fas: a multiple sequence alignment file in FASTA format

tree.newick: a phylogenetic tree in NEWICK format. The input tree is re-rooted such that BIC is minimized.

Parameters:

output_prefix: prefix used for naming output files. Default prefix is alignment_UNRESTselector

Output files:

output_prefix.newick: the output tree in newick format. 

output_prefix.edges: the output tree in edge list format. 

output_prefix.params: the parameters of the selected model

output_prefix.log: log file containing BIC and time elapsed for performing model selection