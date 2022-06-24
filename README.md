# About mst-backbone

mst-backbone is a phylogeny inference program that uses minimum spanning trees to perform constrained search through tree-space. 

Here is the link to the preprint describing the performance of mst-backbone 

https://www.biorxiv.org/content/10.1101/2020.06.30.180315v1

# Instructions for compilation 


```console
git clone https://github.com/prabhavk/mst-backbone-repo.git 
cd mst-backbone-repo/
make

```

# Usage: 

Using a minimum spanning tree framework for constraining search through phylogenetic tree-space

```console
mst-backbone --seq alignment_file.fas --size size_of_subtree --out prefix_for_output_files
```
Input file:
alignment_file.fas: a multiple sequence alignment file in FASTA format
Parameters:
size_of_subtree: size contraint used to select the vertex groups $V_S$ and $V_O$ of the MST, which represent the sequences used in structural EM.
prefix_for_output_files: prefix used for naming output files

The position of the root as inferred under the general Markov model is not always realistic. A more realistic rooting can be obtained by fitting trees under simpler non-reversible models. The model selection procedure described in the paper can be reproduced using the following command. Please note that the model selection procedure is extremely slow!

```console
mst-backbone --seq alignment_file.fas --out prefix_for_output_files --input_tree input_tree_file_name.nwk --perform_model_selection true 
```


Parameters:


size_of_subtree: the size of the vertex sets that induce a subtree and a connected graph in the MST; The sequences of the selected vertices (taxa) are used to infer local phylogenetic trees that are subsequently combined to construct a global phylogenetic tree. Default value of 10 seems adeqate and was tested on multiple empirical and simulated datasets.
