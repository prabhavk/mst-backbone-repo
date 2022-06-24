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
mst-backbone --seq alignment.fas --size size_of_subtree --out output_prefix
```
Input file:

alignment.fas: a multiple sequence alignment file in FASTA format

Parameters:

size_of_subtree: size constraint used to select the vertex groups $V_S$ and $V_O$ of the MST.



The model selection procedure described in the paper can be reproduced using the following command. Please note that the model selection procedure is extremely slow!

```console
mst-backbone --seq alignment_file.fas --out prefix_for_output_files --input_tree input_tree_file_name.nwk --perform_model_selection true 
```
