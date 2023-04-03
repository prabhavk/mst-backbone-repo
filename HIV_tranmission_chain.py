import os
import subprocess as sub
import networkx as n
import matplotlib.pyplot as plt
cwd = os.getcwd()
mst_backbone_path =  cwd + "/mst-backbone"
alignment_path = cwd + "/empirical_results/HIV_1_env/BelgianHIVTransmissionChain_noGaps_noAmbiguousChars.fas"
phylogeny_file = cwd + "/empirical_results/HIV_1_env/HIV_trans.unrooted_newick"
mst_file = cwd + "/empirical_results/HIV_1_env/HIV_trans.initial_MST"
if not (os.path.exists(phylogeny_file)):
    mst_script = mst_backbone_path + " --seq " + alignment_path +  " --out HIV_trans"
#1 build phylogeny using mst-backbone
    sub.call(mst_script,shell=True)
#2 visualize MST
M = n.read_weighted_edgelist(mst_file)
plt.figure()
n.draw_kamada_kawai(M)
plt.show()
# 
#3 visualize phylogeny

