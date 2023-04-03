import os
import subprocess as sub
cwd = os.getcwd()
mst_backbone_path =  cwd + "/mst-backbone"
alignment_path = cwd + "/empirical_results/HIV_1_env/BelgianHIVTransmissionChain_noGaps_noAmbiguousChars.fas"
mst_script = mst_backbone_path + " --seq " + alignment_path +  " --out HIV_trans"
#1 build phylogeny using mst-backbone
sub.call(mst_script,shell=True)
#2 visualize MST
#3 visualize phylogeny

