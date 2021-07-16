import sys
import subprocess as sub

alignment_file_name = sys.argv[1]
mst_backbone_file_name = "path_to_mst-backbone"

mst_backbone_script = ""
for (subtree_size in [10, 20, 40, 50, 100]):
	output_file_name_prefix = alignment_file_name + "subtree_size_" + subtree_size
	mst_backbone_script += mst_backbone_file_name + "\t--seq\t" + alignment_file_name + "\t--sub\t" + subtree_size + "\t--out\t" + output_file_name_prefix + "\n"

sub.call(mst_backbone_script)



