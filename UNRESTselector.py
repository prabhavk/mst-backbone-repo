from re import T
import sys
import os
import subprocess as sub
import argparse

# parse arguments
parser = argparse.ArgumentParser(description = 'Perform model selection using UNRESTselector')
parser.add_argument('--seq', help = "alignment file name in FASTA format")
parser.add_argument('--out', help = "output prefix used for creating file name", default = 'UNREST_selector')
parser.add_argument('--input_tree', help = "Tree file name in NEWICK format")

args = parser.parse_args()
alignment_file_name = args.seq
# print(alignment_file_name)
# os.path.exists(alignment_file_name)

newick_tree_file_name = args.input_tree
output_prefix = args.out

cwd = os.getcwd()
if (os.path.exists(alignment_file_name)):
    pass
else:
    alignment_file_name = str(cwd) + '/' + str(alignment_file_name)[1:-1]
    # alignment_file_name = os.path.join(cwd,alignment_file_name)
    print("alignment_file_name",alignment_file_name)
    assert(os.path.exists(alignment_file_name))
    

if (os.path.exists(newick_tree_file_name)):
    pass
else:
    newick_tree_file_name = cwd + '/' + str(newick_tree_file_name)
    assert(newick_tree_file_name)
    print("newick_file_name",newick_tree_file_name)        

print("output_prefix",output_prefix)

edgeList_file_name = newick_tree_file_name.replace(".newick",".temp_edges")


from fileIO import ReadRootedTree
RT = ReadRootedTree(newick_tree_file_name,"newick")
RT.WriteToFile(edgeList_file_name,fileFormat="edgeList")
edges = RT.edgeLengths.keys()

modelSelectionFileName = cwd + "/" + output_prefix + "_runSelection.sh"
modelSelectionFile = open(modelSelectionFileName,"w")

# Use mxqsub on server
mxqsub_available = os.path.exists('/project/exaptation/')


if mxqsub_available:    
    for edge in edges:
        modelSelectionFile.write('mxqsub -t 12h -m 2G /project/exaptation/Tools/UNRESTselector_perEdge ' + edgeList_file_name +'\'t' + alignment_file_name + '\t' + edge[0] + '\t' + edge[1] + '\n')
        # replacing mxqsub with qsub should work on most Sun Grid Engine systems
else:
    for edge in edges:
        modelSelectionFile.write('UNRESTselector_perEdge ' + edgeList_file_name +'\t' + alignment_file_name + '\t' + edge[0] + '\t' + edge[1] + '\n')
        modelSelectionFile.write('echo #######################################################################################################\n')


modelSelectionFile.close()    

sub.call('chmod +x ' + modelSelectionFileName,shell=True)
print ('Starting model selection')



# Check to see if all jobs are done

def count_finished_jobs():    
    num_jobs = 0
    for edge in edges:
        BIC_file_name = alignment_file_name+".rootedAt_"+edge[0]+"_"+edge[1]+".log"
        if (os.path.exists(BIC_file_name)):
            BIC_file = open(BIC_file_name,'r')
            for line in BIC_file:
                if line.startswith("BIC"):
                    num_jobs += 1                    
            BIC_file.close()
    return num_jobs


num_jobs = count_finished_jobs()
tot_num_jobs = len(edges)

if num_jobs == tot_num_jobs:
    print ("Model selection complete")
else:
    sub.call(modelSelectionFileName, shell=True)

BIC_list = []
for edge in edges:
    BIC_file_name = alignment_file_name+".rootedAt_"+edge[0]+"_"+edge[1]+".log"
    BIC_file = open(BIC_file_name,"r")
    for line in BIC_file:
        if line.startswith("BIC"):
            BIC = float(line.split("BIC:")[1].strip())
            BIC_list.append((BIC,edge))
    BIC_file.close()

# Select rooted tree/model with minimum BIC
BIC_list.sort(key=lambda x: x[0])
print ("Minimum BIC is ", BIC_list[0])
BIC = BIC_list[0][0]
edge = BIC_list[0][1]
file_prefix = alignment_file_name + ".rootedAt_" + edge[0] + "_" + edge[1]
# Paramters
model_parameters_file_name_best_edge = file_prefix + ".modelParameters"
model_parameters_file_name_to_store = cwd+"/"+output_prefix+".params"
sub.call('cp '+model_parameters_file_name_best_edge+"\t"+model_parameters_file_name_to_store,shell=True)
# log file
log_file_name_best_edge = file_prefix + ".log"
log_file_name_to_store = cwd+"/"+output_prefix+".log"
sub.call('cp '+log_file_name_best_edge+"\t"+log_file_name_to_store,shell=True)
# rooted tree newick
newick_tree_file_name_best_edge = file_prefix + ".newick"
newick_tree_file_name_to_store = cwd+"/"+output_prefix+".newick"
sub.call('cp '+newick_tree_file_name_best_edge+"\t"+newick_tree_file_name_to_store,shell=True)
# rooted tree edges
edgeList_tree_file_name_best_edge = file_prefix + ".edgeList"
edgeList_tree_file_name_to_store = cwd+"/"+output_prefix+".edges"
sub.call('cp '+edgeList_tree_file_name_best_edge+"\t"+edgeList_tree_file_name_to_store,shell=True)


# delete all temp files
for edge in edges:
    file_prefix = alignment_file_name + ".rootedAt_" + edge[0] + "_" + edge[1]
    model_parameters_file_name_best_edge = file_prefix + ".modelParameters"
    newick_tree_file_name_best_edge = file_prefix + ".newick"
    edgeList_tree_file_name_best_edge = file_prefix + ".edgeList"
    log_file_name_best_edge = file_prefix + ".log"
    for file_name in [log_file_name_best_edge,edgeList_tree_file_name_best_edge,model_parameters_file_name_best_edge,newick_tree_file_name_best_edge]:
        if os.path.exists(file_name):
            sub.call('rm '+file_name,shell=True)

if os.path.exists(edgeList_file_name):
    sub.call('rm ' + edgeList_file_name, shell=True)

if (os.path.exists(modelSelectionFileName)):
    sub.call('rm ' + modelSelectionFileName, shell=True)