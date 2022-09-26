# Conditions under which systematic bias appears in Grove
import ete3
from fileIO import ReadRootedTree, ReadTree
from MarkovModels import GenerateQForStationaryDistribution, Get11FreeRates
from config import projectPath
import os
import sys
import subprocess as sub
from numpy.random import dirichlet
# parameters
# tree_id = sys.argv[1]
# sequence_length = sys.argv[2]
# gc_diff = sys.argv[3]

def Simulate_ns_evol(tree_id,sequence_length,gc_diff):
    tree_file_name = projectPath + "data/RAxMLGrove/trees/" + tree_id + "/tree_best.newick"
    # T = ete3.Tree(tree_file_name)
    print(tree_file_name)
    T = ReadRootedTree(tree_file_name,"newick")
    # T.ladderize(direction=0)
    # Select nodes that have four grandchildren
    preOrderTraversal = T.GetPreOrderTraversalWithoutLeaves()
    leafNames = T.GetLeafNames()
    branch_length_ratio = []
    for node in preOrderTraversal:        
        # print (node.children)
        if len(node.children) == 2:        
            child_l, child_r = node.children
            max_child_length = max(T.GetEdgeLength(node.name,child_l.name),T.GetEdgeLength(node.name,child_r.name))
            if child_l.name not in leafNames and child_r.name not in leafNames:
                grandchild_l_l, grandchild_l_r = child_l.children
                grandchild_r_l, grandchild_r_r = child_r.children
                min_grandchild_left_length = max(T.GetEdgeLength(child_l.name,grandchild_l_l.name),T.GetEdgeLength(child_l.name,grandchild_l_r.name))
                min_grandchild_right_length = max(T.GetEdgeLength(child_r.name,grandchild_r_l.name),T.GetEdgeLength(child_r.name,grandchild_r_r.name))        
                branch_length_ratio.append((node,max_child_length/min(min_grandchild_left_length,min_grandchild_right_length)))
                # break  
                
    # T.show()        

    branch_length_ratio.sort(key=lambda x: x[1])
    print(branch_length_ratio[0])
    node = branch_length_ratio[0][0]
    child_l, child_r = node.children
    grandchild_l_l, grandchild_l_r = child_l.children
    grandchild_r_l, grandchild_r_r = child_r.children
    change_point_node_names = [grandchild_l_r.name, grandchild_r_l.name]

    at = 0.5 - gc_diff
    # gc = 1 - at    
    pi_at_scaled = [at/2, (1-at)/2, at/2, (1-at)/2]
    pi_gc_scaled = [(1-at)/2, at/2, (1-at)/2, at/2]
    for i in range(4):
        pi_at_scaled[i] *= 50
        pi_gc_scaled[i] *= 50
    pi_at = dirichlet(pi_at_scaled)    
    pi_gc = dirichlet(pi_gc_scaled)
    print(pi_at)    
    print(pi_gc)

    

    Q_at = GenerateQForStationaryDistribution(pi_at)
    Q_gc = GenerateQForStationaryDistribution(pi_gc)

    rateVector_at = Get11FreeRates(Q_at)
    rateVector_gc = Get11FreeRates(Q_gc)
    T.rateVectorForCat[1] = rateVector_at
    T.rateVectorForCat[2] = rateVector_gc

    T.rateCategoryForVertex[change_point_node_names[0]] = 2
    T.rateCategoryForVertex[change_point_node_names[1]] = 2
    root_name = T.root.name
    T.rateCategoryForVertex[root_name] = 1
    excluded_names = change_point_node_names + [root_name]

    for node_name in T.vertices:
        node = T.GetVertex(node_name)
        if node_name not in excluded_names:        
            T.rateCategoryForVertex[node_name] = T.rateCategoryForVertex[node.parent.name]

    T.sequenceLength = sequence_length
    # p_change = 0.5
    # RT.GenerateMarkovModulatedMarkovModel(p_change)

    if not os.path.exists(projectPath + 'data/grove_sequences/'+tree_id):
        sub.call('mkdir '+ projectPath + 'data/grove_sequences/'+tree_id, shell=True)

    sequence_file_suffix = 'sequences_tree_id_'+tree_id+'_GC_diff_'+str(gc_diff)+'_sequence_length_' + str(T.sequenceLength)
    # sequence_file_name_sim = projectPath + 'data/grove_sequences/' + tree_id + '/' + sequence_file_suffix + '.fas'
    T.sequence_file_suffix = sequence_file_suffix

    T.control_file_name = projectPath + 'data/grove_sequences/'+tree_id+'/control.txt'
    T.WriteControlFile()
    T.SimulateEvolutionUsingIndelible()



selected_tree_ids_file = open(projectPath + 'data/selected_grove_tree_ids','r')
tree_ids_num_taxa = []
for line in selected_tree_ids_file:
    tree_id, num_taxa = line.split(",")
    tree_id = tree_id.strip()
    num_taxa = num_taxa.strip()
    tree_ids_num_taxa.append((tree_id,num_taxa))

selected_tree_ids_file.close()
# print (tree_ids[0:5])

tree_ids_num_taxa.sort(key = lambda x: x[1])
tree_ids_selected = [tree_ids_num_taxa[i] for i in range(1,19400,100)]
print (len(tree_ids_selected))


print(tree_ids_selected[0][0])

sequence_length_range = [1000, 2000, 4000, 8000, 16000]
gc_diff_range = [0, 0.1, 0.2, 0.3]
Simulate_ns_evol(tree_ids_selected[0][0], sequence_length = 1000, gc_diff = 0.3)
# run tools
# quantify accuracy


# print(T)

# T = ReadTree(tree_file_name,"newick")
# print(dir(T))
# print(T.vertexCount)

