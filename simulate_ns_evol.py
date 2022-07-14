# Conditions under which systematic bias appears in Grove
import ete3
from fileIO import ReadRootedTree, ReadTree
from MarkovModels import GenerateQForStationaryDistribution, Get11FreeRates
from config import projectPath


tree_id = "0"
tree_file_name = projectPath + "MSTBasedForests/data/RAxMLGrove/trees/" + tree_id + "/tree_best.newick"
# T = ete3.Tree(tree_file_name)
T = ReadRootedTree(tree_file_name,"newick")
# T.ladderize(direction=0)
# Select nodes that have four grandchildren
preOrderTraversal = T.GetPreOrderTraversalWithoutLeaves()
leafNames = T.GetLeafNames()
branch_length_ratio = []
for node in preOrderTraversal:        
    # print (node.children)
    # if not node.is_leaf():
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
branch_length_ratio[0]

at = 0.3
# gc = 1 - at
pi_at = [at/2, (1-at)/2, at/2, (1-at)/2]
print(pi_at)
pi_gc = [(1-at)/2, at/2, (1-at)/2, at/2]
print(pi_gc)

Q_at = GenerateQForStationaryDistribution(pi_at)
Q_gc = GenerateQForStationaryDistribution(pi_gc)

rateVector_at = Get11FreeRates(Q_at)
rateVector_gc = Get11FreeRates(Q_gc)
T.rateVectorForCat[1] = rateVector_at
T.rateVectorForCat[2] = rateVector_gc

for vertex in T.vertices:
    T.rateCategoryForVertex[vertex] = 1

T.rateCategoryForVertex['tax2'] = 2
T.rateCategoryForVertex['tax3'] = 2

sequenceLength = 10000 # 
# p_change = 0.5
# RT.GenerateMarkovModulatedMarkovModel(p_change)
sequence_file_suffix = 'sequences_tree_id_42_num_taxa_4_sequence_length_' + str(sequenceLength)
sequence_file_name_sim = projectPath + 'MSTBasedForests/data/quartet_ns/' + sequence_file_suffix + '.fas'
T.sequence_file_suffix = sequence_file_suffix


T.control_file_name = projectPath + 'MSTBasedForests/data/quartet_ns/control.txt'
T.WriteControlFile()
T.SimulateEvolutionUsingIndelible()




# print(T)

# T = ReadTree(tree_file_name,"newick")
# print(dir(T))
# print(T.vertexCount)

