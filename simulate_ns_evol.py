# Conditions under which systematic bias appears in Grove
import ete3
from fileIO import ReadRootedTree, ReadTree
tree_id = "0"
tree_file_name = "/home/kalaghat/exaptation/Projects/MSTBasedForests/data/RAxMLGrove/trees/" + tree_id + "/tree_best.newick"
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




# print(T)

# T = ReadTree(tree_file_name,"newick")
# print(dir(T))
# print(T.vertexCount)

