
tree_names_file = open("/home/kalaghat/exaptation/Projects/MSTBasedForests/data/selected_grove_tree_ids","r")
tree_names = []
tree_names = tree_names_file.readline().split(",")
print (len(tree_names))
for tree_name in tree_names[:10]:
    print (tree_name)  
  
# Select lineages for evolving 

# Run SEM-M and RAxML-NG
