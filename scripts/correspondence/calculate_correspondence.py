def ReadMST(mst_file_name):
    mst_file = open(mst_file_name,"r")
    for line in mst_file:
        

    mst_file.close()



def Compute_correspondence(mst_file_name,phylogeny_file_name):
    P = ReadPhylogeneticTree(phylogeny_file_name)
    M = ReadMST(mst_file_name)
