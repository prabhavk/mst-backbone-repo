import re
import subprocess as sub
import os
from classDeclarationsAndFunctions import Tree, RootedTree
from numpy import array
from MarkovModels import DNA
from config import projectPath, scriptPath
import numpy as np
import ete3
################################### ALIGNMENT I/O #####################################

def ReadFasta(fileName):    
    fastaFile = open(fileName,'r')
    seq=''
    name=''
    sequenceAlignment={}
    for line in fastaFile:
        if line.startswith('>'):
            if seq != '':
                seq = seq.upper()
                sequenceAlignment[name] = seq
                seq = ''
            name = line.strip().split('>')[1]
        else:
            seq += line.strip()

    sequenceAlignment[name] = seq
    fastaFile.close()
    return (sequenceAlignment)

def ReadAncestralSequences(fileName):
    phylipFile = open(fileName,'r')
    sequenceAlignment={}
    for line in phylipFile:
        splitLine = line.strip().split(' ')
        if len(splitLine)>1:
            name = splitLine[0]
            seq = splitLine[len(splitLine)-1]
        else:
            name, seq = line.strip().split('\t')
        sequenceAlignment[name] = seq
    phylipFile.close()
    return (sequenceAlignment)

def ReadPhylip(fileName):
    phylipFile = open(fileName,'r')
    phylipFile.readline()
    sequenceAlignment={}
    for line in phylipFile:
        splitLine = line.strip().split(' ')
        if len(splitLine)>1:
            name = splitLine[0]
            seq = splitLine[len(splitLine)-1]
        else:
            name, seq = line.strip().split('\t')
        sequenceAlignment[name] = seq
    phylipFile.close()
    return (sequenceAlignment)

def ReadAlignment(fileName):
    alignmentFile = open(fileName,'r')
    firstLine = alignmentFile.readline()
    alignmentFile.close()
    if firstLine.startswith('>'):
        alignment = ReadFasta(fileName)
    else:
        alignment = ReadPhylip(fileName)
    return alignment
    

def WriteAlignment(alignment,fileName,fileFormat="fasta",convertToAlpha=False):
    alignmentFile = open(fileName,'w')
    if fileFormat == "fasta":
        nucList = ['A','C','T','G']
        if convertToAlpha:
            for seqId in sorted(alignment.keys()):
                alignmentFile.write('>'+str(seqId)+'\n')
                for char in alignment[seqId]:
                    alignmentFile.write(nucList[int(char)])
                alignmentFile.write('\n')
        else:
            for seqId in sorted(alignment.keys()):
                alignmentFile.write('>'+str(seqId)+'\n')
                alignmentFile.write(str(alignment[seqId])+'\n')
    elif fileFormat=="phylip":
        numberOfSequences = len(alignment.keys())
        sequenceLength = len(alignment.values()[0])
        alignmentFile.write(str(numberOfSequences)+'\t'+str(sequenceLength)+'\n')
        if convertToAlpha:
            nucList = ['A','C','T','G']
            for seqId in sorted(alignment.keys()):
                alignmentFile.write(str(seqId)+'\t')
                for char in alignment[seqId]:
                    alignmentFile.write(nucList[int(char)])
                alignmentFile.write('\n')
        else:            
            for sequenceName in alignment.keys():
                alignmentFile.write(sequenceName+'\t'+alignment[sequenceName]+'\n')
    alignmentFile.close()

def WriteBootstrapAlignmentForReplicate(sequenceFileName,bootstrapReplicate):
    alignment = ReadAlignment(sequenceFileName)
    sequenceLength = len(alignment.values()[0])
    posList = range(0,sequenceLength)
    bootStrapPos = np.random.choice(posList,size=sequenceLength,replace=True)
    bootstrapAlignment = {}
    for seqName, seq in alignment.items():
        bootstrapSeq = ''
        for pos in bootStrapPos:
            bootstrapSeq += seq[pos]
        bootstrapAlignment[seqName] = bootstrapSeq
    bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
    WriteAlignment(bootstrapAlignment,bootstrapAlignmentFileName)

def WriteBatchScriptForConstructingBootstrapAlignments(sequenceFileName,expName,numberOfBootstrapReplicates = 100):
    batchCommandFile = open(projectPath+"scripts/batchCommandForConstructingBootstrapReplicates_"+expName+".sh","w")
    for bootstrapReplicate in range(1,(numberOfBootstrapReplicates+1)):
        batchCommandFile.write("python2.7\t" + scriptPath + "batch_constructBootstrapAlignment.py\t"+sequenceFileName + "\t" + str(bootstrapReplicate)+"\n")
    batchCommandFile.close()

def WriteBatchScriptForMSTBackboneSEMGMForBootstrapAlignments(sequenceFileName,expName,numberOfBootstrapReplicates = 100):
    batchCommandFile = open(projectPath+"scripts/batchCommandForMSTBackboneSEMGMForBootstrapReplicates_"+expName+".sh","w")
    treeDepth = 10
    for bootstrapReplicate in range(1,(numberOfBootstrapReplicates+1)):
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        batchCommandFile.write(scriptPath+"mst-backbone-SEM\t"+bootstrapAlignmentFileName+"\t"+str(treeDepth)+"\n")
    batchCommandFile.close()
    
def WriteBatchScriptForSEMGMForBootstrapAlignments(sequenceFileName,expName,numberOfBootstrapReplicates = 100):
    batchCommandFile = open(projectPath+"scripts/batchCommandForSEMGMForBootstrapReplicates_"+expName+".sh","w")
    alignment = ReadAlignment(sequenceFileName)
    treeDepth = len(alignment.keys())
    for bootstrapReplicate in range(1,(numberOfBootstrapReplicates+1)):
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        batchCommandFile.write(scriptPath+"mst-backbone-SEM\t"+bootstrapAlignmentFileName+"\t"+str(treeDepth)+"\n")
    batchCommandFile.close()

    
def WriteBatchScriptForModelSelection(sequenceFileName,expName,type="w",groupId=""):
    edgeListFileName = sequenceFileName + ".edgeList"
    edgeListFile = open(edgeListFileName, "r")
    edges = []
    for line in edgeListFile:
        (u_name, v_name, t) = line.strip().split("\t")
        edges.append((u_name, v_name))
    edgeListFile.close()
    mxqsubPrefix=""    
    if 'exaptation' in os.getcwd():
        mxqsubPrefix = "mxqsub -t 12h -m 2G "
        if groupId != "":
            mxqsubPrefix += str(groupId)+" " 
            
    batchCommandFile = open(projectPath+"scripts/batchCommandForModelSelection_"+expName+".sh",type)
    for edge in edges:
        (u_name, v_name) = sorted(edge)
        batchCommandFile.write(mxqsubPrefix+scriptPath+"markovModelSelectorForRootedTree\t"+edgeListFileName+"\t"+sequenceFileName+"\t"+u_name+"\t"+v_name+"\n")        
    batchCommandFile.close()
    
def WriteBatchScriptForModelSelectionForBootstrapAlignments(sequenceFileName, expName, bootstrapReplicates_end, bootstrapReplicates_start = 1):    
    expNameBoot = expName+"_boot"    
    for bootstrapReplicate in range(bootstrapReplicates_start,(bootstrapReplicates_end+1)):
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'        
        if (bootstrapReplicate == bootstrapReplicates_start) :
            WriteBatchScriptForModelSelection(bootstrapAlignmentFileName,expNameBoot,"w")
        else:
            WriteBatchScriptForModelSelection(bootstrapAlignmentFileName,expNameBoot,"a")

def WriteBatchScriptForModelSelectionForIncompleteJobsForBootstrapAlignments(sequenceFileName, expName, numberOfBootstrapReplicates):
    expNameBoot = expName+"_inc_boot"
    for bootstrapReplicate in range(1,(numberOfBootstrapReplicates+1)):
        print ("bootstrapReplicate: ", bootstrapReplicate) 
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'        
        if (bootstrapReplicate == 1) :
            WriteBatchScriptForModelSelectionForIncompleteJobs(bootstrapAlignmentFileName,expNameBoot,"w")
        else:
            WriteBatchScriptForModelSelectionForIncompleteJobs(bootstrapAlignmentFileName,expNameBoot,"a")

    
def  WriteBatchScriptForModelSelectionForIncompleteJobs(sequenceFileName, expName,type="w",groupId=""):
    edgeListFileName = sequenceFileName + ".edgeList"
    edgeListFile = open(edgeListFileName,"r")
    edges = []
    for line in edgeListFile:
        (u_name, v_name, t) = line.strip().split("\t")
        edges.append((u_name, v_name))
    edgeListFile.close()
        
    mxqsubPrefix=""    
    if 'exaptation' in os.getcwd():
        mxqsubPrefix = "mxqsub -t 12h -m 2G "
        if groupId != "":
            mxqsubPrefix += str(groupId)+" " 
    batchCommandFile = open(projectPath+"scripts/batchCommandForModelSelectionTBeCompleted_"+expName+".sh",type)
    for edge in edges:
        (u_name, v_name) = sorted(edge)
        filePrefix_uv = sequenceFileName + ".rootedAt_"+u_name+"_"+v_name
        filePrefix_vu = sequenceFileName + ".rootedAt_"+v_name+"_"+u_name
        if not os.path.isfile(filePrefix_uv+".log") and not os.path.isfile(filePrefix_vu+".log"):
            print ("log file not found for ", filePrefix_uv)  
            batchCommandFile.write(mxqsubPrefix+scriptPath+"markovModelSelectorForRootedTree\t"+edgeListFileName+"\t"+sequenceFileName+"\t"+u_name+"\t"+v_name+"\n")
            print ("---------------")
        elif not os.path.isfile(filePrefix_uv+".edgeList") and not os.path.isfile(filePrefix_vu+".edgeList"):
            print ("edge list not found for ", filePrefix_uv)
            batchCommandFile.write(mxqsubPrefix+scriptPath+"markovModelSelectorForRootedTree\t"+edgeListFileName+"\t"+sequenceFileName+"\t"+u_name+"\t"+v_name+"\n")
            print ("---------------")
        else:
            if os.path.isfile(filePrefix_uv+".log"):
                logFile = open(filePrefix_uv+".log","r")
            else:
                logFile = open(filePrefix_vu+".log","r")
            BICNotFound = True
            for line in logFile: 
                if line.startswith("BIC"):
                    BICNotFound = False
            logFile.close()
            if BICNotFound:
                batchCommandFile.write(mxqsubPrefix+scriptPath+"markovModelSelectorForRootedTree\t"+edgeListFileName+"\t"+sequenceFileName+"\t"+u_name+"\t"+v_name+"\n")            
    batchCommandFile.close()
    
#######################################  DISTANCES I/O #########################################

def ReadDistances(fileName,vertexNameList='',distanceEntryType='line'):
    if distanceEntryType=='matrix':
        distances={}
        # vertexNameList must not be empty
        distFile = open(fileName,'r')
        vertexID_1=0
        for line in distFile:
            vertexName_1 = vertexNameList[vertexID_1]
            vertexID_2=0
            for distance in line.strip().split('\t'):
                vertexName_2 = vertexNameList[vertexID_2]
                if vertexName_1 < vertexName_2:
                    distances[(vertexName_1,vertexName_2)] = float(distance)
                vertexID_2+=1
            vertexID_1+=1
    else:
        distFile = open(fileName,'r')
        lineSplit = distFile.readline().strip().split('\t')
        distFile.close()
        distFile = open(fileName,'r')
        if len(lineSplit)==3:
        # tab seperated file
            distances = {}
            for line in distFile:
                id1, id2, value = line.strip().split('\t')
                distances[tuple(sorted([id1,id2]))] = float(value)
        else:
        # produced by raxml. Each line is "<id1> <id2> \t <distance>"
            distances = {}
            for line in distFile:
                id1, id2,_,value = line.strip().split(' ')
                distances[tuple(sorted([id1,id2]))] = float(value)
    return (distances)

def WriteDistances(distances,distancesFileName,orderedVertices=''):
    if orderedVertices=='':
        vertexNameList = [x[0] for x in distances.keys()]
        vertexNameList+= [x[1] for x in distances.keys()]
        vertexNameList = list(set(vertexNameList))
        vertexNameList.sort()
    else: vertexNameList = orderedVertices
    distanceFile = open(distancesFileName,'w')
    for i in range(len(vertexNameList)):
        for j in range(i+1,len(vertexNameList)):
            distanceFile.write(vertexNameList[i]+'\t'+vertexNameList[j]+'\t'+str(distances[tuple(sorted([vertexNameList[i],vertexNameList[j]]))])+'\n')
    distanceFile.close()

#######################################  TREE I/O #########################################

def DoesLogFileContainBICInfo(logFileName):
    logFile = open(logFileName,"r")
    valToReturn = False
    for line in logFile:
        if line.startswith("BIC"):
            valToReturn = True
            break
    logFile.close()
    return (valToReturn)

def ComputeElapsedCPUTimeForModelSelection(sequenceFileName):
    edgeListFileName = sequenceFileName + ".edgeList"
    edgeListFile = open(edgeListFileName,"r")
    edges = []
    for line in edgeListFile:
        (u_name, v_name, t) = line.strip().split("\t")
        edges.append((u_name, v_name))
    edgeListFile.close()
    nEdges = len(edges)/1.0
    numberOfTimesFound= 0.0
    cumulative_elapsed_time = 0.0
    for edge in edges:
        (u_name, v_name) = edge        
        prefix_uv = sequenceFileName + ".rootedAt_"+u_name+"_"+v_name
        prefix_vu = sequenceFileName + ".rootedAt_"+v_name+"_"+u_name
        if os.path.isfile(prefix_uv+".log"):
                LogFileName = prefix_uv+".log"
        if os.path.isfile(prefix_vu+".log"):
                LogFileName = prefix_vu+".log"
        LogFile = open(LogFileName,"r")
        TimeLine = ""
        for line in LogFile:
            if line.startswith("Total CPU"):
                TimeLine = line                 
        if TimeLine == "" or len(TimeLine.split("is")) != 2:            
            print (TimeLine)
            print (LogFileName)
            print (sequenceFileName)
            print (edge)
        else:
            elapsed_time = float(TimeLine.split("is")[1].split("s")[0].strip())
            numberOfTimesFound += 1.0
            cumulative_elapsed_time += elapsed_time
    cumulative_elapsed_time*= nEdges/numberOfTimesFound# average over files not found
    return cumulative_elapsed_time
    
def StoreTreeSelectedViaModelSelection(sequenceFileName):
    edgeListFileName = sequenceFileName + ".edgeList"
    edgeListFile = open(edgeListFileName,"r")
    edges = []
    selectedPrefix = ""
    prefix_uv = ""
    prefix_vu = ""
    prefix = ""
    logFileContainsBICInfo = False
    for line in edgeListFile:
        (u_name, v_name, t) = line.strip().split("\t")
        edges.append((u_name, v_name))
    edgeListFile.close()
    minBIC = 0
    numberOfFilesOpened = 0
    #optimalEdge = ("","")    
    for edge in edges:
        BICFileName = ""
        BIC = pow(10,10)
        numberOfFilesOpened += 1
        (u_name, v_name) = edge        
        prefix_uv = sequenceFileName + ".rootedAt_"+u_name+"_"+v_name
        prefix_vu = sequenceFileName + ".rootedAt_"+v_name+"_"+u_name
        if os.path.isfile(prefix_uv+".log"):
            logFileContainsBICInfo = DoesLogFileContainBICInfo(prefix_uv+".log")
            if (logFileContainsBICInfo):            
                BICFileName = prefix_uv+".log"
        if os.path.isfile(prefix_vu+".log"):
            logFileContainsBICInfo = DoesLogFileContainBICInfo(prefix_vu+".log")
            if (logFileContainsBICInfo):
                BICFileName = prefix_vu+".log"
        BICFile = open(BICFileName,"r")
        BICLine = ""
        for line in BICFile:
            if line.startswith("BIC"):
                BICLine = line         
#         print len(BICLine.strip().split("BIC:"))
        if BICLine == "" or len(BICLine.strip().split("BIC:")) != 2:
            print (BICLine)
            print (BICFileName)
            print (sequenceFileName)
            print (edge)
        else:
            BIC = float(BICLine.strip().split("BIC:")[1].strip())
        if (numberOfFilesOpened == 1 or minBIC > BIC):
            minBIC = BIC
            selectedPrefix = BICFileName.split(".log")[0] 
            #optimalEdge = edge    
    optimalFilePrefix = sequenceFileName+".modelSelection"
#     print optimalFilePrefix
    sub.call("cp "+selectedPrefix+".log\t"+optimalFilePrefix+".log",shell=True)
    sub.call("cp "+selectedPrefix+".modelParameters\t"+optimalFilePrefix+".modelParameters",shell=True)
    sub.call("cp "+selectedPrefix+".edgeList\t"+optimalFilePrefix+".edgeList",shell=True)
    sub.call("cp "+selectedPrefix+".newick\t"+optimalFilePrefix+".newick",shell=True)

def StoreTreeSelectedViaModelSelectionForBootstrapReplicate(sequenceFileName,numberOfReplicates):
    for bootstrapReplicate in range(1,(1+numberOfReplicates)) :
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        StoreTreeSelectedViaModelSelection(bootstrapAlignmentFileName)
    
def ReadTree(treeFileName,treeFormat='edgeList',experimentName='py'):
    if treeFormat =='edgeList':
        edgeListFile = open(treeFileName,'r')
        T = Tree()        
        for line in edgeListFile:
            u_name, v_name, length = line.strip().split('\t')
            T.AddEdge(u_name, v_name, float(length))
        edgeListFile.close()
    elif treeFormat =='newick':
        devnull=open(os.devnull,'w')
        if os.path.isdir('/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/scripts/'):
            pathForNewickParserInR = '/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/scripts/parseNewickTree.R'
#             pathForRscript = '/TL/opt/bin/Rscript'
            pathForRscript = '/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/scripts/Rscript'
        elif os.path.isdir('/local/home/pk/Projects/MSTBasedForests/scripts/'):
            pathForNewickParserInR = '/local/home/pk/Projects/MSTBasedForests/scripts/parseNewickTree.R'
            pathForRscript = '/usr/local/bin/Rscript'
        elif os.path.isdir('/project/exaptation/Projects/MSTBasedForests/scripts/'):
            pathForNewickParserInR = '/project/exaptation/Projects/MSTBasedForests/scripts/parseNewickTree.R'
            pathForRscript = '/usr/local/package/bin/Rscript'            
        elif os.path.isdir('/home/kalaghat/exaptation/Projects/MSTBasedForests/scripts/'):
            pathForNewickParserInR = '/home/kalaghat/exaptation/Projects/MSTBasedForests/scripts/parseNewickTree.R'
            pathForRscript = '/usr/local/bin/Rscript'
        
        tempTreeFileName=treeFileName+'.tempTree'
        
        RCommandForParsingTrees= pathForRscript+'\t'+pathForNewickParserInR+'\t'+treeFileName+'\t'+tempTreeFileName
        print (RCommandForParsingTrees)
        sub.call(RCommandForParsingTrees,stdout=devnull,shell=True)
#         sub.call(RCommandForParsingTrees,shell=True)
        T = ReadTree(tempTreeFileName,'edgeList')
        sub.call('rm '+tempTreeFileName,stdout=devnull,shell=True)
        devnull.close()
    return T

def WriteConsensusTreeToFile(originalTreeFileName,minCladeFreq=0.7,numberOfBootstrapReplicates=100,rooted=True):
    if rooted:    
        consensusTreeFileName = originalTreeFileName.split(".fas")[0]+"_rooted_consensus_minCladeFreq_"+str(minCladeFreq)+".fas"+originalTreeFileName.split(".fas")[1]
    else:
        consensusTreeFileName = originalTreeFileName.split(".fas")[0]+"_unrooted_consensus_minCladeFreq_"+str(minCladeFreq)+".fas"+originalTreeFileName.split(".fas")[1]
    sumtreeScript = "sumtrees.py "
    if rooted:
        sumtreeScript += "--rooted "
    else:
        sumtreeScript += "--unrooted "
    sumtreeScript += "--min-clade-freq " + str(minCladeFreq)+" "
    for bootstrapReplicate in range(1,numberOfBootstrapReplicates):        
        bootstrapTreeFileName = originalTreeFileName.split(".fas")[0]+"_bootstrapReplicate_"+str(bootstrapReplicate)+".fas"+originalTreeFileName.split(".fas")[1]    
        sumtreeScript += bootstrapTreeFileName + " "
    #sumtreeScript += "-t " + originalTreeFileName + " "
    sumtreeScript += "--output-tree-format newick "
    sumtreeScript += "--suppress-annotations --output-tree-filepath " + consensusTreeFileName
    sub.call(sumtreeScript,shell=True)
 

def     ReadRootedTree(treeFileName,treeFormat='edgeList'):
    if treeFormat == 'edgeList' or treeFormat == 'edge_list':
        RT = RootedTree()
        treeFile = open(treeFileName,"r")
        for line in treeFile:
            parent_name, child_name, length = line.strip().split("\t")
            length = float(length)
            RT.AddDirectedEdge(parent_name, child_name, length)    
        treeFile.close()
        RT.SetRoot()
    elif treeFormat == 'newick':
        tree_ete3 = ete3.Tree(treeFileName)
        RT = RootedTree()
        h_ind = 0
        for node in tree_ete3.traverse('preorder'):            
            if node.name == '':
                node.name = 'h_'+str(h_ind)
                h_ind += 1
        for parent in tree_ete3.traverse('preorder'):
            for child in parent.children:
                branch_length = float(child.dist)                
                RT.AddDirectedEdge(parent.name, child.name, branch_length)
        RT.SetRoot()
        # RT = RootedTree()

        # if os.path.isdir('/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/scripts/'):
        #     pathForNewickParserInR = '/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/scripts/parseNewickTree.R'
        #     pathForRscript = '/TL/opt/bin/Rscript'
        # elif os.path.isdir('/local/home/pk/Projects/MSTBasedForests/scripts/'):
        #     pathForNewickParserInR = '/home/pk/workspace/MSTBasedForests/parseNewickTree.R'
        #     pathForRscript = '/usr/local/bin/Rscript'
        # elif os.path.isdir('/project/exaptation/Projects/MSTBasedForests/scripts/'):
        #     pathForNewickParserInR = '/project/exaptation/Projects/MSTBasedForests/scripts/parseNewickTree.R'
        #     pathForRscript = '/usr/local/package/bin/Rscript'
        # elif os.path.isdir('/home/kalaghat/exaptation/Projects/MSTBasedForests/scripts/'):
        #     pathForNewickParserInR = '/home/kalaghat/exaptation/Projects/MSTBasedForests/scripts/parseNewickTree.R'
        #     pathForRscript = '/usr/local/bin/Rscript'
        
        # tempTreeFileName=treeFileName+'.tempTree'
        
        # RCommandForParsingTrees= pathForRscript+'\t'+pathForNewickParserInR+'\t'+treeFileName+'\t'+tempTreeFileName
        # devnull=open(os.devnull,'w')        
        # # sub.call(RCommandForParsingTrees,shell=True)
        # sub.call(RCommandForParsingTrees,stdout=devnull,shell=True)
        # RT = ReadRootedTree(tempTreeFileName,'edgeList')
        # sub.call('rm '+tempTreeFileName,stdout=devnull,shell=True)
        # devnull.close()
    return RT


def ComputeMST(alignmentFileName):
    from config import scriptPath
    MSTFileName = alignmentFileName+".mstbPy_mst"
    MSTbackbone_script = scriptPath+"computeDistancesAndMST\t"+alignmentFileName+"\t"+MSTFileName+"\n"
    devnull=open(os.devnull,'w')
#     sub.call(MSTbackbone_script, shell=True)
    sub.call(MSTbackbone_script, stdout=devnull, shell=True)
    devnull.close()
    MST = ReadTree(MSTFileName)
    return MST

# def BuildGraphFromEdges(edges):
#     T = ig.Graph(0)
#     T.vs["name"]=[]
#     for vertex_i,vertex_j,length in edges:
#         for vertex in [vertex_i,vertex_j]:
#             if vertex not in T.vs["name"]:
#                 T.add_vertices(1)
#                 T.vs[T.vcount()-1]["name"]=vertex
#         vertexIndex_i = (T.vs["name"].index(vertex_i))
#         vertexIndex_j = (T.vs["name"].index(vertex_j))
#         T.add_edges([(vertexIndex_i,vertexIndex_j)])
#         T.es[T.ecount()-1]["length"] = float(length)
#     return T

# def ConvertToIgraphObject(MST_graphObj):
#     edgeWeightsHash = MST_graphObj.edgeLengths
#     edgesInMST = [[key[0],key[1],edgeWeightsHash[key]] for key in edgeWeightsHash.keys()]
#     return BuildGraphFromEdges(edgesInMST)

def ConvertToGraphObject(graph):
    G = Tree()
    edgeNameAndWeightList = map(lambda e: graph.vs[e.tuple]["name"]+[e["length"]], graph.es)
    map(lambda e: G.AddEdge(e[0], e[1], e[2]),edgeNameAndWeightList)
    return G

# def GetNewickLabelOfLeafLabeledTree(T): 
#     # Add root at the midpoint between a v and its parent
#     degrees = T.degree()
#     T.add_vertices(1)
#     T.vs[T.vcount()-1]["name"] = 'hiddenVertex_root'
#     v = degrees.index(1)
#     parentOfLeaf = T.neighbors(v)[0]
#     length = T.es[ig.Graph.get_eid(T,v,parentOfLeaf)]["length"]
#     T.add_edges([(v,T.vcount()-1)])
#     T.es[T.get_eid(v,T.vcount()-1)]["length"]=length/2  
#     T.add_edges([(parentOfLeaf,T.vcount()-1)])
#     T.es[T.get_eid(parentOfLeaf,T.vcount()-1)]["length"]=length/2
#     T.delete_edges([(v,parentOfLeaf)])
#     degrees.append(2)
#     newickLabels = {} 
#     orderedListOfVertices = []
#     numberOfTimesVertexIsVisited = {}
#     root = T.vcount()-1
#     for vertex in range(0,T.vcount()):
#         if  degrees[vertex]>1:
#             numberOfTimesVertexIsVisited[vertex]=0
#             newickLabels[vertex]='('
#         else:
#             orderedListOfVertices.append(vertex)
#             parentOfVertex = T.get_shortest_paths(root)[vertex][-2]
#             newickLabels[vertex] = T.vs[vertex]["name"]+':'+str(T.es[T.get_eid(parentOfVertex,vertex)]["length"])
#     while len(orderedListOfVertices)>0:
#         vertex = orderedListOfVertices[0]
#         if degrees[vertex]==1:
#             parentOfVertex = T.neighbors(vertex)[0]
#             newickLabels[parentOfVertex]+=newickLabels[vertex]
#             if numberOfTimesVertexIsVisited[parentOfVertex] == degrees[parentOfVertex]-2 and T.vs[parentOfVertex]["name"]!='hiddenVertex_root':
#                 newickLabels[parentOfVertex]+=')'
#                 orderedListOfVertices.append(parentOfVertex)
#             else:
#                 newickLabels[parentOfVertex]+=','
#                 numberOfTimesVertexIsVisited[parentOfVertex]+=1
#         else:
#             pathFromRootToVertex = T.get_shortest_paths(root)[vertex]
#             parentOfVertex = pathFromRootToVertex[-2]
#             newickLabels[vertex]+= T.vs[vertex]["name"]+':'+str(T.es[T.get_eid(parentOfVertex,vertex)]["length"])
#             newickLabels[parentOfVertex]+=newickLabels[vertex]
#             if T.vs[parentOfVertex]["name"]=='hiddenVertex_root' and numberOfTimesVertexIsVisited[parentOfVertex] == degrees[parentOfVertex]-1:
#                 newickLabels[parentOfVertex]+=');'
#             elif numberOfTimesVertexIsVisited[parentOfVertex] == degrees[parentOfVertex]-2:
#                 newickLabels[parentOfVertex]+=')'
#                 orderedListOfVertices.append(parentOfVertex)
#             else:
#                 newickLabels[parentOfVertex]+=','
#                 numberOfTimesVertexIsVisited[parentOfVertex]+=1
#         orderedListOfVertices.remove(vertex)
#     return newickLabels[root]

def ConvertMultifurcatingTreeToBifurcatingTree(T):
    initialDegrees = T.degree()
    numberOfLatentVertices = 1
    verticesToResolve=[]
    for vertex in range(T.vcount()):
        if initialDegrees[vertex]>3:
            verticesToResolve.append(T.vs[vertex]["name"])
        if T.vs[vertex]["name"].startswith("hiddenVertex"):
            numberOfLatentVertices+=1
    for vertexName in verticesToResolve:
        while T.degree()[T.vs["name"].index(vertexName)] > 3:
            v1, v2 = T.vs[T.neighbors(T.vs["name"].index(vertexName))[0:2]]["name"]
            newNode = 'hiddenVertexT' + str(numberOfLatentVertices)
            T.add_vertices(1)
            T.vs[T.vcount()-1]["name"]=newNode
            d_v0_h = T.es[T.get_eid(T.vs["name"].index(v1),T.vs["name"].index(vertexName))]["length"] 
            T.add_edges([(T.vs["name"].index(v1),T.vs["name"].index(newNode))])
            T.es[T.get_eid(T.vs["name"].index(v1),T.vs["name"].index(newNode))]["length"] = d_v0_h
            d_v1_h = T.es[T.get_eid(T.vs["name"].index(v2),T.vs["name"].index(vertexName))]["length"]
            T.add_edges([(T.vs["name"].index(v2),T.vs["name"].index(newNode))])
            T.es[T.get_eid(T.vs["name"].index(v2),T.vs["name"].index(newNode))]["length"] = d_v1_h
            T.add_edges([(T.vs["name"].index(vertexName),T.vs["name"].index(newNode))])
            T.es[T.get_eid(T.vs["name"].index(newNode),T.vs["name"].index(vertexName))]["length"] = 0
            T.delete_edges([(T.vs["name"].index(v1),T.vs["name"].index(vertexName))])
            T.delete_edges([(T.vs["name"].index(v2),T.vs["name"].index(vertexName))])
            numberOfLatentVertices+=1
    return T
 
def ConvertGenerallyLabeledTreeToLeafLabeledTree(T):
    degrees = T.degree()
    internalLabeledVertices=[]
    largestIdOfLatentVertex=1
    for vertex in range(T.vcount()):
        if T.vs[vertex]["name"].startswith("hiddenVertex"):
            largestIdOfLatentVertex+=1
        elif degrees[vertex]>1:
            internalLabeledVertices.append(T.vs[vertex]["name"])
        else:
            pass
    for vertexName in internalLabeledVertices:
        newNode = 'hiddenVertexT'+str(largestIdOfLatentVertex)
        largestIdOfLatentVertex+=1
        T.add_vertices(1)
        T.vs[T.vcount()-1]["name"]=newNode
        neighbors = T.vs[T.neighbors(T.vs["name"].index(vertexName))]["name"]
        for v in neighbors:
            d_v0_h = T.es[T.get_eid(T.vs["name"].index(v),T.vs["name"].index(vertexName))]["length"] 
            T.add_edges([(T.vs["name"].index(v),T.vs["name"].index(newNode))])
            T.es[T.get_eid(T.vs["name"].index(v),T.vs["name"].index(newNode))]["length"] = d_v0_h
            T.delete_edges([(T.vs["name"].index(v),T.vs["name"].index(vertexName))])
        T.add_edges([(T.vs["name"].index(vertexName),T.vs["name"].index(newNode))])
        T.es[T.get_eid(T.vs["name"].index(newNode),T.vs["name"].index(vertexName))]["length"] = 0
    return T

def WriteTree(T,outputFileName,fileFormat='edgeList'):
    treeCopy = T.copy()
    outputFile = open(outputFileName,'w')
    if fileFormat == 'edgeList':
        for e in T.es:
            i,j = e.tuple
            outputFile.write(T.vs[i]["name"]+'\t'+T.vs[j]["name"]+'\t'+str(e["length"])+'\n')
    elif fileFormat =='newick':
        leafLabeledTree = ConvertGenerallyLabeledTreeToLeafLabeledTree(treeCopy)
        binaryTree = ConvertMultifurcatingTreeToBifurcatingTree(leafLabeledTree)
#         treeInNewickFormat = GetNewickLabelOfLeafLabeledTree(binaryTree)
#         treeInNewickFormat = re.sub("hiddenVertex_?l\d+-l\d+-l\d+_\d.\d+_\d*", "", treeInNewickFormat)
#         treeInNewickFormat = re.sub("hiddenVertex_?l\d+-l\d+-l\d+", "", treeInNewickFormat)
#         treeInNewickFormat = re.sub("pseudoHiddenVertexT?_?\d+_?\d*", "", treeInNewickFormat)
#         outputFile.write(treeInNewickFormat)
    outputFile.close()

def AddExpDistEdgeLengthsToTree(inputTreeFileName,outputTreeFileName, meanEdgeLength = 0.1):
    RT = ReadRootedTree(inputTreeFileName, treeFormat = 'edgeList')
    edges = RT.edgeLengths.keys()
    for edge in edges:
        edgeLength = np.random.exponential(meanEdgeLength)
        RT.edgeLengths[edge] = edgeLength
    RT.WriteToFile(outputTreeFileName, fileFormat = 'edgeList')


def WriteRootedTreeInNewickFormat(T,vertexToPlaceRootNextTo_name,fileName):
    vertex = T.GetVertex(vertexToPlaceRootNextTo_name) 
    vertex_name = 'root'
    T.AddVertex(vertex_name)
    root = T.GetVertex(vertex_name)
    T.AddEdge(vertex_name, vertexToPlaceRootNextTo_name, 0)
    neighbors = vertex.neighbors
    for neighbor in neighbors:
        if neighbor.name.startswith('hiddenVertex'):
            w = T.GetEdgeLength(neighbor.name, vertexToPlaceRootNextTo_name)
            T.RemoveEdge(neighbor.name, vertexToPlaceRootNextTo_name)
            T.AddEdge(neighbor.name, vertex_name, w)
    leaves=[]
    for v in T.vertices.values():
        if v.degree==1:
            leaves.append(v)
    verticesInCurrentLevel = leaves
    newickLabel = {}
    for v in leaves:
        newickLabel[v]=v.name
    verticesVisited=set(leaves)
    while(len(verticesInCurrentLevel)>0):
        v = verticesInCurrentLevel[0]
        del verticesInCurrentLevel[0]
        u = list(set(v.neighbors) - verticesVisited)[0]
        w = T.GetEdgeLength(u.name,v.name)
        if u in newickLabel:
            newickLabel[u]+=','+newickLabel[v]+':'+str(w)
            if u.name!='root' and len(set(u.neighbors) - verticesVisited)==1:
                verticesInCurrentLevel.append(u)
                verticesVisited.update([u])
                newickLabel[u]+=')'
        else:
            newickLabel[u]='('+newickLabel[v]+':'+str(w)
        del newickLabel[v]    
    newickLabel[root]+=');'
    newickFile=open(fileName,'w')
    newickFile.write(newickLabel[root])
    newickFile.close()

# Read And Write Q matrix
def WriteQToFile(Q,fileName):
    QFile = open(fileName,"w")
    for i in range(4):
        for j in range(4):
            QFile.write(str(Q[i,j])+"\t")
        QFile.write("\n")

# Write scalar dictionary to file
def WriteScalarDicWithEdgeKeyToFile(scalarDictionary,fileName):
    scalarFile = open(fileName,"w")
    for edge in scalarDictionary.keys():
        scalarFile.write(edge[0]+'\t'+edge[1]+'\t')
        scalarFile.write(str(scalarDictionary[edge])+'\t')
        scalarFile.write('\n')
    scalarFile.close()

# Write scalar dictionary to file
def WriteScalarDicWithVertexKeyToFile(scalarDictionary,fileName):
    scalarFile = open(fileName,"w")
    for vertex_name in scalarDictionary.keys():
        scalarFile.write(vertex_name+'\t')
        scalarFile.write(str(scalarDictionary[vertex_name])+'\t')
        scalarFile.write('\n')
    scalarFile.close()

# Write vector dictionary to file
def WriteVectorDicWithEdgeKeyToFile(vectorDictionary,fileName):
    vectorFile = open(fileName,"w")
    vectorLength = len(vectorDictionary.values()[0])
    for edge in vectorDictionary.keys():
        vectorFile.write(edge[0]+'\t'+edge[1]+'\t')
        for i in range(vectorLength):
            vectorFile.write(str(vectorDictionary[edge][i])+'\t')
        vectorFile.write('\n')
    vectorFile.close()

# Write matrix dictionary to file
def WriteMatrixDicWithEdgeKeyToFile(matrixDictionary,fileName):
    matricesFile = open(fileName,"w")
    matrixDimension = list(matrixDictionary.values())[0].shape[0]
    for edge in matrixDictionary.keys():
        matricesFile.write(edge[0]+'\t'+edge[1]+'\t')
        for i in range(matrixDimension):
            for j in range(matrixDimension):
                matricesFile.write(str(matrixDictionary[edge][i,j])+'\t')
        matricesFile.write('\n')
    matricesFile.close()

def WriteGMM(rootProbability, transitionMatrices, GMMFileName):
    GMMFile = open(GMMFileName,"w")
    GMMFile.write("Root probability\n")
    for i in range(4):
        GMMFile.write("P("+str(DNA[i])+")\t")
    GMMFile.write("\n")
    for i in range(4):
        GMMFile.write(str(rootProbability[i])+"\t")
    GMMFile.write("\n")
    GMMFile.write("Transition matrices\n")
    GMMFile.write("Vertex_from\tVertex_to\t")
    for p in range(4):
        for c in range(4):
            GMMFile.write("P("+str(DNA[p])+"->"+str(DNA[c])+")\t")
    GMMFile.write("\n")
    for (u_name, v_name) in transitionMatrices.keys():
        GMMFile.write(u_name+"\t"+v_name+"\t")
        P = transitionMatrices[(u_name,v_name)]
        for p in range(4):
            for c in range(4):
                GMMFile.write(str(P[p,c])+"\t")
        GMMFile.write("\n")
    GMMFile.close()

def ReadRootProbability(GMMFileName):
    rootProbability = [0.0] * 4
    GMMFile = open(GMMFileName,"r")
    lineNum = 1
    for line in GMMFile:
        if lineNum == 3:
            rootProbability = line.strip().split("\t") 
            break
        lineNum += 1
    GMMFile.close()
    rootProbability = map(lambda x: float(x),rootProbability)
    return (rootProbability)

def ReadTransitionMatrices(GMMFileName):
    transitionMatrices = {}
    GMMFile = open(GMMFileName,"r")
    lineNum = 1
    u_name = ""
    v_name = ""    
    for line in GMMFile:
        if lineNum > 5:
            lineSplit = line.strip().split("\t")
            if (len(lineSplit) != 18):
                print ("check number of elements in split line")
            u_name = lineSplit[0]
            v_name = lineSplit[1]
            P = array([[0.0]*4]*4)            
            P[0,0] = float(lineSplit[2])
            P[0,1] = float(lineSplit[3])
            P[0,2] = float(lineSplit[4])
            P[0,3] = float(lineSplit[5])
            P[1,0] = float(lineSplit[6])
            P[1,1] = float(lineSplit[7])
            P[1,2] = float(lineSplit[8])
            P[1,3] = float(lineSplit[9])
            P[2,0] = float(lineSplit[10])
            P[2,1] = float(lineSplit[11])
            P[2,2] = float(lineSplit[12])
            P[2,3] = float(lineSplit[13])
            P[3,0] = float(lineSplit[14])
            P[3,1] = float(lineSplit[15])
            P[3,2] = float(lineSplit[16])
            P[3,3] = float(lineSplit[17])            
            transitionMatrices[(u_name,v_name)] = P
        lineNum += 1    
    GMMFile.close()
    return (transitionMatrices)

def ReadQFromFile(fileName):
    Q = array([[0.0]*4]*4)
    QFile = open(fileName,"r")
    row = 0
    for line in QFile:
        Qline = line.strip().split("\t") 
        for i in range(4):
            Q[row,i] = float(Qline[i]) 
        row +=1
    return Q


