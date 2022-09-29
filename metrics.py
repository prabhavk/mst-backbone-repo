import copy
from classDeclarationsAndFunctions import RootedTree
from math import log
from fileIO import ReadAlignment, ReadRootedTree
from config import projectPath
def GetLogLikelihood(logFileName,model):
    logLikelihood = 0.0
    if model == "GM":
        stringForSplittingLine = "Log likelihood of fitted general Markov model is"
    elif model == "UNR":
        stringForSplittingLine = "Log-likelihood of the tree:"
    elif model == "GTR":
        stringForSplittingLine = "Final LogLikelihood:"
    logFile = open(logFileName,"r")
    for line in logFile:
        if line.startswith(stringForSplittingLine):
            logLikelihood = float(line.strip().split(stringForSplittingLine)[1].split("(")[0].strip())
            break
    logFile.close()
    return (logLikelihood)

# def GetSplitsForFullyLabeledTree(edgeList):
#     splits = ()
#     vertList1.sort()
#         vertList2.sort()
#         if vertList1 < vertList2:
#             splits.append(''.join(vertList1) + '-' + ''.join(vertList2))
#         else: 
#             splits.append(''.join(vertList2) + '-' + ''.join(vertList1))

def GetScoresForMSTBackbonePaper(sequenceFileName,model):
    if model == "GM":
        logFileName = sequenceFileName + ".mstbackboneSEM_log"
    elif model == "UNR":
        logFileName = sequenceFileName + ".iqtree"
    elif model == "GTR":
        logFileName = sequenceFileName + ".raxml.log"
    logLikelihood = GetLogLikelihood(logFileName, model)
    alignment = ReadAlignment(sequenceFileName)
    nLeaves = len(alignment)
    numSamples = len(alignment.values()[0])
    numFreeParameters = GetNumFreePar(nLeaves, model)
    BIC = GetBIC(logLikelihood, numSamples, numFreeParameters)
    AIC = GetAIC(logLikelihood, numFreeParameters)
    AICc = GetAICc(AIC, numSamples, numFreeParameters)    
    tupleToReturn = (logLikelihood,BIC,AIC,AICc)
    return (tupleToReturn)

def GetNumFreePar(nLeaves,model):
    if model == "GM":
        k = 3 # root prob
        nEdges = (2*nLeaves) - 2
        k += 12 * nEdges
    elif model == "UNR":
        k = 11 # calibrated unrestricted rate matrix
        nEdges = (2*nLeaves) - 2
        k += nEdges
    elif model == "GTR":
        k = 8 # calibrated rate matrix s.t. PiQ is symmetric
        nEdges = (2*nLeaves) - 3
        k += nEdges
    return k

def GetBIC(logLikelihood,numSamples,numFreeParameters):
    n = float(numSamples)
    k = float(numFreeParameters)
    BIC = -(2*float(logLikelihood)) + k*log(n)
    return BIC

def GetAIC(logLikelihood,numFreeParameters):
    k = float(numFreeParameters)
    AIC = -(2*float(logLikelihood)) + (2*k)
    return AIC

def GetAICc(AIC,numSamples,numFreeParameters):
    n = float(numSamples)
    k = float(numFreeParameters)
    AICc = AIC + ((2*k)*(k+1))/(max((n-k-1),1.0))
    
    return (AICc)

def GetPrecisionRecallRFAndBLD(trueTreeOriginal,estimatedTreeOriginal):
    trueTree = copy.deepcopy(trueTreeOriginal)
    estimatedTree = copy.deepcopy(estimatedTreeOriginal)
    for v in range(0,trueTree.vcount()):
        if trueTree.vs[v]["name"].startswith('hiddenVertex'):
            trueTree.vs[v]["name"]="hiddenVertex"
    for v in range(0,estimatedTree.vcount()):
        if estimatedTree.vs[v]["name"].startswith('hiddenVertex'):
            estimatedTree.vs[v]["name"]="hiddenVertex"
    trueSplits2Lengths={}
    for e in trueTree.es:
        graphCopy = copy.copy(trueTree)
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            trueSplits2Lengths[''.join(vertList1) + '-' + ''.join(vertList2)] = e["length"]
        else: 
            trueSplits2Lengths[''.join(vertList2) + '-' + ''.join(vertList1)] = e["length"]
    estimatedSplits2Lengths={}
    for e in estimatedTree.es:
        graphCopy = copy.copy(estimatedTree)
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            estimatedSplits2Lengths[''.join(vertList1) + '-' + ''.join(vertList2)] = e["length"]
        else: 
            estimatedSplits2Lengths[''.join(vertList2) + '-' + ''.join(vertList1)] = e["length"]
    trueSplits = trueSplits2Lengths.keys()
    estimatedSplits = estimatedSplits2Lengths.keys()
    commonSplits = list(set(estimatedSplits)&set(trueSplits))
    splitsOnlyInTrueTree = set(trueSplits)- set(commonSplits)
    splitsOnlyInEstimatedTree = set(estimatedSplits)- set(commonSplits)
    numberOfTrueSplitsRecovered = len(commonSplits)
    precision = numberOfTrueSplitsRecovered/float(len(estimatedSplits))
    recall = numberOfTrueSplitsRecovered/float(len(trueSplits2Lengths))
    RF = 1 - numberOfTrueSplitsRecovered/float(len((set(estimatedSplits)|set(trueSplits))))
    BLD = 0
    for split in splitsOnlyInTrueTree:
        BLD+=trueSplits2Lengths[split]**2
    for split in splitsOnlyInEstimatedTree:
        BLD+=estimatedSplits2Lengths[split]**2
    for split in commonSplits:
        BLD+=(trueSplits2Lengths[split]-estimatedSplits2Lengths[split])**2
    return(precision, recall, RF, BLD)

def GetFractionOfDNAChars(sequence):
    frac = 0;
    sequence.upper()
    for char in sequence:
        if char in ["A","C","G","T"]:
            frac += 1
    frac /= float(len(sequence))
    return (frac)

def GetGCComp(sequence):
    frac = 0;
    sequence.upper()
    nonAmbChars = 0;
    for char in sequence:
        if char in ["C","G"]:
            frac += 1
            nonAmbChars += 1
        elif char in ["A","T"]:
            nonAmbChars += 1
    frac /= float(nonAmbChars)
    return (frac)

def GetMaxGCCompDiff(alignment):
    gcCompRange = [GetGCComp(seq) for seq in alignment.values()]
    return (max(gcCompRange) - min(gcCompRange))


def WriteFullyLabeledRootedPrecisionRecallAndRFForExpPhyloGMForBootstrapReps(modelTree,sequenceFileName,numberOfReplicates):
    fileNameSuffix = sequenceFileName.split("/")[-1].split(".fas")[0]
    rootedReconstuctionAccuracy_fileName = projectPath + 'results/reconstructionAccuracy/bootstrap_rooted_PR_RE_RF_'+fileNameSuffix
    reconAccuracyFile = open(rootedReconstuctionAccuracy_fileName,"w")
    for bootstrapReplicate in range(1,(1+numberOfReplicates)) :
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        treeFileName = bootstrapAlignmentFileName + ".edgeList"
        estimatedTree = ReadRootedTree(treeFileName, "edgeList")
        pr, re, rf = ComputePrecisionRecallAndRFforRootedTrees(modelTree, estimatedTree)
        reconAccuracyFile.write(str(pr)+"\t"+str(re)+"\t"+str(rf)+"\n")
    reconAccuracyFile.close()

def WriteFullyLabeledRootedPrecisionRecallAndRFForModelSelectedTreeForBootstrapReps(modelTree,sequenceFileName,numberOfReplicates):
    fileNameSuffix = sequenceFileName.split("/")[-1].split(".fas")[0]+"_modelSelection"
    rootedReconstuctionAccuracy_fileName = projectPath + 'results/reconstructionAccuracy/bootstrap_rooted_PR_RE_RF_'+fileNameSuffix
    reconAccuracyFile = open(rootedReconstuctionAccuracy_fileName,"w")
    for bootstrapReplicate in range(1,(1+numberOfReplicates)) :
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        treeFileName = bootstrapAlignmentFileName + ".modelSelection.edgeList"
        estimatedTree = ReadRootedTree(treeFileName, "edgeList")
        pr, re, rf = ComputePrecisionRecallAndRFforRootedTrees(modelTree, estimatedTree)
        reconAccuracyFile.write(str(pr)+"\t"+str(re)+"\t"+str(rf)+"\n")
    reconAccuracyFile.close()


def WriteRootedPrecisionRecallAndRFForExpPhyloGMForBootstrapReps(modelTree,sequenceFileName,numberOfReplicates):
    fileNameSuffix = sequenceFileName.split("/")[-1].split(".fas")[0]
    rootedReconstuctionAccuracy_fileName = projectPath + 'results/reconstructionAccuracy/bootstrap_rooted_PR_RE_RF_'+fileNameSuffix
    reconAccuracyFile = open(rootedReconstuctionAccuracy_fileName,"w")
    for bootstrapReplicate in range(1,(1+numberOfReplicates)) :
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        treeFileName = bootstrapAlignmentFileName + ".edgeList"
        estimatedTree = ReadRootedTree(treeFileName, "edgeList")
        pr, re, rf = ComputePrecisionRecallAndRFforNonTrivialClusters(modelTree, estimatedTree)
        reconAccuracyFile.write(str(pr)+"\t"+str(re)+"\t"+str(rf)+"\n")
    reconAccuracyFile.close()

def WriteRootedPrecisionRecallAndRFForModelSelectedTreeForBootstrapReps(modelTree,sequenceFileName,numberOfReplicates):
    fileNameSuffix = sequenceFileName.split("/")[-1].split(".fas")[0]+"_modelSelection"
    rootedReconstuctionAccuracy_fileName = projectPath + 'results/reconstructionAccuracy/bootstrap_rooted_PR_RE_RF_'+fileNameSuffix
    reconAccuracyFile = open(rootedReconstuctionAccuracy_fileName,"w")
    for bootstrapReplicate in range(1,(1+numberOfReplicates)) :
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        treeFileName = bootstrapAlignmentFileName + ".modelSelection.edgeList"
        estimatedTree = ReadRootedTree(treeFileName, "edgeList")
        pr, re, rf = ComputePrecisionRecallAndRFforNonTrivialClusters(modelTree, estimatedTree)
        reconAccuracyFile.write(str(pr)+"\t"+str(re)+"\t"+str(rf)+"\n")
    reconAccuracyFile.close()
    

def WriteFullyLabeledUnrootedPrecisionRecallAndRFForExpPhyloGMForBootstrapReps(modelTree,sequenceFileName,numberOfReplicates):
    fileNameSuffix = sequenceFileName.split("/")[-1].split(".fas")[0]
    rootedReconstuctionAccuracy_fileName = projectPath + 'results/reconstructionAccuracy/bootstrap_unrooted_PR_RE_RF_'+fileNameSuffix
    reconAccuracyFile = open(rootedReconstuctionAccuracy_fileName,"w")
    for bootstrapReplicate in range(1,(1+numberOfReplicates)) :
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        treeFileName = bootstrapAlignmentFileName + ".edgeList"
        estimatedTree = ReadRootedTree(treeFileName, "edgeList")
        pr, re, rf = ComputePrecisionRecallAndRFforAllSplits(modelTree, estimatedTree)
        reconAccuracyFile.write(str(pr)+"\t"+str(re)+"\t"+str(rf)+"\n")
    reconAccuracyFile.close()

def WriteUnrootedPrecisionRecallAndRFForExpPhyloGMForBootstrapReps(modelTree,sequenceFileName,numberOfReplicates):
    fileNameSuffix = sequenceFileName.split("/")[-1].split(".fas")[0]
    rootedReconstuctionAccuracy_fileName = projectPath + 'results/reconstructionAccuracy/bootstrap_unrooted_PR_RE_RF_'+fileNameSuffix
    reconAccuracyFile = open(rootedReconstuctionAccuracy_fileName,"w")
    for bootstrapReplicate in range(1,(1+numberOfReplicates)) :
        bootstrapAlignmentFileName = sequenceFileName.split(".fas")[0] + '_bootstrapReplicate_'+str(bootstrapReplicate)+'.fas'
        treeFileName = bootstrapAlignmentFileName + ".edgeList"
        estimatedTree = ReadRootedTree(treeFileName, "edgeList")
        pr, re, rf = ComputePrecisionRecallAndRFforNonTrivialSplits(modelTree, estimatedTree)
        reconAccuracyFile.write(str(pr)+"\t"+str(re)+"\t"+str(rf)+"\n")
    reconAccuracyFile.close()
    
    
def GetPrecisionRecallAndRF(trueTreeOriginal,estimatedTreeOriginal):
    trueTree = copy.deepcopy(trueTreeOriginal)
    estimatedTree = copy.deepcopy(estimatedTreeOriginal)
    for v in range(0,trueTree.vcount()):
        if trueTree.vs[v]["name"].startswith('hiddenVertex'):
            trueTree.vs[v]["name"]="hiddenVertex"
    for v in range(0,estimatedTree.vcount()):
        if estimatedTree.vs[v]["name"].startswith('hiddenVertex'):
            estimatedTree.vs[v]["name"]="hiddenVertex"
    trueSplits=[]
    for e in trueTree.es:
        graphCopy = copy.copy(trueTree)
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            trueSplits.append(''.join(vertList1) + '-' + ''.join(vertList2))
        else: 
            trueSplits.append(''.join(vertList2) + '-' + ''.join(vertList1))
    estimatedSplits=[]
    for e in estimatedTree.es:
        graphCopy = copy.copy(estimatedTree)
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            estimatedSplits.append(''.join(vertList1) + '-' + ''.join(vertList2))
        else: 
            estimatedSplits.append(''.join(vertList2) + '-' + ''.join(vertList1))
    numberOfTrueSplitsRecovered = len(set(estimatedSplits)&set(trueSplits)) 
    precision = numberOfTrueSplitsRecovered/float(len(estimatedSplits))
    recall = numberOfTrueSplitsRecovered/float(len(trueSplits))
    RF = 1 - numberOfTrueSplitsRecovered/float(len((set(estimatedSplits)|set(trueSplits))))
    return(precision, recall, RF)

# def ComputeKLDivergenceFromDiscreteTimeModelToDiscreteTimeModel(sequences,rootProb_model,rootProb_estimated,RT_model,RT_estimated,transitionMatrices_model,transitionMatrices_estimated):
#     pass
# 
# def ComputeKLDivergenceFromContinuousTimeModelToDiscreteTimeModel(sequences,rootProb_model,rootProb_estimated,RT_model,RT_estimated,transitionMatrices_model,rateMatrix_estimated):
#     pass

def ComputePrecisionRecallAndRFforRootedTrees(RT_model, RT_estimate):
    clusters_model = RT_model.GetAllNonTrivialClusters()
    clusters_estimated = RT_estimate.GetAllNonTrivialClusters()
    # clusters_model = RT_model.GetAllNonRootClusters()
    # clusters_estimate = RT_estimate.GetAllNonRootClusters()
    numberOfClustersInCommon = len(clusters_model & clusters_estimated)
    totalNumberOfDistinctClusters = len(clusters_model | clusters_estimated)
    RF = 1 - numberOfClustersInCommon/float(totalNumberOfDistinctClusters)
    precision = numberOfClustersInCommon/float(len(clusters_estimated))
    recall = numberOfClustersInCommon/float(len(clusters_model))
    return (precision, recall, RF)

def ComputePrecisionRecallAndRFforNonTrivialClusters(RT_model, RT_estimate):
    clusters_model_all = RT_model.GetAllNonRootClusters()
    clusters_model = set([])
    for cluster in clusters_model_all:
        if len(cluster.split(","))>1:
            clusters_model.update([cluster])
    clusters_estimate_all = RT_estimate.GetAllNonRootClusters()
    clusters_estimate = set([])
    for cluster in clusters_estimate_all:
        if len(cluster.split(","))>1:
            clusters_estimate.update([cluster])
    numberOfClustersInCommon = len(clusters_model & clusters_estimate)
    totalNumberOfDistinctClusters = len(clusters_model | clusters_estimate)
    RF = 1 - numberOfClustersInCommon/float(totalNumberOfDistinctClusters)
    precision = numberOfClustersInCommon/float(len(clusters_estimate))
    recall = numberOfClustersInCommon/float(len(clusters_model))
    return (precision, recall, RF)

def ComputePrecisionRecallAndRFforNonTrivialSplits(RT_model, RT_estimate):
    nonTrivialSplits_model = RT_model.GetAllNonTrivialSplits()
    nonTrivialSplits_estimate = RT_estimate.GetAllNonTrivialSplits()
    pr_unrooted, re_unrooted, rf_unrooted = ComputePrecisionRecallAndRFGivenSplits(nonTrivialSplits_model, nonTrivialSplits_estimate)
    return (pr_unrooted, re_unrooted, rf_unrooted)


def ComputePrecisionRecallAndRFforAllSplits(RT_model, RT_estimate):
    allSplits_model = RT_model.GetAllSplits()
    allSplits_estimate = RT_estimate.GetAllSplits()
    pr_unrooted, re_unrooted, rf_unrooted = ComputePrecisionRecallAndRFGivenSplits(allSplits_model, allSplits_estimate)
    return (pr_unrooted, re_unrooted, rf_unrooted)
        
def ComputePrecisionRecallAndRFForUnrootedTrees(T_model,T_estimate):    
    splits_model = set(T_model.GetAllSplitsInTree())
    splits_estimate = set(T_estimate.GetAllSplitsInTree())
    numberOfSplitsInCommon = len(splits_model & splits_estimate)
    totalNumberOfDistinctSplits = len(splits_model | splits_estimate)
    RF = 1 - numberOfSplitsInCommon/float(totalNumberOfDistinctSplits)
    precision = numberOfSplitsInCommon/float(len(splits_estimate))
    recall = numberOfSplitsInCommon/float(len(splits_model))
    return precision, recall, RF

def ComputePrecisionRecallAndRFGivenSplits(splits_model,splits_estimate):
    numberOfSplitsInCommon = len(splits_model & splits_estimate)
    totalNumberOfDistinctSplits = len(splits_model | splits_estimate)
    RF = 1 - numberOfSplitsInCommon/float(totalNumberOfDistinctSplits)
    precision = numberOfSplitsInCommon/float(len(splits_estimate))
    recall = numberOfSplitsInCommon/float(len(splits_model))
    return precision, recall, RF

def GetSplits(RT,returnEdgeToSplitMapping=False):
    treeCopy = copy.deepcopy(RT)
    for v in range(0,treeCopy.vcount()):
        if treeCopy.vs[v]["name"].startswith('hiddenVertex'):
            treeCopy.vs[v]["name"]="hiddenVertex"
    splits=[]
    if returnEdgeToSplitMapping:
        edges = []
    for e in treeCopy.es:
        graphCopy = copy.copy(treeCopy)
        if returnEdgeToSplitMapping:
            edges.append('-'.join(sorted(RT.vs[e.tuple[0],e.tuple[1]]['name'])))
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            splits.append(''.join(vertList1) + '-' + ''.join(vertList2))
        else: 
            splits.append(''.join(vertList2) + '-' + ''.join(vertList1))
    if returnEdgeToSplitMapping:
        return splits, edges
    else:
        return splits

def GetSplitInducedByEdgeInMST(MST,e):
    graphCopy = copy.copy(MST)
    graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
    graphLets = graphCopy.decompose()
    vertList1 = graphLets[0].vs["name"]
    vertList2 = graphLets[1].vs["name"]
    vertList1.sort()
    vertList2.sort()
    if vertList1 < vertList2:
        split = ''.join(vertList1) + '-' + ''.join(vertList2)
    else: 
        split = ''.join(vertList2) + '-' + ''.join(vertList1)
    
    return split

def GetRFDist(T1,T2):
    splitsInT1 = set(T1.GetAllSplitsInTree())
    splitsInT2 = set(T2.GetAllSplitsInTree())
    numberOfSplitsInCommon =  len(splitsInT1 & splitsInT2)
    totalNumberOfUniqueSplits = len(splitsInT1 | splitsInT2)
    return 1 - numberOfSplitsInCommon/float(totalNumberOfUniqueSplits)  
    

def GetRFDistUsingIgraphModule(treeTrueOrig,treeEstimateOrig):
    treeTrue = copy.deepcopy(treeTrueOrig)
    treeEstimate = copy.deepcopy(treeEstimateOrig)
    for v in range(0,treeTrue.vcount()):
        if treeTrue.vs[v]["name"].startswith('hiddenVertex'):
            treeTrue.vs[v]["name"]="hiddenVertex"
    for v in range(0,treeEstimate.vcount()):
        if treeEstimate.vs[v]["name"].startswith('hiddenVertex'):
            treeEstimate.vs[v]["name"]="hiddenVertex"
    splitListTrue=[]
    for e in treeTrue.es:
        graphCopy = copy.copy(treeTrue)
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            splitListTrue.append(''.join(vertList1) + '-' + ''.join(vertList2))
        else:
            splitListTrue.append(''.join(vertList2) + '-' + ''.join(vertList1))
    splitListEstimate=[]
    for e in treeEstimate.es:
        graphCopy = copy.copy(treeEstimate)
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            splitListEstimate.append(''.join(vertList1) + '-' + ''.join(vertList2))
        else:
            splitListEstimate.append(''.join(vertList2) + '-' + ''.join(vertList1))
    jaccardIndex = len(set(splitListEstimate)&set(splitListTrue))/float(len((set(splitListEstimate)|set(splitListTrue))))
    return(1-jaccardIndex)

def CountNumberOfCherriesInBinaryTree(RT):
    degree = RT.degree()
    leaves=[]
    for x in range(RT.vcount()):
        if degree[x]==1:
            leaves.append(x)
    numberOfLeavesConnectedToInternalVertex=[0]*(len(leaves)-2)
    for v in leaves:
        neighbors = RT.neighborhood(v)
        neighbors.remove(v)
        numberOfLeavesConnectedToInternalVertex[neighbors[0]-len(leaves)]+=1
    numberOfCherries=0
    for n in numberOfLeavesConnectedToInternalVertex:
        if n==2:
            numberOfCherries+=1            
    return numberOfCherries

# def ComputeKLDivergenceForMSTBackbone(ModelGMMFileName,EstimatedGMMFileName,leafSequences):
#     RT_model = RootedTree()
#     ModelGMMFile = open(ModelGMMFileName,"r")
#     lineNo = 1
#     for line in ModelGMMFile:
#         if lineNo == 3: 
#             RT_model.rootProb = [float(x) for x in line.strip().split("\t")]
#         if lineNo > 5:
#             lineSplit = line.strip().split("\t")
#             p_name = lineSplit[0]
#             c_name = lineSplit[1]
#             RT_model.AddDirectedEdge(p_name, c_name)
#             P = array([[0.0]*4]*4)
#             ind = 2
#             for dna_from in range(4):
#                 for dna_to in range(4):
#                     P[dna_from][dna_to] = lineSplit[ind]
#                     ind += 1
#             RT_model.transitionMatrices[(p_name,c_name)] = P
#         lineNo += 1
#     ModelGMMFile.close()
#     RT_model.SetRoot()
#     RT_estimated = RootedTree()
#     EstimatedGMMFile = open(EstimatedGMMFileName,"r")
#     lineNo = 1
#     for line in EstimatedGMMFile:
#         if lineNo == 3:
#             RT_estimated.rootProb = [float(x) for x in line.strip().split("\t")]
#         if lineNo > 5:
#             lineSplit = line.strip().split("\t")
#             p_name = lineSplit[0]
#             c_name = lineSplit[1]
#             RT_estimated.AddDirectedEdge(p_name, c_name)
#             P = array([[0.0]*4]*4)
#             ind = 2
#             for dna_from in range(4):
#                 for dna_to in range(4):
#                     P[dna_from][dna_to] = lineSplit[ind]
#                     ind += 1
#             RT_estimated.transitionMatrices[(p_name,c_name)] = P
#         lineNo += 1
#     EstimatedGMMFile.close()
#     RT_estimated.SetRoot()
#     RT_model.AddLeafSequences(leafSequences)
#     RT_estimated.AddLeafSequences(leafSequences)
#     sequenceLength = len(leafSequences.values()[0])
#     KLDivergence = 0    
#     for i in range(sequenceLength):
#         # P(i) Compute likelihood of site i under model tree
#         P = RT_model.ComputeLikelihoodForSite(i)
#         # Q(i) Compute likelihood of site i under estimated tree
#         Q = RT_estimated.ComputeLikelihoodForSite(i)
#         KLDivergence += P*(log(P)-log(Q))
#     return KLDivergence
