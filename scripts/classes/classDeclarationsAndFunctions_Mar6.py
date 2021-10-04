from operator import itemgetter
from numpy.random import choice
from distanceCalculator import hamDist, GetBaseFreq
from numpy import array
from math import log
import cmath as cmath
from scipy.linalg import logm
from scipy.optimize._minimize import minimize_scalar
from scipy.optimize import minimize 
from MarkovModels import NormalizeQMatrix,\
    ComputeQMatrixFromVectorOfOffDiagonalElements,\
    ComputeProbabilityMatrixUsingMatrixExponentiation,\
    ComputeQVectorOfOffDiagonalElementsFromQMatrix,\
    GetMaxLikEstimateOfTransitionMatrix
from config import tempPath

DNA=["A","C","G","T"]

def EMForMapWithFixedQ(RT):
    RT.PerformEMWithFixedQ()
    return RT

def EMForMapWithEstimatedQ(RT):
    RT.PerformEMWithEstimatedQ()
    return RT

def EMForMapForGMM(RT):
    RT.PerformEMForGMMUsingPatternSavings()
    return RT

def poolFuncForComputingLogLikForRooting(RT):
    logLik = RT.ReturnLogLikForRootingAtEdge()
    return logLik

def EMForMapWithFixedStatesAndEdgeLengthsAndEstimatedQ(RT):
    RT.PerformEMWithFixedStatesAndEdgeLengthsAndEstimatedQ()
    return RT

# rep is short of representative
class Vertex:
    def __init__(self,name):
    # union-find
        self.rep = self
        self.comp = name
    # union-find and ML-VR-MST 
        self.rank = 0
    # ML-VR-MST
        self.deltaMax = 0
    # general graph
        self.name = name
        self.degree=0
        self.neighbors = []
        self.verticesInSubtree=set([name])
        self.timesVisited=0
        self.labeled = False
    def AddNeighbor(self,neighbor):
        self.neighbors.append(neighbor)
        self.degree+=1
    def RemoveNeighbor(self,neighbor):
        self.neighbors.remove(neighbor)
        self.degree-=1        
    def GetNeighbors(self):
        return self.neighbors

class Graph:
    def __init__(self):
        self.vertices={}
        self.edgeLengths={}
    def AddVertex(self,u_name):
        self.vertices[u_name] = Vertex(u_name)
    def GetVertex(self,u_name):
        return self.vertices[u_name]
    def AddEdge(self,u_name,v_name,w):
        if u_name not in self.vertices:
            self.AddVertex(u_name)
        if v_name not in self.vertices:
            self.AddVertex(v_name)
        u = self.vertices[u_name]
        v = self.vertices[v_name]
        if Find(u) != Find(v):
            Union(u,v)
        u.AddNeighbor(v)
        v.AddNeighbor(u)
        self.edgeLengths[tuple(sorted((u_name,v_name)))] = w
    def GetEdgeLength(self,u_name,v_name):
        return self.edgeLengths[tuple(sorted((u_name,v_name)))]
    def GetConnectedComponents(self):
        components = {}
        for u in self.vertices.values():
            if Find(u) in components:
                components[Find(u)].append(u.name)
            else:
                components[Find(u)] = [u.name]
        return components        

class Tree:
    def __init__(self):
        self.vertices={}
        self.vertexCount = 0
        self.edgeLengths={}
        self.sequences = {}
        self.sequenceLength = 0
        self.siteArrangedBlockwiseForEachDistinctPattern = []
        self.firstSiteOfEachDistinctPattern = []
        self.patternRepeats = []
        self.treeDepth = 16
        self.QForEdge = {}
        self.labeled_vertices_names = []
        
    def AddVertex(self,v_name):
        if not v_name in self.vertices.keys(): 
            self.vertices[v_name] = Vertex(v_name)
            self.vertexCount += 1
        
    def GetVertex(self,v_name):
        return self.vertices[v_name]
    def ContainsVertex(self,v_name):
        return v_name in self.vertices.keys()
    def RemoveVertex(self,v_name):
        if v_name in self.vertices.keys():
            v = self.vertices[v_name]
            for u in v.neighbors:
                del self.edgeLengths[tuple(sorted((u.name,v_name)))]
                u.RemoveNeighbor(v)
            del self.vertices[v_name]
            self.vertexCount -= 1
            if v_name in self.sequences.keys():
                del self.sequences[v_name]
        else:
            return 'Error: The vertex that is to be removed is not in the graph'
    def AddPendantEdgesToLabeledAncestors(self): 
        for v_name in self.vertices.keys():
            v = self.vertices[v_name]:
            if v.degree == 1:
                pass

    def GetSplits(self):
        leaves = self.GetLeaves()[:]
        splits = {}
        for v_name in self.vertices.keys():
            v = self.vertices[v_name]
            pass # TO DO

    def RenameNonLeafVertices(self,h_ind):
        nonLeafVertices = self.GetNonLeafVertices()
        oldNameToNewName = {}
        for v in nonLeafVertices:
            h_ind += 1
            oldNameToNewName[v.name] = "h_"+str(h_ind)                            
        for v_name in oldNameToNewName.keys():
            v = self.GetVertex(v_name)
            v.name = oldNameToNewName[v_name]
            self.vertices[v.name] = v
        oldEdgesNames = self.edgeLengths.keys()
        for u_name, v_name in oldEdgesNames:
            t = self.GetEdgeLength(u_name, v_name)
            self.RemoveEdge(u_name, v_name)
            u = self.GetVertex(u_name)
            v = self.GetVertex(v_name)
            self.AddEdge(u.name, v.name,t)                        
        for v_name in oldNameToNewName.keys():
            del self.vertices[v_name]        
        # rename sequences dictionary
        for v_name in oldNameToNewName.keys():
            if v_name in self.sequences.keys():
                self.sequences[oldNameToNewName[v_name]] = self.sequences[v_name] 
                del self.sequences[v_name]
                 
    def AddEdge(self,u_name,v_name,w):
        if u_name not in self.vertices.keys():
            self.AddVertex(u_name)
        if v_name not in self.vertices.keys():
            self.AddVertex(v_name)
        u = self.vertices[u_name]
        v = self.vertices[v_name]
        # set labeled vertices
        if not u_name.startswith("h_"):
            u.labeled = True
        if not v_name.startswith("h_"):
            v.labeled = True
        u.AddNeighbor(v)
        v.AddNeighbor(u)
        self.edgeLengths[tuple(sorted((u_name,v_name)))] = w
    
    def AddEdges(self,edgeLengthsDic):
        for (u_name, v_name), t in edgeLengthsDic.items():
            self.AddEdge(u_name, v_name, t)
            
    def ContainsEdge(self,u_name,v_name):
        return tuple(sorted((u_name,v_name))) in self.edgeLengths
    
    def GetEdgeLength(self,u_name,v_name):
        return self.edgeLengths[tuple(sorted((u_name,v_name)))]
    
    def UpdateEdgeLength(self,u_name,v_name,length):
        self.edgeLengths[tuple(sorted((u_name,v_name)))] = length
        
    def RemoveEdge(self,u_name,v_name):
        del self.edgeLengths[tuple(sorted((u_name,v_name)))]
        u = self.vertices[u_name]
        v = self.vertices[v_name]
        u.RemoveNeighbor(v)
        v.RemoveNeighbor(u)
        
    def GetLeaves(self):
        leaves=[]
        for v in self.vertices.values():
            if v.degree == 1:
                leaves.append(v)
        return leaves
    
    def GetNonLeafVertices(self):
        nonLeaves=[]
        for v in self.vertices.values():
            if v.degree != 1:
                nonLeaves.append(v)
        return nonLeaves
    
    def GetLeafNames(self):
        leaves=[]
        for v in self.vertices.values():
            if v.degree == 1:
                leaves.append(v.name)
        return leaves
    
    def AddLeafSequences(self,sequences):
        if len(sequences) > 0:
            self.sequences.update(sequences)
            self.sequenceLength = len(sequences.values()[0])
            self.IdentifyRepeatedSitePatterns()
    
    def AddSequences(self,sequences):        
        if len(sequences) > 0:
            self.sequences.update(sequences)
            if self.sequenceLength == 0:
                self.sequenceLength = len(sequences.values()[0])            
        
    def GetSplitsForGenerallyLabeledTree(self):
        verticesToVisit = self.GetLeaves()[:]
        allVertexNames = set([v.name for v in self.vertices.values()])
        for v in self.vertices.values():
            if v.labeled:
                v.verticesInSubtree = set([v.name])
            else:
                v.verticesInSubtree = set([])
        while len(verticesToVisit)>0:
            v = verticesToVisit.pop()
            v.timesVisited+=1
            for n in v.neighbors:
                n.timesVisited += 1
                if n.degree - n.timesVisited > 1:
                    n.verticesInSubtree.update(v.verticesInSubtree)
                elif n.degree - n.timesVisited == 1:
                    n.verticesInSubtree.update(v.verticesInSubtree)
                    verticesToVisit.append(n)
        splitsPerEdge = []
        for (u_name, v_name) in self.edgeLengths.keys():
            u = self.GetVertex(u_name)
            v = self.GetVertex(v_name)
            leftSet = set(u.verticesInSubtree)&set(v.verticesInSubtree)
            leftList = sorted(list(leftSet))
            rightList = sorted(list(allVertexNames-leftSet))
            left = ','.join(leftList)
            right = ','.join(rightList)
            split = "-".join(sorted([left,right]))
            splitsPerEdge.append(split)
        return set(splitsPerEdge)
    def ComputeRFDistance(self,tree):
        self_splits = self.GetSplitsForGenerallyLabeledTree()
        tree_splits = tree.GetSplitsForGenerallyLabeledTree()
        assert(len(tree_splits) == len(self_splits))        
        
        list_1 = self_splits[0].split("-")
        list_2 = list_1[0].split(",") + list_1[1].split(",")
        self_labeled_names = set(list_2)
        list_1 = tree_splits[0].split("-")
        list_2 = list_1[0].split(",") + list_1[1].split(",")
        tree_labeled_names = set(list_2)
        assert(self_labeled_names == tree_labeled_names)
        RF_dist = len(self_splits - tree_splits) + len(tree_splits - self_splits)
        
    def IdentifyRepeatedSitePatterns(self):
        sequences = self.sequences
        seqNames = sequences.keys()
        firstSiteOfPatternToSitesWherePatternRepeats = {}
        distinctPatternsToFirstSiteWherePatternIsDetected = {}        
        for site in range(self.sequenceLength):
            sitePattern = ""
            for seqName in seqNames:
                sitePattern += sequences[seqName][site]
            if sitePattern in distinctPatternsToFirstSiteWherePatternIsDetected.keys():
                firstSite = distinctPatternsToFirstSiteWherePatternIsDetected[sitePattern]
                firstSiteOfPatternToSitesWherePatternRepeats[firstSite].append(site)
            else:
                distinctPatternsToFirstSiteWherePatternIsDetected[sitePattern] = site
                firstSiteOfPatternToSitesWherePatternRepeats[site] = []
        for distinctPattern in distinctPatternsToFirstSiteWherePatternIsDetected.keys():
            firstSite = distinctPatternsToFirstSiteWherePatternIsDetected[distinctPattern]
            self.firstSiteOfEachDistinctPattern.append(firstSite)
            self.siteArrangedBlockwiseForEachDistinctPattern.extend([firstSite]+firstSiteOfPatternToSitesWherePatternRepeats[firstSite])
            self.patternRepeats.append(len(firstSiteOfPatternToSitesWherePatternRepeats[firstSite])+1)
        return None
                        
    def ContainsSiblings(self,u_name,v_name):
        u = self.vertices[u_name]
        v = self.vertices[v_name]
        if u.degree==1 and v.degree==1: 
            if u.neighbors[0]==v.neighbors[0]:
                return True
        return False    
    
    def GetAllSplitsInTree(self,trivial=False):
            
        def GetSplitLabel(side,leafNames):
            side_1 = ','.join(sorted(side))
            side_2 = ','.join(sorted(list(set(leafNames)-set(side))))
            if side_1 < side_2:
                splitLabel = side_1+'-'+side_2
            else:
                splitLabel = side_2+'-'+side_1
            return splitLabel
        
        leaves = self.GetLeaves()
        leafNames = [l.name for l in leaves]
        smallerSplitSide={v:[v.name] for v in leaves}
        activeVertices=leaves
        verticesVisited = set(leaves)
        numberOfVerticesNotVisited = 2*len(leaves)-2
        while numberOfVerticesNotVisited>2:
            u = activeVertices[0]
            del activeVertices[0]
            numberOfVerticesNotVisited-=1
            v = list(set(u.neighbors)-verticesVisited)[0]
            verticesThatAreTwoEdgesAwayFromV= list(set(u.neighbors)-set([v]))
            for u_farAway in verticesThatAreTwoEdgesAwayFromV:
                verticesVisited.remove(u_farAway)
            if v in smallerSplitSide:
                smallerSplitSide[v]=smallerSplitSide[v]+smallerSplitSide[u]
                activeVertices.append(v)
                verticesVisited.add(v)
            else:
                smallerSplitSide[v]=smallerSplitSide[u]                 
            if len(smallerSplitSide[u])==1:
                if not trivial:
                    del smallerSplitSide[u]
        splits = map(lambda side: GetSplitLabel(side,leafNames), smallerSplitSide.values())
        return set(splits)        
    
    def GetRootedTree(self,u_name,v_name,t_u,t_v):
        RT = RootedTree()
        unvisitedVertices = self.GetLeaves()
        timesVisited = {v : 0 for v in self.vertices.values() if v.degree > 1}
        if self.GetVertex(v_name) in unvisitedVertices:
            unvisitedVertices.remove(self.GetVertex(v_name))
        elif self.GetVertex(u_name) in unvisitedVertices:
            unvisitedVertices.remove(self.GetVertex(u_name))
        while len(unvisitedVertices) > 0 :
            x = unvisitedVertices.pop()
            neighborsOfx = x.neighbors
            for y in neighborsOfx:
                if RT.ContainsEdge(x.name,y.name):
                    pass
                else:
                    t = self.edgeLengths[tuple(sorted((x.name,y.name)))]
                    RT.AddDirectedEdge(y.name, x.name, t)
                    timesVisited[y] += 1
                    if y.degree - timesVisited[y] == 1 and y.name not in [u_name,v_name]: 
                        unvisitedVertices.append(y)                                                                                    
        RT.AddDirectedEdge("h_root", u_name, t_u)
        RT.AddDirectedEdge("h_root", v_name, t_v)
        RT.SetRootByVertexName("h_root")
        if len(self.sequences) > 0 :
            RT.sequenceLength = self.sequenceLength
            RT.sequences = self.sequences.copy()
            RT.firstSiteOfEachDistinctPattern = self.firstSiteOfEachDistinctPattern[:]
            RT.patternRepeats = self.patternRepeats[:]
            RT.siteArrangedBlockwiseForEachDistinctPattern = self.siteArrangedBlockwiseForEachDistinctPattern[:]
        return RT
    
    def GetPathLengthsForEachLeafPair(self):
        u_name, v_name = self.edgeLengths.keys()[0]
        t_uv = self.edgeLengths[(u_name, v_name)]   
        RT = self.GetRootedTree(u_name, v_name, t_uv/2.0, t_uv/2.0)
        pathLengths = RT.ComputePathLengthsForAllLeafPairs()
        return pathLengths

    def GetPathLengthsForEachVertexPair(self):
        u_name, v_name = self.edgeLengths.keys()[0]
        t_uv = self.edgeLengths[(u_name, v_name)] 
        RT = self.GetRootedTree(u_name, v_name, t_uv/2.0, t_uv/2.0)
        pathLengths = RT.ComputePathLengthsForAllVertexPairs()
        return pathLengths
    
    def GetGraphDistance(self,vertexA_name,vertexB_name):
        vertexA = self.GetVertex(vertexA_name)
        vertexB = self.GetVertex(vertexB_name)
        if vertexA_name == vertexB_name:
            return 0
        verticesVisited=set([vertexA])
        verticesInCurrentLevel=[vertexA]
        verticesInNextLevel=[]
        continueSearch=True
        levelsSearched=1
        while continueSearch and len(verticesInCurrentLevel)>0:
            for u in verticesInCurrentLevel:
                unvisitedNeighborsOfu = list(set(u.neighbors)-verticesVisited)
                for v in unvisitedNeighborsOfu:
                    if v.name==vertexB.name:            
                        continueSearch=False
                        return levelsSearched
                    if v.degree > 1:
                        verticesInNextLevel.append(v)
                    verticesVisited.add(v)
            levelsSearched+=1
            verticesInCurrentLevel=verticesInNextLevel[:]
            verticesInNextLevel=[]
        return levelsSearched
    
    def GrowSubtreesUptoSizeK(self,k,subtreesSmallerThanK={}):
        if len(subtreesSmallerThanK)==0:
            leaves = self.GetLeaves()
            for v in leaves:
                subtreesSmallerThanK[v]=[v.name]
        verticesInCurrentLevel=[]
        for v in subtreesSmallerThanK.keys():
            if len(set(map(lambda l:l.name, v.neighbors))-set(subtreesSmallerThanK[v]))==1:
                verticesInCurrentLevel.append(v)
        verticesInCurrentLevel = sorted(verticesInCurrentLevel,key=lambda v:len(subtreesSmallerThanK[v]))    
        subtreesJustLargerThanK = {}
        while (len(verticesInCurrentLevel)>0):
            u = verticesInCurrentLevel[0]
            del verticesInCurrentLevel[0]
            unvisitedNeighborList_named = list(set(map(lambda l:l.name, u.neighbors))-set(subtreesSmallerThanK[u]))
            unvisitedNeighborList = map(lambda v_name: self.GetVertex(v_name),unvisitedNeighborList_named)
            if len(unvisitedNeighborList)==0:
                return subtreesSmallerThanK, {}
            else:
                v_name = list(set(map(lambda l:l.name, u.neighbors))-set(subtreesSmallerThanK[u]))[0]
                v = self.GetVertex(v_name)
            if v in subtreesSmallerThanK.keys():
                subtreesSmallerThanK[v].extend(subtreesSmallerThanK[u])
            else:
                subtreesSmallerThanK[v]=[v.name]+subtreesSmallerThanK[u]
            del subtreesSmallerThanK[u]
            if len(set(map(lambda l:l.name, v.neighbors))-set(subtreesSmallerThanK[v]))==1:
                if len(subtreesSmallerThanK[v]) <= k:
                    verticesInCurrentLevel.append(v)            
        for root in subtreesSmallerThanK.keys():
            if len(subtreesSmallerThanK[root]) > k and len(set(map(lambda l:l.name, root.neighbors))-set(subtreesSmallerThanK[root]))==1:
                subtreesJustLargerThanK[root] = subtreesSmallerThanK[root]
                del subtreesSmallerThanK[root]
        return subtreesJustLargerThanK, subtreesSmallerThanK
    
    def CanSubtreesBeExtractedForGivenTreeDepth(self,treeDepth):
        MSTSubtreesJustLargerThanK, _ = self.GrowSubtreesUptoSizeK(treeDepth, {})
        if len(MSTSubtreesJustLargerThanK) == 0:
            return False
        elif self.vertexCount == len(MSTSubtreesJustLargerThanK.values()[0]):
            return False
        else:
            return True
    
    def GetListOfBigVertexSetAndNestedVertexSets(self,treeDepth):
#         treeDepth = self.treeDepth
        vertexSetsList = []    
        MSTSubtreesJustLargerThanK, _ = self.GrowSubtreesUptoSizeK(treeDepth, {})
        for rootInMST in MSTSubtreesJustLargerThanK.keys():
            MSTSubtreesToGrow = {}
            for v_name in MSTSubtreesJustLargerThanK[rootInMST]:
                v = self.GetVertex(v_name)
                if v.degree == 1:
                    MSTSubtreesToGrow[v] = [v.name]            
            d = 2
            smallSubtreeSize = treeDepth/float(d)
            MSTSubtreesJustLargerThanKByf, _ = self.GrowSubtreesUptoSizeK(smallSubtreeSize, MSTSubtreesToGrow.copy())    
            if len(MSTSubtreesJustLargerThanKByf.values()) < 1:
                print ("Unable to find nested list")
#                 self.WriteToFile(tempPath+"unableToFindNestedList")
            while len(MSTSubtreesJustLargerThanK[rootInMST]) == len(MSTSubtreesJustLargerThanKByf.values()[0]) and smallSubtreeSize > 4:
                d += 1
                smallSubtreeSize = max(4,treeDepth/float(d))
                MSTSubtreesJustLargerThanKByf, _ = self.GrowSubtreesUptoSizeK(smallSubtreeSize, MSTSubtreesToGrow.copy())
            if len(MSTSubtreesJustLargerThanK[rootInMST]) > len(MSTSubtreesJustLargerThanKByf.values()[0]):
                vertexSetsList.append([MSTSubtreesJustLargerThanK[rootInMST],MSTSubtreesJustLargerThanKByf])
        return vertexSetsList
    
    def GetMLRootedTreeForGMM(self,pool):
        RTList = []
        for (u_name, v_name) in self.edgeLengths.keys():
            rt = self.GetRootedTree(u_name, v_name, 1, 1)
            RTList.append(rt)
        if pool == None:
            RTList = map(EMForMapForGMM,RTList)
        else:
            RTList = pool.map(EMForMapForGMM,RTList)
        maxLogLik = max([rt.logLik for rt in RTList])
        
        MLRTs = []    
        for rt in RTList:
            if rt.logLik == maxLogLik:
                MLRTs.append(rt)
        
        MLRT = choice(MLRTs,size=1)[0]
        return MLRT
    
    def WriteToFile(self,fileName):
        edgeListFile = open(fileName,'w')
        for u_name, v_name in self.edgeLengths.keys():
            edgeListFile.write(u_name+'\t'+v_name+'\t'+str(self.edgeLengths[tuple(sorted((u_name,v_name)))])+'\n')
        edgeListFile.close()

        
class VertexInRootedGraph:
    
    def __init__(self,name):
        self.name = name
        self.inDegree = 0
        self.outDegree = 0
        self.numberOfDescendants = 0
        self.parent = None
        self.children = []
        
    def GetParent(self):
        return self.parent
    
    def AddParent(self,parent):
        self.parent = parent
        self.inDegree = 1
        
    def RemoveParent(self):
        self.parent = None
        self.inDegree = 0
        
    def AddChild(self,child):
        self.children.append(child)
        self.outDegree += 1
        
    def RemoveChild(self,child):
        self.children.remove(child)
        self.outDegree -= 1
        
    def GetChildren(self):
        return self.children
    
    def GetDescendantNames(self):
        children = self.children[:]
        descendants = []
        while len(children) > 0:
            c = children.pop()
            if c.outDegree > 0:
                children.extend(c.children)
            else:
                descendants.append(c.name)
        return descendants

    def GetLabeledDescendantNames(self):
        children = self.children[:]
        labeledDescendantNames = []
        if not self.name.startswith("h_"):
                labeledDescendantNames.append(self.name)
        while len(children) > 0:
            c = children.pop()
            if not c.name.startswith("h_"):
                labeledDescendantNames.append(c.name)
            if c.outDegree > 0:
                children.extend(c.children)            
        return labeledDescendantNames

    def GetAllDescendantNames(self):
        children = self.children[:]
        descendants = []
        while len(children) > 0:
            c = children.pop()
            if c.outDegree > 0:
                children.extend(c.children)
                descendants.append(c.name)
            else:
                descendants.append(c.name)
        return descendants


class RootedTree():
    def __init__(self):
        
        self.root = None
        self.vertices = {}
        self.vertexCount = 0        
        self.edgeLengths = {}
        self.sequenceLength = 0
        self.sequences = {}        
        self.rootProb = array([0.0]*4)
        self.logLik = 0
        self.transitionMatrices = {}
        self.rateMatrices = {}
        self.location = {}
        self.time = {}
        self.host = {}        
        
    def AddVertex(self,u_name):
        if not u_name in self.vertices.keys(): 
            self.vertices[u_name] = VertexInRootedGraph(u_name)
            self.vertexCount += 1
          
    def GetVertex(self,u_name):
        return self.vertices[u_name]
    
    def ContainsVertex(self,u_name):        
        return u_name in self.vertices.keys()
    
    def RemoveVertex(self,u_name):
        u = self.vertices[u_name]
        for v in u.children:
            u.RemoveChild(v)
        u.RemoveParent()
        del self.vertices[u_name]
        self.vertexCount -= 1
        if u_name in self.sequences.keys():
            del self.sequences[u_name]    
        
    def AddTransitionMatrices(self,transitionMatricesToAdd):
        self.transitionMatrices.update(transitionMatricesToAdd)                    
    
    def RenameNonLeafVertices_defunct(self, h_ind):
        preOrderVertices = self.GetPreOrderTraversalWithoutLeaves()
        oldNameToNewName = {}
        for p in preOrderVertices:
            h_ind += 1
            oldNameToNewName[p.name] = "h_"+str(h_ind)                            
        for v_name in oldNameToNewName.keys():
            v = self.GetVertex(v_name)
            v.name = oldNameToNewName[v_name]
            self.vertices[v.name] = v        
        oldEdgesNames = self.edgeLengths.keys()
        for p_name, c_name in oldEdgesNames:
            t = self.GetEdgeLength(p_name, c_name)
            self.RemoveDirectedEdge(p_name, c_name)            
            if c_name in oldNameToNewName:
                self.AddDirectedEdge(oldNameToNewName[p_name], oldNameToNewName[c_name], t)
            else:
                self.AddDirectedEdge(oldNameToNewName[p_name], c_name, t)            
        for v_name in oldNameToNewName.keys():
            del self.vertices[v_name]        
        # rename sequences dictionary
        for v_name in oldNameToNewName.keys():
            if v_name in self.sequences.keys():
                self.sequences[oldNameToNewName[v_name]] = self.sequences[v_name] 
                del self.sequences[v_name]
                
    def GetDescendants(self,v_name):
        v = self.GetVertex(v_name)
        return v.GetDescendantNames()

    def GetAllDescendants(self,v_name):
        v = self.GetVertex(v_name)
        return v.GetAllDescendantNames()
        
    def GetAllNonTrivialClusters(self):
        clusters = set([])
        orderedVertices = self.GetPreOrderTraversalWithoutLeaves()[1:]
        for v in orderedVertices:
            clusters.update([','.join(sorted(self.GetDescendants(v.name)))])
        return clusters
    
    def GetAllNonRootClusters(self):
        clusters = set([])        
        orderedVertices = self.GetPostOrderTraversalWithoutRoot()
        for v in orderedVertices:
            clusters.update([','.join(sorted(v.GetLabeledDescendantNames()))])
        return clusters    
    
    def GetAllSplits(self):
        clusters = self.GetAllNonRootClusters()
        labeledVertexNames = set(self.GetLabeledVertexNames())        
        splits = set([])
        for cluster in clusters:
            side = set(cluster.split(','))
            if len(side) < len(labeledVertexNames) -1:        
                otherSide =  labeledVertexNames - side
                side_a = ','.join(sorted(list(side)))
                side_b = ','.join(sorted(list(otherSide)))
                if side_a < side_b:
                    splits.update([side_a+'-'+side_b])
                else:
                    splits.update([side_b+'-'+side_a])        
        return splits
    
    def GetAllNonTrivialSplits(self):
        clusters = self.GetAllNonTrivialClusters()
        leafNames = set(self.GetLeafNames())
        splits = set([])
        for cluster in clusters:
            side = set(cluster.split(','))
            if len(side) < len(leafNames) -1:        
                otherSide =  leafNames - side
                side_a = ','.join(sorted(list(side)))
                side_b = ','.join(sorted(list(otherSide)))
                if side_a < side_b:
                    splits.update([side_a+'-'+side_b])
                else:
                    splits.update([side_b+'-'+side_a])        
        return splits

    
    def SetRootByVertexName(self,rootName):
        if rootName not in self.vertices:
            self.AddVertex(rootName)
            self.vertexCount += 1
        self.root = self.vertices[rootName]
        
    def SetRoot(self):
        if self.root != None:
            pass
        else:
            for v in self.vertices.values():
                if v.inDegree == 0:
                    self.root = v
                    return v
            if self.root == None:
                print ("There is no vertex with in-degree 0")
                
    def GetRoot(self):
        if self.root != []:
            return self.root
        else:
            self.SetRoot()
            return self.root
           
    def GetUnrootedTree(self):
        T = Tree()
        for u_name, v_name in self.edgeLengths.keys():
            w = self.edgeLengths[u_name,v_name]
            T.AddEdge(u_name, v_name, w)
        rootToRemove = T.GetVertex(self.root.name)
        if rootToRemove.degree == 2:
            u, v = rootToRemove.neighbors
            t_uv = T.GetEdgeLength(u.name, rootToRemove.name) + T.GetEdgeLength(v.name, rootToRemove.name)
            T.AddEdge(u.name, v.name, t_uv)
        else:
            print (rootToRemove.degree)
            neighbors = rootToRemove.neighbors
            n_sel = choice(neighbors,size=1)[0]
            t_sel_root = T.GetEdgeLength(n_sel.name, rootToRemove.name)
            for n in neighbors:
                if n.name != n_sel.name:
                    t_n_root = T.GetEdgeLength(n_sel.name, rootToRemove.name)
                    t_sel_n = t_n_root + t_sel_root
                    T.AddEdge(n.name, n_sel.name, t_sel_n)
                                
        T.RemoveVertex(rootToRemove.name)
        if len(self.sequences) > 0 :
            T.sequences = {v.name : self.sequences[v.name] for v in T.vertices.values()}
            T.sequenceLength = self.sequenceLength
            T.firstSiteOfEachDistinctPattern = self.firstSiteOfEachDistinctPattern[:]
            T.patternRepeats = self.patternRepeats[:]
            T.siteArrangedBlockwiseForEachDistinctPattern = self.siteArrangedBlockwiseForEachDistinctPattern[:]
        
        return T
    
    def IdentifyRepeatedSitePatterns(self):
        sequences = self.sequences
        seqNames = sequences.keys()
        firstSiteOfPatternToSitesWherePatternRepeats = {}
        distinctPatternsToFirstSiteWherePatternIsDetected = {}
        self.siteArrangedBlockwiseForEachDistinctPattern = []
        self.firstSiteOfEachDistinctPattern = []
        self.patternRepeats = []
        for site in range(self.sequenceLength):
            sitePattern = ""
            for seqName in seqNames:
                sitePattern += sequences[seqName][site]
            if sitePattern in distinctPatternsToFirstSiteWherePatternIsDetected.keys():
                firstSite = distinctPatternsToFirstSiteWherePatternIsDetected[sitePattern]
                firstSiteOfPatternToSitesWherePatternRepeats[firstSite].append(site)
            else:
                distinctPatternsToFirstSiteWherePatternIsDetected[sitePattern] = site
                firstSiteOfPatternToSitesWherePatternRepeats[site] = []
        for distinctPattern in distinctPatternsToFirstSiteWherePatternIsDetected.keys():
            firstSite = distinctPatternsToFirstSiteWherePatternIsDetected[distinctPattern]
            self.firstSiteOfEachDistinctPattern.append(firstSite)
            self.siteArrangedBlockwiseForEachDistinctPattern.extend([firstSite]+firstSiteOfPatternToSitesWherePatternRepeats[firstSite])
            self.patternRepeats.append(len(firstSiteOfPatternToSitesWherePatternRepeats[firstSite])+1)
        return None
    
    def AddLeafSequences(self,sequences):
        if len(sequences) > 0:
            self.sequences.update(sequences)
            self.sequenceLength = len(sequences.values()[0])
            self.IdentifyRepeatedSitePatterns()
            
    def AddSequences(self,sequences):
        if len(sequences) > 0:
            self.sequences.update(sequences)
            if self.sequenceLength == 0:
                self.sequenceLength = len(sequences.values()[0])
                 
    def OrderCharsWrtSitePatternRepeatBlocks(self,sequence):    
        char_pos_tuple = [(sequence[i], self.siteArrangedBlockwiseForEachDistinctPattern[i]) for i in range(self.sequenceLength)]            
        char_pos_tuple = sorted(char_pos_tuple, key=itemgetter(1))
        return "".join(char_pos[0] for char_pos in char_pos_tuple)
        
    def GetSequence(self,v_name):
        return self.sequences[v_name]
    
    def hamDist_pattern(self,seq_u,seq_v):
        self.firstSiteOfEachDistinctPattern
        self.patternRepeats
        dist = 0.0
        for siteInd in range(len(self.firstSiteOfEachDistinctPattern)):
            site = self.firstSiteOfEachDistinctPattern[siteInd]
            numberOfRepeats = self.patternRepeats[siteInd]
            if seq_u[site] != seq_v[site]:
                dist += numberOfRepeats
        dist /= float(sum(self.patternRepeats))
        return dist           
        
    def AddDirectedEdge(self,parent_name,child_name,length=-1):
        if parent_name not in self.vertices:
            self.AddVertex(parent_name)
        if child_name not in self.vertices:
            self.AddVertex(child_name)
        parent=self.GetVertex(parent_name)
        child=self.GetVertex(child_name)
        parent.AddChild(child)
        child.AddParent(parent)
        self.transitionMatrices[parent_name,child_name] = array([[0.0]*4]*4)
        if length != -1:
            self.edgeLengths[parent_name,child_name] = length
            
    def AddDirectedEdges(self,edgeLengthsDic):
        for (u_name, v_name), t in edgeLengthsDic.items():
            self.AddDirectedEdge(u_name, v_name, t) 
    
       
    def AddEdgeLength(self,parent_name,child_name,length):
        self.edgeLengths[parent_name,child_name] = length
        
    def UpdateEdgeLength(self,parent_name,child_name,length):
        self.edgeLengths[parent_name,child_name] = length
        
    def GetEdgeLength(self,parent_name,child_name):
        return self.edgeLengths[parent_name,child_name]
    
    def GetEdgeLengths(self):
        return self.edgeLengths
    
    def ContainsEdge(self,parent_name,child_name):
        return (parent_name,child_name) in self.edgeLengths
    
    def RemoveDirectedEdge(self,parent_name,child_name):
        parent=self.GetVertex(parent_name)
        child=self.GetVertex(child_name)
        parent.RemoveChild(child)
        child.RemoveParent()
        del self.edgeLengths[parent_name,child_name]
        
    def ContractAShortestEdgeIncidentToTheRoot(self):
        root = self.root        
        c_l, c_r = root.children        
        t_l = self.GetEdgeLength(root.name, c_l.name)
        t_r = self.GetEdgeLength(root.name, c_r.name)
        self.RemoveDirectedEdge(root.name, c_l.name)
        self.RemoveDirectedEdge(root.name, c_r.name)        
        if t_l < t_r:
            self.AddDirectedEdge(c_l.name, c_r.name, t_l + t_r)            
        else:
            self.AddDirectedEdge(c_r.name, c_l.name, t_l + t_r)
        
        self.RemoveVertex(root.name)        
        return None

    def ContractAShortestEdgeIncidentToTheRootForGMM(self):
        root = self.root        
        c_l, c_r = root.children        
        t_l = self.GetEdgeLength(root.name, c_l.name)
        t_r = self.GetEdgeLength(root.name, c_r.name)
        self.RemoveDirectedEdge(root.name, c_l.name)
        
        del self.transitionMatrices[(root.name,c_l.name)]
        self.RemoveDirectedEdge(root.name, c_r.name)        
        del self.transitionMatrices[(root.name,c_r.name)]        
        
        if t_l < t_r:
            parent_name = c_l.name
            child_name = c_r.name            
        else:
            parent_name = c_r.name
            child_name = c_l.name
            
        seq_parent = self.sequences[parent_name]
        seq_child = self.sequences[child_name]
        P = GetMaxLikEstimateOfTransitionMatrix(seq_parent, seq_child)
        t = hamDist(seq_parent, seq_child)
        self.AddDirectedEdge(parent_name, child_name, t)
        self.transitionMatrices[(parent_name, child_name)] = P
        self.RemoveVertex(root.name)        
        return None

    def GetDistance(self,u_name,v_name):
        label_u = self.vertexLabels[u_name]
        label_v = self.vertexLabels[v_name]
        d = 0
        for i in range(len(label_u)):
            if label_u[i] != label_v[i]:
                d += 1
        return d
    
    def GetLabeledVertexNames(self):
        labeledVertexNames = []
        for v_name in self.vertices.keys():
            if v_name.startswith("h_"):
                pass
            else:
                labeledVertexNames.append(v_name)
        return labeledVertexNames
    
    def GetLeafNames(self):
        leafNames = []
        for v in self.vertices.values():
            if v.outDegree + v.inDegree == 1:
                leafNames.append(v.name)
        return leafNames
    
    def GetLeaves(self):
        leaves = []
        for v in self.vertices.values():
            if v.outDegree == 0:
                leaves.append(v)
        return leaves
    
    def GetNonLeafVertices(self):
        nonLeaves=[]
        for v in self.vertices.values():
            if v.degree != 1:
                nonLeaves.append(v)
        return nonLeaves
    
    def GetLeavesAndNonLeafVertices(self):
        leaves = []
        nonLeaves = []
        for v in self.vertices.values():
            if v.outDegree == 0:
                leaves.append(v)
            else:
                nonLeaves.append(v)
        return leaves, nonLeaves
    
    def ComputePathLengthsForAllLeafPairs(self):
        root = self.GetRoot()
        rootName = root.name
        edgeWeights = self.edgeLengths
        leaves = self.GetLeaves()
        distances = {}
        def ComputePathLength(path):
            l = 0.0
            for i in range(len(path)-1):
                l += edgeWeights[path[i+1],path[i]]
            return l
                    
        def ComputePathLengthForLeafPair(leaf_1,leaf_2):
            lcaNotFound = True
            ancestorsOfLeaf_1 = [leaf_1.name]
            ancestorsOfLeaf_2 = [leaf_2.name]
            child_1 = leaf_1
            child_2 = leaf_2
            while lcaNotFound:                
                if child_1.name!=rootName:    
                    parent_1 = child_1.parent
                    ancestorsOfLeaf_1.append(parent_1.name)
    
                if child_2.name!=rootName:    
                    parent_2 = child_2.parent
                    ancestorsOfLeaf_2.append(parent_2.name)
                
                if parent_1.name in ancestorsOfLeaf_2:
                    lcaNotFound = False
                    l_1 = ComputePathLength(ancestorsOfLeaf_1)
                    if parent_1.name == rootName or parent_1.name == parent_2.name:
                        l_2 = ComputePathLength(ancestorsOfLeaf_2)
                    else:
                        l_2 = ComputePathLength(ancestorsOfLeaf_2[:ancestorsOfLeaf_2.index(parent_1.name)+1])                    
                if parent_2.name in ancestorsOfLeaf_1:
                    lcaNotFound = False
                    l_2 = ComputePathLength(ancestorsOfLeaf_2)
                    if parent_2.name == rootName:
                        l_1 = ComputePathLength(ancestorsOfLeaf_1)
                    else:
                        l_1 = ComputePathLength(ancestorsOfLeaf_1[:ancestorsOfLeaf_1.index(parent_2.name)+1])
                child_1 = parent_1
                child_2 = parent_2
            return (l_1 + l_2)            
                
        for leaf_1 in leaves:
            for leaf_2 in leaves[(leaves.index(leaf_1)+1):]:
                distances[tuple(sorted([leaf_1.name,leaf_2.name]))] = ComputePathLengthForLeafPair(leaf_1,leaf_2)    
        
        return distances
    
    def ComputePathLengthsForAllVertexPairs(self, includeRoot = False):
        root = self.GetRoot()
        rootName = root.name
        edgeWeights = self.edgeLengths
        vertices = self.GetPostOrderTraversalWithoutRoot()
        if includeRoot:
            vertices.append(root)
        distances = {}
        def ComputePathLength(path):
            l = 0.0
            for i in range(len(path)-1):
                l += edgeWeights[path[i+1],path[i]]
            return l
                    
        def ComputePathLengthForVertexPair(v_1,v_2):
            lcaNotFound = True
            ancestorsOfLeaf_1 = [v_1.name]
            ancestorsOfLeaf_2 = [v_2.name]
            child_1 = v_1
            child_2 = v_2
            while lcaNotFound:                
                if child_1.name!=rootName:    
                    parent_1 = child_1.parent
                    ancestorsOfLeaf_1.append(parent_1.name)
    
                if child_2.name!=rootName:    
                    parent_2 = child_2.parent
                    ancestorsOfLeaf_2.append(parent_2.name)
                
                if parent_1.name in ancestorsOfLeaf_2:
                    lcaNotFound = False
                    l_1 = ComputePathLength(ancestorsOfLeaf_1)
                    if parent_1.name == rootName or parent_1.name == parent_2.name:
                        l_2 = ComputePathLength(ancestorsOfLeaf_2)
                    else:
                        l_2 = ComputePathLength(ancestorsOfLeaf_2[:ancestorsOfLeaf_2.index(parent_1.name)+1])                    
                try:
                    parent_2.name in ancestorsOfLeaf_1
                except:
                    print (v_1.name, v_2.name)
                
                if parent_2.name in ancestorsOfLeaf_1:
                    lcaNotFound = False
                    l_2 = ComputePathLength(ancestorsOfLeaf_2)
                    if parent_2.name == rootName:
                        l_1 = ComputePathLength(ancestorsOfLeaf_1)
                    else:
                        l_1 = ComputePathLength(ancestorsOfLeaf_1[:ancestorsOfLeaf_1.index(parent_2.name)+1])
                child_1 = parent_1
                child_2 = parent_2
            return (l_1 + l_2)            
                
        for v_1 in vertices:
            for v_2 in vertices[(vertices.index(v_1)+1):]:
                distances[tuple(sorted([v_1.name,v_2.name]))] = ComputePathLengthForVertexPair(v_1,v_2)    
        
        return distances

# Tree traversal operations
    def GetPostOrderTraversalWithoutRoot(self):
        unvisitedVertices = self.GetLeaves()
        orderedVertices = []
        timesVisited = {v : 0 for v in self.vertices.values() if v.outDegree > 0}
        while len(unvisitedVertices) > 0:
            c = unvisitedVertices.pop()
            orderedVertices.append(c)
            p = c.parent
            if p.inDegree == 1:
                timesVisited[p] += 1
                if p.outDegree - timesVisited[p] == 0:           
                    unvisitedVertices.append(p)
        return orderedVertices
    
    def GetPreOrderTraversalWithoutLeaves(self):
        root = self.GetRoot()        
        orderedVertices = [root]
        parents = root.children[:]
        while len(parents) > 0:
            p = parents.pop()
            if len(p.children) > 0:
                orderedVertices.append(p)
                parents.extend(p.children)
        return orderedVertices

# Initialization   
    
    def InitializeNonLeafSequencesUsingMPHartigan(self):
        allSequences = self.sequences            
        leaves = self.GetLeaves()
        root = self.GetRoot()
        postOrderVerticesWithRoot = self.GetPostOrderTraversalWithoutRoot()
        postOrderVerticesWithRoot.append(root)
        preOrderNonLeafVertices = self.GetPreOrderTraversalWithoutLeaves()
        nonLeafSequences = {v.name : "" for v in preOrderNonLeafVertices}
        preOrderNonLeafVerticesWithoutRoot = preOrderNonLeafVertices[1:]    
        def GetStatesForSite(site):
            V = {v : set([allSequences[v.name][site]]) for v in leaves}
            VU = {v : set([allSequences[v.name][site]]) for v in leaves}        
            for p in postOrderVerticesWithRoot:
                if p.outDegree > 0:
                    VU[p]=set([])
                    statesInProgenyOfP = {}
                    for c in p.children:
                        stateInc = list(VU[c])[0]
                        if stateInc in statesInProgenyOfP:
                            statesInProgenyOfP[stateInc] += 1
                        else:
                            statesInProgenyOfP[stateInc] = 0
                    sorted_states = sorted(statesInProgenyOfP.items(), key=itemgetter(1), reverse=True)
                    maxCount = 0
                    for state, count in sorted_states:
                        if maxCount <= count:
                            maxCount = count
                            VU[p].update(set([state]))     
            V[root] = set(choice(list(VU[root]),size=1))
            nonLeafSequences[root.name] += list(V[root])[0]       
            for p in preOrderNonLeafVerticesWithoutRoot:
                if V[p.parent] in VU[p]:
                    V[p] = V[p.parent]
                else:
                    V[p] = set(choice(list(VU[p]),size=1))
                nonLeafSequences[p.name] += list(V[p])[0]
            return None
        
        map(GetStatesForSite,range(self.sequenceLength))
        self.sequences.update(nonLeafSequences)
        return None
    
    def InitializeNonLeafSequencesUsingMPUsingPatternSavings(self):
        numberOfDistinctPatterns = len(self.firstSiteOfEachDistinctPattern)
        allSequences = self.sequences
        leaves = self.GetLeaves()
        root = self.GetRoot()
        postOrderVerticesWithRoot = self.GetPostOrderTraversalWithoutRoot()
        postOrderVerticesWithRoot.append(root)
        preOrderNonLeafVertices = self.GetPreOrderTraversalWithoutLeaves()
        nonLeafSequences = {v.name : "" for v in preOrderNonLeafVertices}
        preOrderNonLeafVerticesWithoutRoot = preOrderNonLeafVertices[1:]
        def GetStatesForSiteOfDistinctPattern(siteInd):
            site = self.firstSiteOfEachDistinctPattern[siteInd]
            numberOfRepeats = self.patternRepeats[siteInd]
            V = {v : set([allSequences[v.name][site]]) for v in leaves}
            VU = {v : set([allSequences[v.name][site]]) for v in leaves}    
            for p in postOrderVerticesWithRoot:
                if p.outDegree > 0:
                    VU[p]=set([])
                    statesInProgenyOfP = {}
                    for c in p.children:
                        stateInc = list(VU[c])[0]
                        if stateInc in statesInProgenyOfP:
                            statesInProgenyOfP[stateInc] += 1
                        else:
                            statesInProgenyOfP[stateInc] = 0
                    sorted_states = sorted(statesInProgenyOfP.items(), key=itemgetter(1), reverse=True)
                    maxCount = 0
                    for state, count in sorted_states:
                        if maxCount <= count:
                            maxCount = count
                            VU[p].update(set([state]))     
            V[root] = set(choice(list(VU[root]),size=1))
            nonLeafSequences[root.name] += "".join(list(V[root])*numberOfRepeats)      
            for p in preOrderNonLeafVerticesWithoutRoot:
                if V[p.parent] in VU[p]:
                    V[p] = V[p.parent]
                else:
                    V[p] = set(choice(list(VU[p]),size=1))
                nonLeafSequences[p.name] += "".join(list(V[p])*numberOfRepeats)
            return None
        
        map(GetStatesForSiteOfDistinctPattern,range(numberOfDistinctPatterns))       
        for v_name in nonLeafSequences.keys():
            nonLeafSequences[v_name] = self.OrderCharsWrtSitePatternRepeatBlocks(nonLeafSequences[v_name])        
        self.sequences.update(nonLeafSequences)
        return None
    
    def InitializeEdgeLengths(self):
        allSequences = self.sequences
        postOrderVerticesWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        for c in postOrderVerticesWithoutRoot:
            p = c.parent
            seq_c = allSequences[c.name]
            seq_p = allSequences[p.name]
            self.UpdateEdgeLength(p.name, c.name, max(self.hamDist_pattern(seq_c, seq_p),pow(10,-5)))
        return None

    def InitializeRootProbUsingPatternSavings(self):
        root = self.GetRoot()
        rootSeq = self.sequences[root.name]
        rootProb = [0.0] * 4
        for siteInd in range(len(self.firstSiteOfEachDistinctPattern)):
            rootProb[DNA.index(rootSeq[self.firstSiteOfEachDistinctPattern[siteInd]])] += self.patternRepeats[siteInd]  
        for i in range(4):
            rootProb[i] /= float(len(rootSeq))
        self.rootProb = rootProb
        return None
    
    def InitializeTransitionMatricesAndEdgeLengths(self):
        for (u_name, v_name) in self.transitionMatrices.keys():            
            seq_u = self.sequences[u_name]
            seq_v = self.sequences[v_name]
            P = GetMaxLikEstimateOfTransitionMatrix(seq_u, seq_v)
            self.transitionMatrices[(u_name,v_name)] = P
            t = hamDist(seq_u,seq_v)
            self.edgeLengths[(u_name,v_name)] = t
                
# Expectation    
    def ComputeExpectedAncestralStatesAndUpdateLogLikUsingPatternSavingForTransitionMatrices(self):
        ancestralStates = {}
        root = self.GetRoot()
        leaves = []
        nonLeafVertices = []
        for v in self.vertices.values():
            if v.outDegree == 0:
                leaves.append(v)
            else:
                nonLeafVertices.append(v)
                ancestralStates[v.name]=""
        priorProbs = {}
        priorProbs[root] = self.rootProb        
        allSequences = self.sequences
        postOrderVerticesWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        preOrderTraversalWithoutLeaves = self.GetPreOrderTraversalWithoutLeaves()
        self.logLik = 0
        rootProb = self.rootProb
        
        def ComputeExpectedStatesAndUpdateLogLikPerSite(siteInd):
            site = self.firstSiteOfEachDistinctPattern[siteInd]
            numberOfRepeats = self.patternRepeats[siteInd]
            conditionalLikelihoods = {}        
            for v in leaves:
                conditionalLikelihoods[v] = [0.0]*4
                conditionalLikelihoods[v][DNA.index(allSequences[v.name][site])] = 1.0 
            for v in nonLeafVertices:
                conditionalLikelihoods[v] = [1.0]*4
            for c in postOrderVerticesWithoutRoot:                              
                p = c.parent
                P = self.transitionMatrices[(p.name,c.name)]         
                for nuc_p in range(4):
                    partialLikelihood_c_nuc_p = 0.0
                    for nuc_c in range(4):
                        partialLikelihood_c_nuc_p += P[nuc_p,nuc_c]*conditionalLikelihoods[c][nuc_c]
                    conditionalLikelihoods[p][nuc_p] *= partialLikelihood_c_nuc_p                    
            
            lik = 0
            for nuc_root in range(4):
                lik += rootProb[nuc_root]*conditionalLikelihoods[root][nuc_root]            
            self.logLik += log(max(lik,pow(10,-10)))*numberOfRepeats
                
            for p in preOrderTraversalWithoutLeaves:                
                maxProb = -1
                stateWithMaxProb = "" 
                for nuc in range(4):
                    prob_nuc = priorProbs[p][nuc]*conditionalLikelihoods[p][nuc]
                    if prob_nuc > maxProb:
                        stateWithMaxProb = DNA[nuc]
                        maxProb = prob_nuc
                if stateWithMaxProb == "":
                    print ("Error")                
                ancestralStates[p.name] += stateWithMaxProb*numberOfRepeats
                for c in p.children:
                    if c.outDegree > 0:
                        P = self.transitionMatrices[(p.name,c.name)]
                        priorProbs[c] = P[DNA.index(stateWithMaxProb),:]                    
            return None
        
        map(ComputeExpectedStatesAndUpdateLogLikPerSite,range(len(self.firstSiteOfEachDistinctPattern)))
        for v_name in ancestralStates.keys():                        
            self.sequences[v_name] = self.OrderCharsWrtSitePatternRepeatBlocks(ancestralStates[v_name])
        return None

# Maximization
    def OptimizeRootProbability(self):
        root = self.GetRoot()
        rootSeq = self.sequences[root.name]
        if len(rootSeq)==0:
            print (len(rootSeq))
        rootProb = [0.0] * 4
        for siteInd in range(len(self.firstSiteOfEachDistinctPattern)):
            site = self.firstSiteOfEachDistinctPattern[siteInd]
            numberOfRepeats = self.patternRepeats[siteInd]
            rootProb[DNA.index(rootSeq[site])] += float(numberOfRepeats)    
        for i in range(4):
            rootProb[i] /= float(sum(self.patternRepeats))
        self.rootProb = rootProb
        return None
    
    def OptimizeTransitionMatrices(self):
        self.InitializeTransitionMatricesAndEdgeLengths()    
                
    def PerformEMForGMMUsingPatternSavings(self,maxIterations = 20, logLikImprovementThreshold = 1.0, verbose = False):
        self.InitializeNonLeafSequencesUsingMPUsingPatternSavings()
        # Initialized ancestral states using MP
        self.InitializeRootProbUsingPatternSavings()
        # "Initialized root probability"
        self.InitializeTransitionMatricesAndEdgeLengths()
        # "Initialized transition matrices"
        logLikConvergenceNotReached = True
        iteration = 1
        currentLogLik = -1*pow(10,10)   
        current_rootSeq = self.sequences[self.root.name]
        while logLikConvergenceNotReached and iteration <= maxIterations:        
            # Expectation step
            self.ComputeExpectedAncestralStatesAndUpdateLogLikUsingPatternSavingForTransitionMatrices()
            changeInRootSeq = hamDist(current_rootSeq, self.sequences[self.root.name])
            current_rootSeq = self.sequences[self.root.name]
            if verbose:
                print ("Expected states computed")            
#                 print changeInRootSeq
            # Maximization steps
            # Given root sequence, optimize rootProb            
            self.OptimizeRootProbability()
            if verbose:
                print ("Optimized root probability")                                   
            # Given t and rootProb, optimize Q
            self.OptimizeTransitionMatrices()
            if verbose:
                print ("Optimized transition matrices ")
                print ("iteration: "+str(iteration)+". Change in root seq is " + str(changeInRootSeq))                
            iteration += 1
            logLikConvergenceNotReached = self.logLik - currentLogLik > logLikImprovementThreshold            
            currentLogLik = self.logLik
        
    def ContractShortEdgesIncidentToVerticesWithIdenticalExpectedStates(self):
        orderedVertices = self.GetPreOrderTraversalWithoutLeaves()
        edgesToContract = []
        for u in orderedVertices:
            for c in u.children:
                if self.sequences[u.name] == self.sequences[c.name]:
                    edgesToContract.append([u.name,c.name])
                                                           
        for u_name, v_name in edgesToContract:
            t_uv = self.GetEdgeLength(u_name, v_name)
            v = self.GetVertex(v_name)
            v_children = v.children
            for c in v_children:
                t_vc = self.GetEdgeLength(v_name, c.name)
                self.AddDirectedEdge(u_name, c.name, t_uv + t_vc)
                self.RemoveDirectedEdge(v_name, c.name)
            self.RemoveVertex(v_name)
            
        return None    
            
    def ReturnLogLikForRootingAtEdge(self):        
        preOrderVerticesWithoutRoot = self.GetPreOrderTraversalWithoutLeaves()[1:]
        logLik = 0.0            
        for u in preOrderVerticesWithoutRoot:
            seq_u = self.sequences[u.name]
            for v in u.children:
                seq_v = self.sequences[v.name]
                if (u.name,v.name) in self.transitionMatrices:                
                    P = self.transitionMatrices[(u.name,v.name)]
                else:
                    P = self.transitionMatrices[(v.name,u.name)]
                try:
                    logLik += sum(map(lambda site: log(P[DNA.index(seq_u[site]),DNA.index(seq_v[site])]), range(self.sequenceLength)))
                except:
                    print (P)
        return logLik
    
    def FindEdgesForRootingForGMM(self,pool):
        T = Tree()
        T.AddEdges(self.edgeLengths)
        T.AddSequences(self.sequences)
        RTList = []
        edgeList = self.edgeLengths.keys()
        for edge in edgeList:
            t = T.GetEdgeLength(edge[0], edge[1])
            RT = T.GetRootedTree(edge[0], edge[1], t/2.0, t/2.0)
            RT.transitionMatrices = self.transitionMatrices.copy()
            RTList.append(RT)
        if pool == None:
            logLikList = map(poolFuncForComputingLogLikForRooting,RTList)
        else:
            logLikList = pool.map(poolFuncForComputingLogLikForRooting,RTList)
        tup_edge_logLik = zip(edgeList,logLikList)
        tup_edge_logLik.sort(key= itemgetter(1),reverse = True)
        maxLogLik = tup_edge_logLik[0][1]
        edgesToReturn = [tup_edge_logLik[0][0]]
        for edge, logLik in tup_edge_logLik[1:]:
            if logLik == maxLogLik:
                edgesToReturn.append(edge)           
        return edgesToReturn

    def GetNewickLabel(self):
        labelAtVertex = {v : "" for v in self.vertices.values() if v.outDegree > 0}
        timesVisited = {v : 0 for v in labelAtVertex.keys() if v.outDegree > 0}
        labelAtVertex.update({v : v.name for v in self.GetLeaves()})    
        postOrderTraversalWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        for c in postOrderTraversalWithoutRoot:            
            p = c.parent
            timesVisited[p] += 1
            if timesVisited[p] == 1:
                labelAtVertex[p] += "("+labelAtVertex[c]+":"+str(self.GetEdgeLength(p.name, c.name))
            else :
                labelAtVertex[p] += ","+labelAtVertex[c]+":"+str(self.GetEdgeLength(p.name, c.name))
            if p.outDegree - timesVisited[p] == 0:
                labelAtVertex[p] += ")"
        labelAtRoot = labelAtVertex[self.GetRoot()] + ";"
        return labelAtRoot

    def WriteToFile(self,fileName,fileFormat="edgeList"):
        treeFile = open(fileName,'w')
        if fileFormat == "edgeList":            
            for parent_name, child_name in self.edgeLengths.keys():
                treeFile.write(parent_name+'\t'+child_name+'\t'+ str('%1.9f'% self.edgeLengths[parent_name,child_name])+'\n')            
        elif fileFormat == "newick":
            treeFile.write(self.GetNewickLabel()+"\n")
        else:
            print ("file format not recognized")
        treeFile.close()

    def WriteTransitionMatricesToFile(self,fileName):
        transitionMatricesFile = open(fileName)
        for (p_name, c_name), P in self.transitionMatrices.items():
            transitionMatricesFile.write(p_name+'\t'+c_name)
            for row in range(4):
                for col in range(4):
                    transitionMatricesFile.write('\t'+P[row,col])
            transitionMatricesFile.write('\n')
        transitionMatricesFile.close()
        
    def ComputeLikelihoodForSite(self,site):                
        root = self.GetRoot()        
        leaves = []
        nonLeafVertices = []
        for v in self.vertices.values():
            if v.outDegree == 0:
                leaves.append(v)
            else:
                nonLeafVertices.append(v)                        
        allSequences = self.sequences
        postOrderVerticesWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        rootProb = self.rootProb 
        conditionalLikelihoods = {}        
        for v in leaves:
            conditionalLikelihoods[v] = [0.0]*4
            conditionalLikelihoods[v][DNA.index(allSequences[v.name][site])] = 1.0 
        for v in nonLeafVertices:
            conditionalLikelihoods[v] = [1.0]*4
        for c in postOrderVerticesWithoutRoot:                              
            p = c.parent
            P = self.transitionMatrices[(p.name,c.name)]         
            for nuc_p in range(4):
                partialLikelihood_c_nuc_p = 0.0
                for nuc_c in range(4):
                    partialLikelihood_c_nuc_p += P[nuc_p,nuc_c]*conditionalLikelihoods[c][nuc_c]
                conditionalLikelihoods[p][nuc_p] *= partialLikelihood_c_nuc_p                            
        lik = 0
        for nuc_root in range(4):
            lik += rootProb[nuc_root]*conditionalLikelihoods[root][nuc_root]            
        return lik
    
def FindRoot(u):
    if u.inDegree == 0:
        return u
    else:
        p = u.parent
        while p.inDegree != 0:
            p = p.parent
        return p    

def Find(u):
    if u.rep != u:
        u.rep = Find(u.rep)
    return u.rep

def Union(u,v):
    repOfu = Find(u)
    repOfv = Find(v)
    if repOfu == repOfv: 
        pass
    else:
        if repOfu.rank < repOfv.rank:
            repOfu.rep = repOfv 
        elif repOfu.rank > repOfv.rank:
            repOfv.rep = repOfu            
        else:
            repOfu.rep = repOfv            
            repOfv.rank+=1

def UnionForJGAA(u,v,unionCount,edgeCount):
    repOfu = Find(u)
    repOfv = Find(v)
    if repOfu == repOfv: 
        pass
    else:
        if repOfu.rank < repOfv.rank:
            repOfu.rep = repOfv 
            repToKeep = repOfv
            repToRemove = repOfu
        elif repOfu.rank > repOfv.rank:
            repOfv.rep = repOfu            
            repToKeep = repOfu
            repToRemove = repOfv
        else:
            repOfu.rep = repOfv            
            repToKeep = repOfv
            repToRemove = repOfu
            repOfv.rank+=1
        if repToKeep not in unionCount:            
            unionCount[repToKeep]=1
            edgeCount[repToKeep]=1
        else:
            unionCount[repToKeep]+=1
            edgeCount[repToKeep]+=1
        if repToRemove in unionCount:
            unionCount[repToKeep]+=unionCount[repToRemove]
            edgeCount[repToKeep]+=edgeCount[repToRemove]
            del unionCount[repToRemove]
            del edgeCount[repToRemove]
            
def GetSortedDistancesAndVertexPairList(distances):
    tuplesToSort=[]
    for (vertex_i, vertex_j),d in distances.iteritems():
        tuplesToSort.append([d,vertex_i,vertex_j])              
    tuplesToSort.sort(key=lambda x:x[0])
    sortedDistances = []
    vertexPairList = []
    for tup in tuplesToSort:
        d,vertex_i,vertex_j = tup
        sortedDistances.append(d)
        vertexPairList.append(sorted((vertex_i,vertex_j)))    
    return sortedDistances, vertexPairList

def ComputeCompNeighboursUnionCountAndEdgeCount(E_w):
    compNeighbors={}
    unionCount={}
    edgeCount={}
    for u,v,oldRep_u,oldRep_v in E_w:
        if u not in compNeighbors:
            compNeighbors[u]=set([oldRep_v])
        else:
            compNeighbors[u].update(set([oldRep_v]))
        if v not in compNeighbors:
            compNeighbors[v]=set([oldRep_u])
        else:
            compNeighbors[v].update(set([oldRep_u]))
        if Find(u)!= Find(v):
            UnionForJGAA(u, v, unionCount, edgeCount)
        else:
            edgeCount[Find(u)]+=1
    return compNeighbors, unionCount, edgeCount

def GetEdgesInMSTBasedOnVertexRank(structuredEdges,unionCount,edgeCount):
    edgesToReturn=[]
    edge_compA_compB_minRank_maxRank = []
    for (u,v,compA,compB) in structuredEdges:
        minRank = min(u.rank,v.rank)
        maxRank = max(u.rank,v.rank)
        compA.rep = compA
        compB.rep = compB
        edge_compA_compB_minRank_maxRank.append((u,v,compA,compB,minRank,maxRank))
        edge_compA_compB_minRank_maxRank = sorted(edge_compA_compB_minRank_maxRank,key=itemgetter(4,5))
        edge_compA_compB_minRank_maxRank = sorted(edge_compA_compB_minRank_maxRank,key=itemgetter(4,5))
    for (u,v,compA,compB,_,_) in edge_compA_compB_minRank_maxRank:
        if Find(compA)!=Find(compB):
            UnionForJGAA(compA,compB,unionCount,edgeCount)
            edgesToReturn.append([u,v])
    return edgesToReturn

def GetKeyBasedOnDeltaMaxAndVertexName(u):
    return (u.deltaMax,u.name)

def GetMLVRMSTFromDistances(distanceMatrix):
    M = Tree()
    sortedDistances, vertexPairList = GetSortedDistancesAndVertexPairList(distanceMatrix)
    E_w = []
    edgesInCycles = []
    w_old = sortedDistances[0]
    for edgeInd in range(len(sortedDistances)):
        w = sortedDistances[edgeInd]
        if w > w_old:
            compNeighbors, unionCount, edgeCount = ComputeCompNeighboursUnionCountAndEdgeCount(E_w)
            E_temp = []
            for u,v,comp_a,comp_b in E_w:
                if Find(u) in edgeCount and edgeCount[Find(u)] == unionCount[Find(u)]:
                    del unionCount[Find(u)]
                    del edgeCount[Find(u)]
                    M.AddEdge(u.name, v.name, w_old)
                else:
                    E_temp.append([u,v,comp_a,comp_b])
            edgesInCycles.append(E_temp)
            # Compute delta_max
            for u in compNeighbors:
                u.deltaMax+=len(compNeighbors[u])
            E_w=[]
            w_old = w
        u_name,v_name = vertexPairList[edgeInd]
        if not M.ContainsVertex(u_name):
            M.AddVertex(u_name)
        if not M.ContainsVertex(v_name):
            M.AddVertex(v_name)
        u = M.GetVertex(u_name)
        v = M.GetVertex(v_name)   
        if Find(u)!=Find(v):
            E_w.append((u,v,Find(u),Find(v)))
    compNeighbors, unionCount, edgeCount = ComputeCompNeighboursUnionCountAndEdgeCount(E_w)
    E_temp = []
    for u,v,comp_a,comp_b in E_w:
        if Find(u) in edgeCount and edgeCount[Find(u)] == unionCount[Find(u)]:
            del unionCount[Find(u)]
            del edgeCount[Find(u)]
            M.AddEdge(u.name, v.name, w)
        else:
            E_temp.append([u,v,comp_a,comp_b])
    edgesInCycles.append(E_temp)
    for u in compNeighbors:
        u.deltaMax+=len(compNeighbors[u])
    # Assign ranks based on deltaMax and vertex names
    vertexObjList = M.vertices.values()
    vertexObjList = sorted(vertexObjList,key=GetKeyBasedOnDeltaMaxAndVertexName)
    for i in range(len(vertexObjList)):
        vertexObjList[i].rank = i+1
    for structuredEdges in edgesInCycles:
        edgesToReport = GetEdgesInMSTBasedOnVertexRank(structuredEdges,unionCount,edgeCount)
        for u,v in edgesToReport:
            M.AddEdge(u.name, v.name, sortedDistances[vertexPairList.index(sorted((u.name, v.name)))])
    return M