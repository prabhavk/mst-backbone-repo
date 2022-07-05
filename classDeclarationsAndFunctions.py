from operator import itemgetter
from numpy.random import choice, sample
from distanceCalculator import hamDist, GetBaseFreq
from numpy import array
from math import log
import subprocess as sub
import cmath as cmath
from scipy.linalg import logm
from scipy.optimize._minimize import minimize_scalar
from scipy.optimize import minimize
import os 
from MarkovModels import GenerateProbabilityDistributionWithMinVal, GenerateQForStationaryDistribution, GenerateRandomQ, Get11FreeRates, NormalizeQMatrix,\
    ComputeQMatrixFromVectorOfOffDiagonalElements,\
    ComputeProbabilityMatrixUsingMatrixExponentiation,\
    ComputeQVectorOfOffDiagonalElementsFromQMatrix,\
    GetMaxLikEstimateOfTransitionMatrix, GetSubstitutionCountMatrix, GetStationaryDistribution
# from config import tempPath
from itertools import repeat
from numpy.random import uniform

# DNA=["A","C","G","T"]
DNA=["T","C","A","G"] # INDELIBLE INDEX

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

def poolFuncForComputingLogLikOfRootingAtV(Tv_tuple):
    T, v = Tv_tuple
    logLik = T.ComputeLogLikelihoodForRootingAtVertex(v.name)
    return logLik

def poolFuncForComputingLogLikOfRootingAlongEdge(Tedge_tuple):
    T, (u_name, v_name) = Tedge_tuple
    try:
        logLik = T.ComputeLogLikelihoodForRootingAlongEdge(u_name,v_name)
    except:
        print (u_name, v_name)
        return None
    return logLik

def EMForMapWithFixedStatesAndEdgeLengthsAndEstimatedQ(RT):
    RT.PerformEMWithFixedStatesAndEdgeLengthsAndEstimatedQ()
    return RT

# rep is short of representative
class Vertex:
    def __init__(self,name):
        self.rateMatrix = array([[0.0]*4]*4)
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
        self.vertices = {}
        self.vertexCount = 0
        self.edgeLengths = {}
        self.sequences = {}
        self.edgeLogLikelihoods = {}
        self.sequenceLength = 0
        self.indelible_replicates = 1
        self.siteArrangedBlockwiseForEachDistinctPattern = []
        self.firstSiteOfEachDistinctPattern = []
        self.patternRepeats = []
        self.treeDepth = 16
        self.QForEdge = {}        
        
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
        u.AddNeighbor(v)
        v.AddNeighbor(u)
        self.edgeLengths[tuple(sorted((u_name,v_name)))] = w
    
    def AddQForEdge(self,u_name,v_name,Q):
        self.QForEdge[(u_name,v_name)] = Q 
    
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
        
    def AddEdgeLogLikelihoods(self,edgeLogLikelihoodsDic):
        self.edgeLogLikelihoods.update(edgeLogLikelihoodsDic)
        
    def GetEdgeLogLikelihoods(self,p_name,c_name):
        return self.edgeLogLikelihoods[(p_name,c_name)]
    
    def ComputeEdgeBasedSitePatterns(self):        
        def ComputeSitePatternForEdge(u_name,v_name):
            v_smaller_name, v_larger_name = sorted(tuple((u_name,v_name)))  
            seq_s = self.sequences[v_smaller_name]
            seq_l = self.sequences[v_larger_name]
            sitePatternCount = array([[0]*4]*4)
            for site in range(self.sequenceLength):
                sitePatternCount[DNA.index(seq_s[site])][DNA.index(seq_l[site])] += self.siteWeight[site]
                    
            self.edgeBasedSitePatterns[(v_smaller_name, v_larger_name)] = sitePatternCount.copy() 
        
        map(ComputeSitePatternForEdge,self.edgeLengths.keys())
        
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
    
    def GetRootedTreeUsingSiteWeights(self,u_name,v_name,t_u,t_v):
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
            RT.siteWeights = self.siteWeights[:]
        return RT
    
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
    
    def GetTreeRootedAtVertex(self, vertex_name):
        RT = RootedTree()        
        verticesToVisit = [self.GetVertex(vertex_name)]
        while len(verticesToVisit) > 0:
            p = verticesToVisit.pop()
            neighborsOfP = p.neighbors
            for c in neighborsOfP:
                if not RT.ContainsVertex(c.name):
                    t = self.GetEdgeLength(c.name, p.name)
                    RT.AddDirectedEdge(p.name, c.name, t)
                    verticesToVisit.append(c)
        RT.SetRootByVertexName(vertex_name)
        return RT
    
    def ComputeLogLikelihoodForRootingAtVertex(self, vertex_name):
        logLik = 0
        verticesVisited = []        
        verticesToVisit = [self.GetVertex(vertex_name)]
        while len(verticesToVisit) > 0:
            p = verticesToVisit.pop()
            neighborsOfP = p.neighbors
            for c in neighborsOfP:
                if not c.name in verticesVisited:
                    logLik += self.GetEdgeLogLikelihoods(p.name, c.name)
                    verticesVisited.append(c.name)
                    verticesToVisit.append(c)
        return logLik
    
    def ComputeDyadLogLikelihoodForRootingAlongEdgeUsingSiteWeights(self, edge, maxIterations = 20):        
        # E-M for root state, root probability, and edge transition matrices
        u_name, v_name = edge
        seq_u = self.sequences[u_name]
        seq_v = self.sequences[v_name]
        # Initialize root state
        seq_r = ''
        for site in range(len(seq_u)):
            if seq_u[site] == seq_v[site]:
                seq_r += seq_u[site]
            else:
                if uniform() > 0.5:
                    seq_r += seq_v[site]
                else:
                    seq_r += seq_u[site]
        logLikForRooting = -100000000
        logLikForRooting_prev = logLikForRooting - 10
        iteration = 0
        while logLikForRooting - logLikForRooting_prev > 1.0 and iteration < maxIterations:
            iteration += 1
            logLikForRooting_prev = logLikForRooting
            logLikForRooting = 0
            # M-step root probability
            p_r = self.GetBaseFreqUsingSiteWeights(seq_r)
            # M-step transition matrices
            P_u = self.GetMaxLikEstimateOfTransitionMatrixUsingSiteWeights(seq_r, seq_u)
            P_v = self.GetMaxLikEstimateOfTransitionMatrixUsingSiteWeights(seq_r, seq_v)                 
            # E-step root state
            rootCharList = list(seq_r)
            for site in range(len(seq_u)):
                logLik_state_tup = []
                for nuc_code in range(4):
                    if p_r[nuc_code] != 0:
                        try:
                            stateLogLik = log(p_r[nuc_code]) + log(P_u[nuc_code,DNA.index(seq_u[site])]) + log(P_v[nuc_code,DNA.index(seq_v[site])])
                        except:
                            print ('root prob', p_r[nuc_code])
                            print ('transition matrix 1')
                            print (P_u)
                            print ('transition matrix 2')
                            print (P_v )                               
                    else:
                        stateLogLik = 0.0
                    logLik_state_tup.append((stateLogLik,nuc_code))
                logLik_state_tup.sort(key=itemgetter(0), reverse=True)
                rootCharList[site] = DNA[logLik_state_tup[0][1]]
                logLikForRooting +=  logLik_state_tup[0][0]*self.siteWeights[site]
            seq_r = ''.join(rootCharList)
        self.dyadLogLikelihood[sorted(tuple((u_name,v_name)))]=logLikForRooting    
        return logLikForRooting
    
    def ComputeMaximumLikelihoodEdgeForRootingUsingSiteWeights(self):
        self.ComputeEdgeLogLikelihoodsUsingSiteWeights()
        edge = self.edgeLengths.keys()[0]
        rt = self.GetTreeRootedAtVertex(edge[0])
        self.logLikelihoodSum = {e:0 for e in self.edgeLogLikelihoods.keys()}
        for v in rt.GetPostOrderTraversalWithoutRoot():
            for c in v.children:
                self.logLikelihoodSum[(v.name,v.parent.name)] += self.logLikelihoodSum[(c.name,v.name)] + self.edgeLogLikelihoods[(c.name,v.name)]
        
        for v in rt.GetPreOrderTraversalWithoutLeaves():
            neighborsOfV = [u.name for u in self.GetVertex(v.name).neighbors]    
            for c in v.children:
                remainingNeighborsOfV = list(set(neighborsOfV)-set([c.name]))
                for n in remainingNeighborsOfV:
                    self.logLikelihoodSum[(v.name,c.name)] += self.logLikelihoodSum[(n,v.name)] + self.edgeLogLikelihoods[(n,v.name)]
                
        map(lambda edge: self.ComputeDyadLogLikelihoodForRootingAlongEdgeUsingSiteWeights(edge),self.edgeLengths.keys())
        logLikOfRootingAtEdge = [(edge,self.logLikelihoodSum[(edge[0],edge[1])]+self.logLikelihoodSum[(edge[1],edge[0])]+self.d(edge))+self.dyadLogLikelihood[sorted((tuple(edge[0],edge[1])))] for edge in self.edgeLengths.keys()]            
        logLikOfRootingAtEdge.sort(key=itemgetter(1), reverse=True)
        return logLikOfRootingAtEdge[0]
    
    def ComputeLogLikelihoodForRootingAlongEdge(self, u_name, v_name, maxIterations = 20):
        logLik = 0
        verticesVisited = [u_name,v_name]        
        u = self.GetVertex(u_name)
        v = self.GetVertex(v_name)
        verticesToVisit = [u,v]
        while len(verticesToVisit) > 0:
            p = verticesToVisit.pop()
            neighborsOfP = p.neighbors
            for c in neighborsOfP:
                if not c.name in verticesVisited:
                    logLik += self.GetEdgeLogLikelihoods(p.name, c.name)
                    verticesVisited.append(c.name)
                    verticesToVisit.append(c)
        # E-M for root state, root probability, and edge transition matrices
        seq_u = self.sequences[u_name]
        seq_v = self.sequences[v_name]
        # Initialize root state
        seq_r = ''
        for pos in range(len(seq_u)):
            if seq_u[pos] == seq_v[pos]:
                seq_r += seq_u[pos]
            else:
                if uniform() > 0.5:
                    seq_r += seq_v[pos]
                else:
                    seq_r += seq_u[pos]
        logLikForRooting = -100000000
        logLikForRooting_prev = logLikForRooting - 10
        iteration = 0
        while logLikForRooting - logLikForRooting_prev > 1.0 and iteration < maxIterations:
            iteration += 1
            logLikForRooting_prev = logLikForRooting
            logLikForRooting = 0
            # M-step root probability
            p_r = GetBaseFreq(seq_r)
            # M-step transition matrices
            P_u = GetMaxLikEstimateOfTransitionMatrix(seq_r, seq_u)
            P_v = GetMaxLikEstimateOfTransitionMatrix(seq_r, seq_v)     
            # E-step root state
            rootCharList = list(seq_r)
            for pos in range(len(seq_u)):
                if seq_u[pos] != seq_v[pos]:
                    logLik_state_tup = []
                    for nuc_code in range(4):
                        if p_r[nuc_code] != 0:
                            try:
                                stateLogLik = log(p_r[nuc_code]) + log(P_u[nuc_code,DNA.index(seq_u[pos])]) + log(P_v[nuc_code,DNA.index(seq_v[pos])])
                            except:
                                print ('root prob', p_r[nuc_code])
                                print ('transition matrix 1')
                                print (P_u)
                                print ('transition matrix 2')
                                print (P_v)                               
                        else:
                            stateLogLik = 0.0
                        logLik_state_tup.append((stateLogLik,nuc_code))
                    logLik_state_tup.sort(key=itemgetter(0), reverse=True)
                    rootCharList[pos] = DNA[logLik_state_tup[0][1]]
                    logLikForRooting +=  logLik_state_tup[0][0]
                else:
                    rootCharList[pos] = seq_u[pos]
            seq_r = ''.join(rootCharList)
        logLik += logLikForRooting
        return logLik
    
    def UpdateMSTWithEdgeWeights(self,edgeWeights,signifDigits=8,debug=False):
        edgeTuple = sorted(edgeWeights.items(), key=itemgetter(1))        
        for (u_name, v_name), w in edgeTuple:
            if not self.ContainsVertex(v_name): 
                self.AddVertex(v_name)
            if not self.ContainsVertex(u_name):
                self.AddVertex(u_name)                
            u = self.GetVertex(u_name)
            v = self.GetVertex(v_name)
            if debug:
                if u_name in ['l_29','l_45']:      
                    print ("rep of", u_name, "is", u.rep.name, "and has rank", u.rep.rank) 
                    print ("edge added is", (u_name, v_name))
                    print ("rep of", v_name, "is", v.rep.name, "and has rank", v.rep.rank)
                    print ("---")
                elif v_name in ['l_29','l_45']:
                    print ("rep of", v_name, "is", v.rep.name, "and has rank", v.rep.rank)
                    print ("edge added is", (u_name, v_name))
                    print ("rep of", u_name, "is", u.rep.name, "and has rank", u.rep.rank)
                    print ("---")
            if Find(u)!=Find(v):
                self.AddEdge(u_name, v_name, w)
                Union(u, v)               
    
    def GetMLTreeRootedAlongEdge(self, u_name, v_name, maxIterations = 20):
        RT = RootedTree()
        verticesVisited = [u_name,v_name]        
        u = self.GetVertex(u_name)
        v = self.GetVertex(v_name)
        verticesToVisit = [u,v]
        while len(verticesToVisit) > 0:
            p = verticesToVisit.pop()
            neighborsOfP = p.neighbors
            for c in neighborsOfP:
                if not c.name in verticesVisited:
                    length = self.GetEdgeLength(p.name, c.name)
                    RT.AddDirectedEdge(p.name, c.name, length)
                    verticesVisited.append(c.name)
                    verticesToVisit.append(c)
        
        # E-M for root state, root probability, end edge transition matrices
        seq_u = self.sequences[u_name]
        seq_v = self.sequences[v_name]        
        # Initialize root state
        seq_r = ''
        for pos in range(len(seq_u)):
            if seq_u[pos] == seq_v[pos]:
                seq_r += seq_u[pos]
            else:
                if uniform() > 0.5:
                    seq_r += seq_v[pos]
                else:
                    seq_r += seq_u[pos]
        logLikForRooting = -100000000
        logLikForRooting_prev = logLikForRooting - 10
        iteration = 0
        while logLikForRooting - logLikForRooting_prev > 1.0 and iteration < maxIterations:
            iteration += 1
            logLikForRooting_prev = logLikForRooting
            logLikForRooting = 0
            # M-step root probability
            p_r = GetBaseFreq(seq_r)
            # M-step transition matrices   
            P_u = GetMaxLikEstimateOfTransitionMatrix(seq_r, seq_u)
            P_v = GetMaxLikEstimateOfTransitionMatrix(seq_r, seq_v)     
            # E-step root state
            rootCharList = list(seq_r)
            for pos in range(len(seq_u)):
                logLik_state_tup = []
                for nuc_code in range(4):
                    if p_r[nuc_code] != 0:
                        try:
                            stateLogLik = log(p_r[nuc_code]) + log(P_u[nuc_code,DNA.index(seq_u[pos])]) + log(P_v[nuc_code,DNA.index(seq_v[pos])])
                        except:
                            print ('root prob', p_r[nuc_code])
                            print ('transition matrix 1')
                            print (P_u)
                            print ('transition matrix 2')
                            print (P_v)                                
                    else:
                        stateLogLik = 0.0                    
                    logLik_state_tup.append((stateLogLik,nuc_code))
                logLik_state_tup.sort(key=itemgetter(0), reverse=True)
                rootCharList[pos] = DNA[logLik_state_tup[0][1]]
                logLikForRooting +=  logLik_state_tup[0][0]
            seq_r = ''.join(rootCharList)
        RT.AddDirectedEdge('h_root', u_name, hamDist(seq_r, seq_u))
        RT.AddDirectedEdge('h_root', v_name, hamDist(seq_r, seq_v))
        RT.SetRootByVertexName('h_root')
        
        return RT
    
    def GetMaximumLikelihoodRootedTrees(self, pool=None):
#         pool = None
        edgeList = self.edgeLengths.keys()
        if pool != None:
            logLikList = pool.map(poolFuncForComputingLogLikOfRootingAlongEdge,zip(repeat(self),edgeList))
        else:
            logLikList = map(poolFuncForComputingLogLikOfRootingAlongEdge,zip(repeat(self),edgeList))
        print (max(logLikList), min(logLikList))
        tup_edge_logLik = zip(edgeList,logLikList)
        tup_edge_logLik.sort(key= itemgetter(1), reverse = True)
        maxLogLik = tup_edge_logLik[0][1]
        edgesToRootAlong = [tup_edge_logLik[0][0]]
        for edge, logLik in tup_edge_logLik[1:]:
            if logLik == maxLogLik:
                edgesToRootAlong.append(edge)
        RTList = [self.GetMLTreeRootedAlongEdge(edge[0],edge[1]) for edge in edgesToRootAlong]
        print ("Number of ML trees is", len(RTList))
        return RTList
                
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
    
    def GetGraphDistance(self, vertexA_name, vertexB_name):
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
    
    def GetForestInducedByLeaves(self, leafNames):
        rootProb = RootedForest()
        hiddenVerticesVisitedOnce=[]        
        map(rootProb.AddVertex,leafNames)
        verticesToVisit=leafNames
        verticesVisited=leafNames
        while len(verticesToVisit)>0:
            v_name = verticesToVisit[0]
            del verticesToVisit[0]
            v = self.GetVertex(v_name)
            h = v.neighbors[0]
            w = self.GetEdgeLength(h.name, v_name)
            rootProb.AddEdge(h.name, v_name, w)
            verticesVisited.append(h.name)
            if h.name in hiddenVerticesVisitedOnce:
                verticesToVisit.append(h.name)
                hiddenVerticesVisitedOnce.remove(h.name)
            else:
                hiddenVerticesVisitedOnce.append(h.name)
                         
        return rootProb
    
    def GrowSubtreesUptoSizeK(self, k, subtreesSmallerThanK={}):
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
    
    def CanSubtreesBeExtractedForGivenTreeDepth(self, treeDepth):
        MSTSubtreesJustLargerThanK, _ = self.GrowSubtreesUptoSizeK(treeDepth, {})
        if len(MSTSubtreesJustLargerThanK) == 0:
            return False
        elif self.vertexCount == len(MSTSubtreesJustLargerThanK.values()[0]):
            return False
        else:
            return True
    
    def GetListOfBigVertexSetAndNestedVertexSets(self, treeDepth):
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
                print (len(MSTSubtreesJustLargerThanK))
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
    
    def GetMLRootedTreeForGMMUsingMultiThreadingForSiteLikelihoods(self,pool):
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
    
    def GetMLRootedTreeWithFixedQ(self,pool,Q):    
        maxLogLik = -10000000        
        RTList = []
        for (u_name, v_name) in self.edgeLengths.keys():
            t_uv = self.edgeLengths[tuple(sorted([u_name,v_name]))]
            rt = self.GetRootedTree(u_name, v_name, t_uv/2, t_uv/2)
            rt.Q = Q
            RTList.append(rt)
        
        RTList = pool.map(EMForMapWithFixedQ,RTList)
            
        for rt in RTList:
            if maxLogLik < rt.logLik:
                bestRT = rt
                maxLogLik = rt.logLik
        return bestRT             
    
    def GetMLRootedTreeWithFixedStatesAndEdgeLengthsWithEstimatedQ(self,pool):
        maxLogLik = -1*pow(10,10)
        RTList = []
        
        for (u_name, v_name) in self.edgeLengths.keys():
            t_uv = self.edgeLengths[tuple(sorted([u_name,v_name]))]
            rt = self.GetRootedTree(u_name, v_name, t_uv/2, t_uv/2)
            RTList.append(rt)
        RTList = pool.map(EMForMapWithFixedStatesAndEdgeLengthsAndEstimatedQ,RTList)
        
        for rt in RTList:
            if maxLogLik < rt.logLik:
                bestRT = rt
                maxLogLik = rt.logLik
        return bestRT
    
    def GetMLRootedTreeWithEstimatedQ(self,pool):
        RTList = []
        for (u_name, v_name) in self.edgeLengths.keys():
            t_uv = self.edgeLengths[tuple(sorted([u_name,v_name]))]
            rt = self.GetRootedTree(u_name, v_name, t_uv/2, t_uv/2)
            RTList.append(rt)
        if pool == None:
            RTList = map(EMForMapWithEstimatedQ,RTList)
        else:
            RTList = pool.map(EMForMapWithEstimatedQ,RTList)
        
        maxLogLik = max([rt.logLik for rt in RTList])
        
        MLRTs = []    
        for rt in RTList:
            if rt.logLik == maxLogLik:
                MLRTs.append(rt)
        
        MLRT = choice(MLRTs,size=1)[0]
        return MLRT
    
    def GetBaseFreqUsingSiteWeights(self,seq):
        baseFreq=[0.0]*4
        for site in range(len(seq)):
            baseFreq[DNA.index(seq[site])]+=float(self.siteWeights[site])
        return map(lambda x: x/float(sum(baseFreq)), baseFreq)
    
    def GetMaxLikEstimateOfTransitionMatrixUsingSiteWeights(self,seq_parent,seq_child):
        P = array([[pow(10,-5)]*4]*4)
        for site in range(len(seq_parent)):
            P[DNA.index(seq_parent[site]),DNA.index(seq_child[site])] += float(self.siteWeights[site])
        for row in range(4):
            rowSum = sum(P[row,:]) 
            for col in range(4):
                P[row,col] /= rowSum
        return P
    
    def GetSubstitutionCountMatrixUsingSiteWeights(self,seq_parent,seq_child):
        C = array([[0]*4]*4)
        for site in range(len(seq_parent)):
            C[DNA.index(seq_parent[site]),DNA.index(seq_child[site])] += self.siteWeights[site]
        return C
    
    def ComputeEdgeLogLikelihoodsViaEM(self,u_name,v_name):
        self.sequenceLength = len(self.sequences.values()[0])
        if u_name == "" and v_name == "":
            u_name, v_name = choice(self.edgeLengths.keys())
        edge_curr = sorted(tuple((u_name,v_name)))
        rt_curr = self.GetRootedTreeUsingSiteWeights(edge_curr[0], edge_curr[1], 0.01, 0.01)                
        rt_curr.PerformEMForGMMUsingSiteWeights()
        logLik_curr = rt_curr.logLik
        self.sequences = {v_name:rt_curr.sequences[v_name] for v_name in self.vertices.keys()}        
        edge_reRooting, logLik_reRooting = self.ComputeMaximumLikelihoodEdgeForRootingUsingSiteWeights() 
        while logLik_reRooting > logLik_curr and edge_reRooting != edge_curr:
            logLik_curr = logLik_reRooting
            edge_curr = edge_reRooting[:]  
            rt_curr = self.GetRootedTreeUsingSiteWeights(edge_curr[0], edge_curr[1], 0.01, 0.01)
            rt_curr.PerformEMForGMMUsingSiteWeights()
            self.sequences = {v_name:rt_curr.sequences[v_name] for v_name in self.vertices.keys()}        
            edge_reRooting, logLik_reRooting = self.ComputeMaximumLikelihoodEdgeForRootingUsingSiteWeights() 
        return None
    
    def GetMLRootedTreeUsingEM_v2(self,u_name,v_name):
        self.sequenceLength = len(self.sequences.values()[0])
        edge_curr = sorted(tuple((u_name,v_name)))
        rt_curr = self.GetRootedTreeUsingSiteWeights(edge_curr[0], edge_curr[1], 0.01, 0.01)                
        rt_curr.PerformEMForGMMUsingSiteWeights()
        logLik_curr = rt_curr.logLik
        self.sequences = {v_name:rt_curr.sequences[v_name] for v_name in self.vertices.keys()}        
        edge_reRooting, logLik_reRooting = self.ComputeMaximumLikelihoodEdgeForRootingUsingSiteWeights() 
        while logLik_reRooting > logLik_curr and edge_reRooting != edge_curr:
            logLik_curr = logLik_reRooting
            edge_curr = edge_reRooting[:]  
            rt_curr = self.GetRootedTreeUsingSiteWeights(edge_curr[0], edge_curr[1], 0.01, 0.01)
            rt_curr.PerformEMForGMMUsingSiteWeights()
            self.sequences = {v_name:rt_curr.sequences[v_name] for v_name in self.vertices.keys()}        
            edge_reRooting, logLik_reRooting = self.ComputeMaximumLikelihoodEdgeForRootingUsingSiteWeights() 
        return rt_curr
    
    def ComputeEdgeLogLikelihoodsUsingSiteWeights(self):
        def ComputeLogLikelihoodsPerEdge(edge):
            p_name, c_name = edge
            seq_p = self.sequences[p_name]
            seq_c = self.sequences[c_name]
            P_pc = self.GetMaxLikEstimateOfTransitionMatrixUsingSiteWeights(seq_p, seq_c)
            P_cp = self.GetMaxLikEstimateOfTransitionMatrixUsingSiteWeights(seq_c, seq_p)
            C_pc = self.GetSubstitutionCountMatrixUsingSiteWeights(seq_p, seq_c)
            C_cp = self.GetSubstitutionCountMatrixUsingSiteWeights(seq_c, seq_p)
            logLik_pc = 0
            logLik_cp = 0
            for i in range(4):
                for j in range(4):
                    logLik_pc += C_pc[i,j]*log(P_pc[i,j])
                    logLik_cp += C_cp[i,j]*log(P_cp[i,j]) 
            self.edgeLogLikelihoods[(p_name,c_name)] = logLik_pc
            self.edgeLogLikelihoods[(c_name,p_name)] = logLik_cp
        map(ComputeLogLikelihoodsPerEdge,self.edgeLengths.keys())
    
    def ComputeEdgeLogLikelihoodsUsingTwoDirectionalTransitionMatrices(self):
        def ComputeLogLikelihoodsPerEdge(edge):
            p_name, c_name = edge
            seq_p = self.sequences[p_name]
            seq_c = self.sequences[c_name]
            P_pc = GetMaxLikEstimateOfTransitionMatrix(seq_p, seq_c)
            P_cp = GetMaxLikEstimateOfTransitionMatrix(seq_c, seq_p)
            C_pc = GetSubstitutionCountMatrix(seq_p, seq_c)
            C_cp = GetSubstitutionCountMatrix(seq_c, seq_p)
            logLik_pc = 0
            logLik_cp = 0  
            for i in range(4):
                for j in range(4):
                    logLik_pc += C_pc[i,j]*log(P_pc[i,j])
                    logLik_cp += C_cp[i,j]*log(P_cp[i,j]) 
            self.edgeLogLikelihoods[(p_name,c_name)] = logLik_pc
            self.edgeLogLikelihoods[(c_name,p_name)] = logLik_cp
        map(ComputeLogLikelihoodsPerEdge,self.edgeLengths.keys())                

    def GetMLRootedGLFWithEstimatedQ(self,pool):    
        maxLogLik = -1*pow(10,10)    
        RTList = []
        for (u_name, v_name) in self.edgeLengths.keys():
            t_uv = self.edgeLengths[tuple(sorted([u_name,v_name]))]
            rt = self.GetRootedTree(u_name, v_name, t_uv/2, t_uv/2)
            RTList.append(rt)
        RTList = pool.map(EMForMapWithEstimatedQ,RTList)
                    
        for rt in RTList:
            if maxLogLik < rt.logLik:
                bestRT = rt
                maxLogLik = rt.logLik
                
        # Compute GLF by pruning long edges and contracting short edges s.t. BIC improves.
        
        return bestRT
            
    def WriteToFile(self,fileName):
        edgeListFile = open(fileName,'w')
        for u_name, v_name in self.edgeLengths.keys():
            edgeListFile.write(u_name+'\t'+v_name+'\t'+str(self.edgeLengths[tuple(sorted((u_name,v_name)))])+'\n')
        edgeListFile.close()


def RemoveHiddenVerticesOfDeg1And2(t,namesOfVerticesToInclude):
    hiddenVerticesOfDeg1=[]
    hiddenVerticesOfDeg2=[]
    for v in t.vertices.values():
        if v.name.startswith('hiddenVertex') and v.name not in namesOfVerticesToInclude:
            if v.degree==1:
                hiddenVerticesOfDeg1.append(v)
            elif v.degree==2:
                hiddenVerticesOfDeg2.append(v) 
    while len(hiddenVerticesOfDeg1)>0:
        h = hiddenVerticesOfDeg1[0]
        del hiddenVerticesOfDeg1[0]
        u = h.neighbors[0]
        t.RemoveVertex(h.name)
        if u.name.startswith('hiddenVertex') and u.name not in namesOfVerticesToInclude:
            if u.degree==1:
                hiddenVerticesOfDeg1.append(u)
                if u in hiddenVerticesOfDeg2:
                    hiddenVerticesOfDeg2.remove(u)
            elif u.degree==2 and u not in hiddenVerticesOfDeg2 and u.name not in namesOfVerticesToInclude:
                hiddenVerticesOfDeg2.append(u)
    # contract path containing degree 2 hidden vertices
    for h in hiddenVerticesOfDeg2:
        u, v = h.neighbors
        w1 = t.GetEdgeLength(u.name,h.name)
        w2 = t.GetEdgeLength(v.name,h.name)
        t.RemoveEdge(u.name,h.name)
        t.RemoveEdge(v.name,h.name)
        t.AddEdge(u.name,v.name,w1+w2)
    return t

        
class VertexInRootedGraph:
    
    def __init__(self,name):
        self.name = name
        self.inDegree = 0
        self.outDegree = 0
        self.numberOfDescendants = 0
        self.parent = None
        self.children = []
        self.numerOfVerticesInPathToRootIncludingRoot = 1
        self.distanceToClosestLeaf = 0
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
        self.indelible_replicates = 1
        self.sequences = {}
        self.Q = array([[0.0]*4]*4)
        self.QForEdge = {}
        self.rootProb = array([0.0]*4)
        self.logLik = 0
        self.edgeLogLikelihoods = {}
        self.transitionMatrices = {}        
        self.rateVectorForEdge = {}
        self.rateCategoryForEdge = {}
        self.rateCategoryForVertex = {}
        self.totalNumberOfRateCategories = 0
        self.rateVectorForCat = {}         
        self.p_change = 0
        
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
    
    def AddQForEdge(self,u_name,v_name,Q):
        self.QForEdge[(u_name,v_name)] = Q
        
    def AddTransitionMatrices(self,transitionMatricesToAdd):
        self.transitionMatrices.update(transitionMatricesToAdd)     
        
    def AddEdgeLogLikelihoods(self,edgeLogLikelihoodsDic):
        self.edgeLogLikelihoods.update(edgeLogLikelihoodsDic)
    
    def GetEdgeLogLikelihoods(self,p_name,c_name):
        return self.edgeLogLikelihoods[(p_name,c_name)]             
    
    def RenameNonLeafVertices(self, h_ind):
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
    
    def RenameRoot(self, h_ind):
        root = self.root        
        oldRootName = root.name
        newRootName = 'h_'+str(h_ind+1) 
        root.name = newRootName
        self.vertices[root.name] = root
        rootChildrenCopy = root.children[:]       
        for c in rootChildrenCopy:        
            t = self.GetEdgeLength(oldRootName, c.name)
            self.transitionMatrices[(newRootName,c.name)] = self.transitionMatrices[(oldRootName,c.name)].copy()
            self.RemoveDirectedEdge(oldRootName, c.name)
            self.AddDirectedEdge(newRootName, c.name, t)
                        
        del self.vertices[oldRootName]        
        # rename sequences dictionary       
        if oldRootName in self.sequences.keys():
            self.sequences[newRootName] = self.sequences[oldRootName] 
            del self.sequences[oldRootName]
         
    def GetDescendants(self,v_name):
        v = self.GetVertex(v_name)
        return v.GetDescendantNames()

    def GetAllDescendants(self,v_name):
        v = self.GetVertex(v_name)
        return v.GetAllDescendantNames()
        
    def SetNumberOfDescendants(self):
        for l in self.GetLeaves():
            l.numberOfDescendants = 1            
        for v in self.GetPostOrderTraversalWithoutRoot()+[self.root]:
            for d in v.children:
                v.numberOfDescendants += d.numberOfDescendants                

    def GetSackinIndex(self, normalize = "max"):
        if self.root.numberOfDescendants == 0:
            self.SetNumberOfDescendants()
        n = float(self.root.numberOfDescendants)
        sackinIndex = float(sum([v.numberOfDescendants for v in self.GetPreOrderTraversalWithoutLeaves()])) 
        if normalize == "numberOfLeaves":
            sackinIndex /= n
        elif normalize == "max":
            sackinIndex /= ((n*(n+1)/2)-1)
        return sackinIndex

    def SelectChangeNodePairs(self):
        self.change_node_pairs = []
        for v in self.vertices:
            if v.outDegree == 2:
                c_l, c_r = v.children
                if c_l.outDegree == 2 and c_r.outDegree == 2:
                    self.change_node_pairs.append([c_l,c_r])

    def GenerateConvergentlyEvolvingMarkovModel(self,change_node_left, change_node_right, GC_diff):
        # assert (child_left.parent == child_right.parent)
        # parent_of_change_nodes = child_left.parent
        # # child_left, child_right = parent_of_change_nodes.children
        # change_node_left =  sample(child_left.children,1)
        # change_node_right =  sample(child_right.children,1)
        condition_not_satisfied = True

        while condition_not_satisfied:
            # pi_at = uniform(0.2,0.5 + GC_diff)
            pi_a = uniform(0.1,0.25 + GC_diff/2)
            pi_t = 0.5 + GC_diff - pi_a
            pi_g = uniform(0.1,0.25 - GC_diff/2)
            pi_c = 0.5 - GC_diff - pi_g
            if pi_t > 0.1 and pi_c > 0.1:
                condition_not_satisfied = False
            # pi_r = GenerateProbabilityDistributionWithMinVal(0.1)
            pi_r = [0]*4
            pi_r[DNA.index("A")] = pi_a
            pi_r[DNA.index("T")] = pi_t
            pi_r[DNA.index("G")] = pi_g
            pi_r[DNA.index("C")] = pi_c
            pi_conv = pi_r[:]
            change_a = uniform(0, GC_diff/2)
            change_t = GC_diff - change_a
            change_g = uniform(0, GC_diff/2)
            change_c = GC_diff - change_g
            pi_conv[DNA.index("A")] -= change_a
            pi_conv[DNA.index("T")] -= change_t
            pi_conv[DNA.index("G")] += change_g
            pi_conv[DNA.index("C")] += change_c     
            for pi_i in pi_conv:
                if pi_i > 0.1 and pi_i < 1:
                    condition_not_satisfied = False
                else:
                    condition_not_satisfied = True

        
        Q_root = GenerateQForStationaryDistribution(pi_r)
        # GC increases via parallel evolution
        Q_conv = GenerateQForStationaryDistribution(pi_conv)
        pi_1 = GetStationaryDistribution(Q_root)
        pi_2 = GetStationaryDistribution(Q_conv)
        gc_diff = pi_2[DNA.index("G")] + pi_2[DNA.index("C")] - pi_1[DNA.index("G")] - pi_1[DNA.index("C")]
        at_diff = pi_2[DNA.index("A")] + pi_2[DNA.index("T")] - pi_1[DNA.index("A")] - pi_1[DNA.index("T")]
        print("gc_diff is ", gc_diff)
        print("at_diff is ", at_diff)

        # Q_root, Q_conv        
        self.rateCategoryForVertex = {}        
        self.rateCategoryForEdge = {}
        self.rateVectorForEdge = {}
        self.totalNumberOfRateCategories = 0        
        self.rateCategoryForVertex = {self.root.name: 1, change_node_left.name: 2, change_node_right.name: 2}                
        self.rateVectorForCat[1] = Get11FreeRates(Q_root)
        self.rateVectorForCat[2] = Get11FreeRates(Q_conv)      
        parents = [self.root]
        while len(parents) > 0:
            p = parents.pop()
            if len(p.children) > 0:
                for c in p.children:
                    if c.name != change_node_left.name and c.name != change_node_right.name:
                        self.rateCategoryForVertex[c.name] = self.rateCategoryForVertex[p.name]                    
                    parents.append(c)
                    
    def GenerateConvergentlyEvolvingMarkovModel(self,change_node_left, change_node_right, GC_diff):
        # assert (child_left.parent == child_right.parent)
        # parent_of_change_nodes = child_left.parent
        # # child_left, child_right = parent_of_change_nodes.children
        # change_node_left =  sample(child_left.children,1)
        # change_node_right =  sample(child_right.children,1)
        condition_not_satisfied = True

        while condition_not_satisfied:
            # pi_at = uniform(0.2,0.5 + GC_diff)
            pi_a = uniform(0.1,0.25 + GC_diff/2)
            pi_t = 0.5 + GC_diff - pi_a
            pi_g = uniform(0.1,0.25 - GC_diff/2)
            pi_c = 0.5 - GC_diff - pi_g
            if pi_t > 0.1 and pi_c > 0.1:
                condition_not_satisfied = False
            # pi_r = GenerateProbabilityDistributionWithMinVal(0.1)
            pi_r = [0]*4
            pi_r[DNA.index("A")] = pi_a
            pi_r[DNA.index("T")] = pi_t
            pi_r[DNA.index("G")] = pi_g
            pi_r[DNA.index("C")] = pi_c
            pi_conv = pi_r[:]
            change_a = uniform(0, GC_diff/2)
            change_t = GC_diff - change_a
            change_g = uniform(0, GC_diff/2)
            change_c = GC_diff - change_g
            pi_conv[DNA.index("A")] -= change_a
            pi_conv[DNA.index("T")] -= change_t
            pi_conv[DNA.index("G")] += change_g
            pi_conv[DNA.index("C")] += change_c     
            for pi_i in pi_conv:
                if pi_i > 0.1 and pi_i < 1:
                    condition_not_satisfied = False
                else:
                    condition_not_satisfied = True

        
        Q_root = GenerateQForStationaryDistribution(pi_r)
        # GC increases via parallel evolution
        Q_conv = GenerateQForStationaryDistribution(pi_conv)
        pi_1 = GetStationaryDistribution(Q_root)
        pi_2 = GetStationaryDistribution(Q_conv)
        gc_diff = pi_2[DNA.index("G")] + pi_2[DNA.index("C")] - pi_1[DNA.index("G")] - pi_1[DNA.index("C")]
        at_diff = pi_2[DNA.index("A")] + pi_2[DNA.index("T")] - pi_1[DNA.index("A")] - pi_1[DNA.index("T")]
        print("gc_diff is ", gc_diff)
        print("at_diff is ", at_diff)

        # Q_root, Q_conv        
        self.rateCategoryForVertex = {}        
        self.rateCategoryForEdge = {}
        self.rateVectorForEdge = {}
        self.totalNumberOfRateCategories = 0        
        self.rateCategoryForVertex = {self.root.name: 1, change_node_left.name: 2, change_node_right.name: 2}                
        self.rateVectorForCat[1] = Get11FreeRates(Q_root)
        self.rateVectorForCat[2] = Get11FreeRates(Q_conv)      
        parents = [self.root]
        while len(parents) > 0:
            p = parents.pop()
            if len(p.children) > 0:
                for c in p.children:
                    if c.name != change_node_left.name and c.name != change_node_right.name:
                        self.rateCategoryForVertex[c.name] = self.rateCategoryForVertex[p.name]                    
                    parents.append(c)

    def GenerateMarkovModulatedMarkovModel(self,p_change):
        rateMatricesForVertices = {}
        rateMatricesForEdges = {}
        self.rateCategoryForVertex = {}        
        self.rateCategoryForEdge = {}
        self.rateVectorForEdge = {}
        self.totalNumberOfRateCategories = 0
        Q = GenerateRandomQ()
        cat = 1
        self.rateCategoryForVertex = {self.root.name: cat}
        rateMatricesForVertices[self.root.name] = Q
        self.rateVectorForCat[cat] = Get11FreeRates(Q)
        parents = [self.root]        
        while len(parents) > 0:
            p = parents.pop()
            if len(p.children) > 0:
                for c in p.children:
                    p_sampled = uniform(0,1)
                    if p_change > p_sampled:                
                        Q = GenerateRandomQ()
                        cat += 1
                        self.totalNumberOfRateCategories += 1
                        self.rateVectorForCat[cat] = Get11FreeRates(Q)
                    else:
                        Q = rateMatricesForVertices[p.name][:]
                    rateMatricesForVertices[c.name] = Q
                    rateMatricesForEdges[(p.name,c.name)] = Q
                    self.rateCategoryForVertex[c.name] = cat                    
                    self.rateCategoryForEdge[(p.name,c.name)] = cat                    
                    self.rateVectorForEdge[(p.name,c.name)] = Get11FreeRates(Q)
                    parents.append(c)
    
    def GetBranchModelForINDELIBLE(self):
        labelAtVertex = {v : "" for v in self.vertices.values() if v.outDegree > 0}
        timesVisited = {v : 0 for v in labelAtVertex.keys() if v.outDegree > 0}
        labelAtVertex.update({v : v.name for v in self.GetLeaves()})    
        postOrderTraversalWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        for c in postOrderTraversalWithoutRoot:            
            p = c.parent
            timesVisited[p] += 1
            if timesVisited[p] == 1:
                labelAtVertex[p] += "("+labelAtVertex[c]+"#M"+str(self.rateCategoryForVertex[c.name])
            else :
                labelAtVertex[p] += ","+labelAtVertex[c]+"#M"+str(self.rateCategoryForVertex[c.name])
            if p.outDegree - timesVisited[p] == 0:
                labelAtVertex[p] += ")"
        labelAtRoot = labelAtVertex[self.GetRoot()]+"#M" + str(self.rateCategoryForVertex[self.root.name]) + ";"
        return(labelAtRoot)        
    
    def WriteControlFile(self):
        controlFile = open(self.control_file_name,"w")
        sequence_file_suffix = self.sequence_file_suffix
        # print(sequence_file_name)
        # TYPE block
        controlFile.write('[TYPE] NUCLEOTIDE 1\n')
        
        # Model block        
        for rateCat in self.rateVectorForCat.keys():
            rateVector = self.rateVectorForCat[rateCat]
            model_name = "M"+str(rateCat)
            controlFile.write('\n[MODEL] ' + model_name + '\n')
            controlFile.write(' [submodel] UNREST ')
            for rate in rateVector:
                controlFile.write(str(rate) + " ")
            controlFile.write("\n")
        
        # Tree block
        controlFile.write('\n[TREE] T ')
        newick_label = self.GetNewickLabel()                
        controlFile.write(newick_label + '\n')
        
        # Branch model block
        controlFile.write('\n[BRANCHES] B ')
        branch_model = self.GetBranchModelForINDELIBLE()
        controlFile.write(branch_model + '\n')

        # Partitions block
        controlFile.write('\n[PARTITIONS] P [T  B ' + str(self.sequenceLength) + ']\n')

        # Evolve block
        # controlFile.write('[EVOLVE] P  1 ' + sequence_file_name + '\n')
        controlFile.write('\n[EVOLVE] P '+ str(self.indelible_replicates) + '\t' + sequence_file_suffix + '\n')

        controlFile.close()
    def SimulateEvolutionUsingIndelible(self):
        devnull=open(os.devnull,'w')
        self.WriteControlFile()
        controlDirectory = '/'.join(self.control_file_name.split('/')[:-1])
        # sub.call("cd " + controlDirectory,shell=True)
        sub.call("indelible",cwd = controlDirectory, shell=True)
        # sub.call("indelible",cwd = controlDirectory, shell=True,stdout=devnull)
        devnull.close()
    
    # For multifurcating trees, count no. of steps to vertex with outdegree k as log(k,base=2)
    # For a bifurcating tree, each step will have size log(2,2) which equals to one.
    def GetSackinIndexUsingInformationTheoreticMaxNLevels(self, normalize="max"):
        leaves = self.GetLeaves()
        n = float(len(leaves))
        rootName = self.root.name
        def ComputeDistanceToRoot(v):
            distanceToRoot = 0
            while v.name != rootName:
                v = v.parent
                distanceToRoot += log(v.outDegree,2) 
            return float(distanceToRoot)
        sackinIndex = sum(map(ComputeDistanceToRoot,leaves))  
        if normalize == "numberOfLeaves":
            sackinIndex /= n
        elif normalize == "max":
            sackinIndex /= ((n*(n+1)/2)-1)
        return sackinIndex
        
#     def GetCollessIndex(self):
#         self.SetNumberOfDescendants()
#         for v in self.GetPreOrderTraversalWithoutLeaves():       
            
    def GetDistancesFromRootToLeaves(self):
        distancesFromRoot = {l:0 for l in self.GetLeaves()}
        for leaf in distancesFromRoot.keys():
            c = leaf
            while c.inDegree != 0:                          
                distancesFromRoot[leaf] += self.GetEdgeLength(c.parent.name, c.name)
                c = c.parent
        return distancesFromRoot
            
        
    def GetAllNonTrivialClusters(self):
        clusters = set([])
        orderedVertices = self.GetPreOrderTraversalWithoutLeaves()[1:]
        for v in orderedVertices:
            clusters.update([','.join(sorted(self.GetDescendants(v.name)))])
        return clusters
    
    def GetAllNonTrivialClustersForFullyLabeledTree(self):
        clusters = set([])
        orderedVertices = self.GetPreOrderTraversalWithoutLeaves()[1:]
        for v in orderedVertices:
            clusters.update([','.join(sorted(self.GetAllDescendants(v.name)))])
        return clusters
    
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
    
    def GetUnrootedTree_v2(self):
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
#             print rootToRemove.degree
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
            T.siteWeights = self.siteWeights[:]                        
        
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
    
    def ReRoot(self,u_name,v_name): # don't use this
        root = self.root
        u = self.vertices[u_name]
        v = self.vertices[v_name]
        ancestor = u.parent
        pathFromNewRootToOldRoot = [ancestor]
        while ancestor != root:
            ancestor = ancestor.parent
            pathFromNewRootToOldRoot.append(ancestor)
        for ind in range(len(pathFromNewRootToOldRoot)-2):
            b = pathFromNewRootToOldRoot[ind]
            a = pathFromNewRootToOldRoot[ind+1]
            t_ab = self.GetEdgeLength(a.name, b.name)
            self.RemoveDirectedEdge(a.name, b.name)
            self.AddDirectedEdge(b.name, a.name, t_ab)
        print (map (lambda x: x.name, root.children))
        i, j = root.children
        t_ij = self.GetEdgeLength(root.name, i.name) + self.GetEdgeLength(root.name, j.name)
        self.RemoveDirectedEdge(root.name, i.name)
        self.RemoveDirectedEdge(root.name, j.name)
        t_uv = self.GetEdgeLength(u.name, v.name)
        self.RemoveDirectedEdge(u.name, v.name)
        self.AddDirectedEdge(root.name, u.name,t_uv/float(2))
        self.AddDirectedEdge(root.name, v.name,t_uv/float(2))
        if i.parent == root:
            self.AddDirectedEdge(i.name, j.name, t_ij)
        elif j.parent == root:
            self.AddDirectedEdge(j.name, i.name, t_ij)
        elif i in pathFromNewRootToOldRoot:
            self.AddDirectedEdge(i.name, j.name, t_ij)
        else:
            self.AddDirectedEdge(j.name, i.name, t_ij) 
        
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
        parent = self.GetVertex(parent_name)
        child = self.GetVertex(child_name)
        parent.RemoveChild(child)
        child.RemoveParent()
        del self.edgeLengths[parent_name,child_name]
        if (parent_name,child_name) in self.transitionMatrices:
            del self.transitionMatrices[(parent_name,child_name)]
        
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
    
    def SuppressAllRootLocations(self):
        pass
        
    def ContractAShortestEdgeIncidentToTheRootForGMM(self):
        root = self.root        
        c_l, c_r = root.children        
        t_l = self.GetEdgeLength(root.name, c_l.name)
        t_r = self.GetEdgeLength(root.name, c_r.name)
        self.RemoveDirectedEdge(root.name, c_l.name)
        self.RemoveDirectedEdge(root.name, c_r.name)        
        
        if t_l < t_r:
            parent_name = c_l.name
            child_name = c_r.name            
        else:
            parent_name = c_r.name
            child_name = c_l.name
            
        seq_parent = self.sequences[parent_name]
        seq_child = self.sequences[child_name]
        P = self.GetMaxLikEstimateOfTransitionMatrixUsingSiteWeights(seq_parent, seq_child)
        t = self.hamDistUsingSiteWeights(seq_parent,seq_child)
        self.AddDirectedEdge(parent_name, child_name, t)
        self.transitionMatrices[(parent_name, child_name)] = P
        self.RemoveVertex(root.name)        
        return None

    def NumberOfConflictingEdges(self, edge):
        edgesInOriginalTree = set(self.edgeLengths.keys())-set([edge])
        T = Tree()
        T.AddEdges(self.edgeLengths.copy())
        u_name, v_name = edge
        t = T.GetEdgeLength(u_name, v_name)
        RT_edge = T.GetRootedTree(u_name, v_name, t/2.0, t/2.0)
        RT_edge.ContractAShortestEdgeIncidentToTheRoot()
        edgesInReorientedTree = set(RT_edge.edgeLengths.keys())-set([edge])
        return len(edgesInOriginalTree-edgesInReorientedTree)

    def GetTopNLeastConflictingEdges(self,N):
        tup_edge_noOfConflictingEdges = [(edge,self.NumberOfConflictingEdges(edge)) for edge in self.edgeLengths.keys()]
        tup_edge_noOfConflictingEdges.sort(key=itemgetter(1))
        edgesToReturn = [tup_edge_noOfConflictingEdges[i][0] for i in range(N)]
        return edgesToReturn
        

    def GetDistance(self,u_name,v_name):
        label_u = self.vertexLabels[u_name]
        label_v = self.vertexLabels[v_name]
        d = 0
        for i in range(len(label_u)):
            if label_u[i] != label_v[i]:
                d += 1
        return d
    
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

    def GetMPStateChanges(self,leafSequences):
        sequenceLength = len(leafSequences.values()[0])
        leaves = self.GetLeaves()        
        postOrderVerticesWithoutRoot = self.GetPostOrderTraversalWithoutRoot()                
        def CountMinStateChangesForSite(site):
            minStateChanges = 0  
            statesForSite = {v : set([leafSequences[v.name][site]]) for v in leaves}
            for v in postOrderVerticesWithoutRoot:
                if v.parent in statesForSite:
                    if len(statesForSite[v.parent] & statesForSite[v]) == 0:
                        statesForSite[v.parent].update(statesForSite[v])
                        minStateChanges += 1
                    else:
                        statesForSite[v.parent] = statesForSite[v.parent] & statesForSite[v]
                else:
                    statesForSite[v.parent] = set(list(statesForSite[v]))
            return minStateChanges
        numberOfStateChanges = sum(map(CountMinStateChangesForSite,range(sequenceLength)))
        return numberOfStateChanges

# Initialization   
    def EstimateNonLeafSequencesUsingMPFitch(self):
        allSequences = self.sequences            
        leaves = self.GetLeaves()
        root = self.GetRoot()
        postOrderVerticesWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        preOrderNonLeafVertices = self.GetPreOrderTraversalWithoutLeaves()
        nonLeafSequences = {v.name : "" for v in preOrderNonLeafVertices}
        
        def GetStatesForSite(site):
            statesForSite = {v : set([allSequences[v.name][site]]) for v in leaves}
            for v in postOrderVerticesWithoutRoot:            
                if v.parent in statesForSite:
                    if len(statesForSite[v.parent] & statesForSite[v]) == 0:
                        statesForSite[v.parent].update(statesForSite[v])
                    else:
                        statesForSite[v.parent] = statesForSite[v.parent] & statesForSite[v]
                else:
                    statesForSite[v.parent] = set(list(statesForSite[v]))
            statesForSite[root] = choice(list(statesForSite[root]),1)[0]
            nonLeafSequences[root.name]+=statesForSite[root]
            for v in preOrderNonLeafVertices[1:]:            
                if statesForSite[v.parent] in statesForSite[v]:
                    statesForSite[v] = statesForSite[v.parent]
                else:
                    statesForSite[v] = choice(list(statesForSite[v]),1)[0]
                nonLeafSequences[v.name]+=statesForSite[v]
            return 0
        map(GetStatesForSite,range(self.sequenceLength))
        self.sequences.update(nonLeafSequences)
        return None
    
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
        
        map(GetStatesForSite,range(len(self.sequences.values()[0])))
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
    
    def OptimizeRootProbUsingSiteWeights(self):
        root = self.GetRoot()
        rootSeq = self.sequences[root.name]        
        rootProb = [0.0] * 4
        for siteInd in range(len(rootSeq)):
            rootProb[DNA.index(rootSeq[siteInd])] += self.siteWeights[siteInd]
            
        rootProb_sum = float(sum(rootProb))
        for i in range(4):
            rootProb[i] /= rootProb_sum
        self.rootProb = rootProb
        return None
            
    def InitializeQ(self,skipRoot=False):
        allSequences = self.sequences
        substitutionMatrix = array([[0.0]*4]*4)
        orderedVertices = self.GetPreOrderTraversalWithoutLeaves()
        if skipRoot:
            orderedVertices = orderedVertices[1:][:]
        for parent in orderedVertices:
            seq_p = allSequences[parent.name]
            for c in parent.children:                
                seq_c = allSequences[c.name]                 
                for siteInd in range(len(self.firstSiteOfEachDistinctPattern)):
                    char_p = seq_p[self.firstSiteOfEachDistinctPattern[siteInd]]
                    char_c = seq_c[self.firstSiteOfEachDistinctPattern[siteInd]]
                    substitutionMatrix[DNA.index(char_p),DNA.index(char_c)] += self.patternRepeats[siteInd]
        for nuc in range(4):
            substitutionMatrix[nuc,:] /= float(sum(substitutionMatrix[nuc,:]))
        Q = logm(substitutionMatrix).real
        for i in range(4):
            for j in range(4):
                if i!=j:
                    Q[i,j] = max(pow(10,-5),Q[i,j])
        self.Q = NormalizeQMatrix(Q)
        return None
    
    
    def InitializeTransitionMatricesAndEdgeLengths(self):
        for (u_name, v_name) in self.transitionMatrices.keys():            
            seq_u = self.sequences[u_name]
            seq_v = self.sequences[v_name]
            P = GetMaxLikEstimateOfTransitionMatrix(seq_u, seq_v)
            self.transitionMatrices[(u_name,v_name)] = P
            t = hamDist(seq_u,seq_v)
            self.edgeLengths[(u_name,v_name)] = t
    
    def EstimateTransitionMatrices(self):
        def EstimateTransitionMatrix(u_name,v_name):
            seq_u = self.sequences[u_name]
            seq_v = self.sequences[v_name]
            P = GetMaxLikEstimateOfTransitionMatrix(seq_u, seq_v)
            self.transitionMatrices[(u_name,v_name)] = P
        
        map(EstimateTransitionMatrix,self.transitionMatrices.keys())
    
    def EstimateEdgeLengths(self):        
        def ComputeAndStoreEdgeLengths(u_name,v_name):
            seq_u = self.sequences[u_name]
            seq_v = self.sequences[v_name]
            t = hamDist(seq_u,seq_v)
            self.edgeLengths[(u_name,v_name)] = t
        
        map(ComputeAndStoreEdgeLengths,self.transitionMatrices.keys())
        
# Expectation
    def ComputeExpectedAncestralStatesAndUpdateLogLik(self):
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
        Q = self.Q
        postOrderVerticesWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        preOrderTraversalWithoutLeaves = self.GetPreOrderTraversalWithoutLeaves()
        self.logLik = 0
        rootProb = self.rootProb
        
        def ComputeExpectedStatesAndUpdateLogLikPerSite(site):
            conditionalLikelihoods = {}        
            for v in leaves:
                conditionalLikelihoods[v] = [0.0]*4
                conditionalLikelihoods[v][DNA.index(allSequences[v.name][site])] = 1.0 
            for v in nonLeafVertices:
                conditionalLikelihoods[v] = [1.0]*4
            for c in postOrderVerticesWithoutRoot:                              
                p = c.parent
                t = self.GetEdgeLength(p.name, c.name)
                P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q,t)         
                for nuc_p in range(4):
                    partialLikelihood_c_nuc_p = 0.0
                    for nuc_c in range(4):
                        partialLikelihood_c_nuc_p += P[nuc_p,nuc_c]*conditionalLikelihoods[c][nuc_c]
                    conditionalLikelihoods[p][nuc_p] *= partialLikelihood_c_nuc_p                    
            lik = 0
            for nuc_root in range(4):
                lik += rootProb[nuc_root]*conditionalLikelihoods[root][nuc_root]            
            self.logLik += log(max(lik,pow(10,-10))) # multiply by site weight
                
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
                ancestralStates[p.name] += stateWithMaxProb
                for c in p.children:
                    if c.outDegree > 0:
                        t = self.GetEdgeLength(p.name, c.name)
                        P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q,t)
                        priorProbs[c] = P[DNA.index(stateWithMaxProb),:]                    
            return None
        
        map(ComputeExpectedStatesAndUpdateLogLikPerSite,range(self.sequenceLength))
        for v_name in ancestralStates.keys():
            self.sequences[v_name] = ancestralStates[v_name]
        return None

    def ComputeExpectedAncestralStatesAndUpdateLogLikUsingPatternSaving(self):
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
        Q = self.Q
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
                t = self.GetEdgeLength(p.name, c.name)
                P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q,t)         
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
                        t = self.GetEdgeLength(p.name, c.name)
                        P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q,t)
                        priorProbs[c] = P[DNA.index(stateWithMaxProb),:]                    
            return None
        
        map(ComputeExpectedStatesAndUpdateLogLikPerSite,range(len(self.firstSiteOfEachDistinctPattern)))
        for v_name in ancestralStates.keys():                        
            self.sequences[v_name] = self.OrderCharsWrtSitePatternRepeatBlocks(ancestralStates[v_name])
        return None
    
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
    
    def ComputeExpectedAncestralStatesAndUpdateLogLikUsingSiteWeights(self):        
        root = self.GetRoot()
        leaves = []
        nonLeafVertices = []
        for v in self.vertices.values():
            if v.outDegree == 0:
                leaves.append(v)
            else:
                nonLeafVertices.append(v)
                self.sequences[v.name]=""
        priorProbs = {}
        priorProbs[root] = self.rootProb        
        postOrderVerticesWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        preOrderTraversalWithoutLeaves = self.GetPreOrderTraversalWithoutLeaves()
        self.logLik = 0
        rootProb = self.rootProb
        
        def ComputeExpectedStatesAndUpdateLogLikPerSite(site):
            conditionalLikelihoods = {}        
            for v in leaves:
                conditionalLikelihoods[v] = [0.0]*4
                conditionalLikelihoods[v][DNA.index(self.sequences[v.name][site])] = 1.0 
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
            self.logLik += log(max(lik,pow(10,-10)))*self.siteWeights[site]
                
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
                self.sequences[p.name] += stateWithMaxProb
                for c in p.children:
                    if c.outDegree > 0:
                        P = self.transitionMatrices[(p.name,c.name)]
                        priorProbs[c] = P[DNA.index(stateWithMaxProb),:]                    
            return None
        
        map(ComputeExpectedStatesAndUpdateLogLikPerSite,range(len(self.siteWeights)))        
        return None
    
    
    def ComputeExpectedAncestralStatesAndUpdateLogLikUsingPatternSavingForTransitionMatricesUsingPool(self,pool):
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
        
        pool.map(ComputeExpectedStatesAndUpdateLogLikPerSite,range(len(self.firstSiteOfEachDistinctPattern)))
        for v_name in ancestralStates.keys():                        
            self.sequences[v_name] = self.OrderCharsWrtSitePatternRepeatBlocks(ancestralStates[v_name])
        return None

# Maximization
    def OptimizeEdgeLengths(self):
        Q = self.Q
        allSequences = self.sequences
        def OptimizeEdgeLength(edge):
            u_name, v_name = edge
            sequence_u = allSequences[u_name]
            sequence_v = allSequences[v_name]        
            def ComputeNegLogLik(t):
                P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q,t)            
                def GetLogLikPerSite(siteInd):
                    site = self.firstSiteOfEachDistinctPattern[siteInd]
                    numberOfRepeats = self.patternRepeats[siteInd]                    
                    ind_u = DNA.index(sequence_u[site])
                    ind_v = DNA.index(sequence_v[site])
                    return log(P[ind_u,ind_v])*numberOfRepeats
                return -1*sum(map(GetLogLikPerSite,range(len(self.firstSiteOfEachDistinctPattern))))
            t_opt = minimize_scalar(ComputeNegLogLik, bounds = (pow(10,-5), pow(10,1)), method="bounded").x
            self.UpdateEdgeLength(u_name, v_name, t_opt)
        map(OptimizeEdgeLength,self.GetEdgeLengths().keys())
        return None

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

    def GetEuclideanDistanceBetweenRateMatrices(self,Q_1,Q_2):
        d = 0
        for i in range(4):
            for j in range(4):
                d += pow(Q_1[i,j]-Q_2[i,j],2)
        d = cmath.sqrt(d)
        d /= 4.0
        return d
       
    def OptimizeQ(self):
        allSequences = self.sequences
        rootName = self.root.name
        Q_vector_0 = ComputeQVectorOfOffDiagonalElementsFromQMatrix(self.Q)
        edgesToUse= []
        for edge in self.edgeLengths.keys():
            if edge[0] != rootName:
                edgesToUse.append(edge)
        def ComputeNegLogLik(Q_vector_offDiag):
            Q = ComputeQMatrixFromVectorOfOffDiagonalElements(Q_vector_offDiag)
            def ComputeLogLikPerEdge(edge):
                u_name, v_name = edge            
                t = self.GetEdgeLength(u_name, v_name)
                P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q,t)
                sequence_u = allSequences[u_name]
                sequence_v = allSequences[v_name]
    
                def ComputeLogLikPerEdgePerSite(siteInd):
                    site = self.firstSiteOfEachDistinctPattern[siteInd]
                    numberOfRepeats = self.patternRepeats[siteInd]
                    ind_u = DNA.index(sequence_u[site])
                    ind_v = DNA.index(sequence_v[site])
                    return log(P[ind_u,ind_v])*numberOfRepeats          
    
                return sum(map(ComputeLogLikPerEdgePerSite,range(len(self.firstSiteOfEachDistinctPattern))))
            negLogLik = -1*sum(map(ComputeLogLikPerEdge,edgesToUse))
            # rootProb not included since likelihood contribution from rootProb does not change as Q is optimized.
            return negLogLik
        optimizationResult = minimize(ComputeNegLogLik, Q_vector_0, method="L-BFGS-B", bounds = [(pow(10,-3), None)]*len(Q_vector_0), options={'ftol': 1e-3, 'disp': False})                               
        Q_vector_opt = optimizationResult.x
        Q_opt = ComputeQMatrixFromVectorOfOffDiagonalElements(Q_vector_opt)        
        self.Q = Q_opt
        return None
    
    def OptimizeTransitionMatrices(self):
        self.InitializeTransitionMatricesAndEdgeLengths()    
    
    def GetMaxLikEstimateOfTransitionMatrixUsingSiteWeights(self,seq_parent,seq_child):
        P = array([[pow(10,-5)]*4]*4)
        for site in range(len(seq_parent)):
            P[DNA.index(seq_parent[site]),DNA.index(seq_child[site])] += float(self.siteWeights[site])
        for row in range(4):
            rowSum = sum(P[row,:]) 
            for col in range(4):
                P[row,col] /= rowSum
        return P
    
    def OptimizeTransitionMatricesUsingSiteWeights(self):
        for (u_name, v_name) in self.edgeLengths.keys():            
            seq_u = self.sequences[u_name]
            seq_v = self.sequences[v_name]
            P = self.GetMaxLikEstimateOfTransitionMatrixUsingSiteWeights(seq_u, seq_v)
            self.transitionMatrices[(u_name,v_name)] = P        
                    
    def PerformEMForGMMUsingPatternSavings(self,maxIterations = 20, logLikImprovementThreshold = 1.0, verbose = False):
        self.InitializeNonLeafSequencesUsingMPUsingPatternSavings()
        # Initialized ancestral states using MP
        self.InitializeRootProbUsingPatternSavings()
        # Initialized root probability
        self.InitializeTransitionMatricesAndEdgeLengths()
        # Initialized transition matrices
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
        # store edgeLogLikelihoods
    
    def PerformEMForGMMUsingSiteWeights(self,maxIterations = 20, logLikImprovementThreshold = 1.0, verbose = False):
        # Initialize ancestral states using MP
        self.InitializeNonLeafSequencesUsingMPHartigan()
        # Initialize root probability
        self.OptimizeRootProbUsingSiteWeights()
        # Initialize transition matrices
        self.OptimizeTransitionMatricesUsingSiteWeights()        
        logLikConvergenceNotReached = True
        iteration = 1
        currentLogLik = -1*pow(10,10)   
        current_rootSeq = self.sequences[self.root.name]
        while logLikConvergenceNotReached and iteration <= maxIterations:
            # Expectation step
            self.ComputeExpectedAncestralStatesAndUpdateLogLikUsingSiteWeights()         
            if verbose:
                print ("Expected states computed")
                changeInRootSeq = hamDist(current_rootSeq, self.sequences[self.root.name])
                current_rootSeq = self.sequences[self.root.name]
            # Maximization steps
            # Given root sequence, optimize rootProb            
            self.OptimizeRootProbUsingSiteWeights()
            if verbose:
                print ("Optimized root probability")                                   
            # Given t and rootProb, optimize Q
            self.OptimizeTransitionMatricesUsingSiteWeights()
            if verbose:
                print ("Optimized transition matrices ")
                print ("iteration: "+str(iteration)+". Change in root seq is " + str(changeInRootSeq))                
            iteration += 1
            logLikConvergenceNotReached = self.logLik - currentLogLik > logLikImprovementThreshold            
            currentLogLik = self.logLik
            
    
    def PerformEMForGMMUsingMultiThreadingForSiteLikelihoods(self,pool,maxIterations = 20, logLikImprovementThreshold = 1.0, verbose = False):
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
            self.ComputeExpectedAncestralStatesAndUpdateLogLikUsingPatternSavingForTransitionMatricesUsingPool(pool)
            if verbose:
                changeInRootSeq = hamDist(current_rootSeq, self.sequences[self.root.name])
                current_rootSeq = self.sequences[self.root.name]
                print ("Expected states computed")
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
        # store edgeLogLikelihoods
    
    def ComputeEdgeLogLikelihoods(self):
        def ComputeLogLikelihoodsPerEdge(edge):
            p_name, c_name = edge
            seq_p = self.sequences[p_name]
            seq_c = self.sequences[c_name]
            P = GetMaxLikEstimateOfTransitionMatrix(seq_p, seq_c)            
            C_pc = GetSubstitutionCountMatrix(seq_p, seq_c)
            C_cp = GetSubstitutionCountMatrix(seq_c, seq_p)            
            logLik_pc = 0
            logLik_cp = 0
            for i in range(4):
                for j in range(4):
                    logLik_pc += C_pc[i,j]*log(P[i,j])
                    logLik_cp += C_cp[i,j]*log(P[i,j])
            self.edgeLogLikelihoods[(p_name,c_name)] = logLik_pc
            self.edgeLogLikelihoods[(c_name,p_name)] = logLik_cp
        map(ComputeLogLikelihoodsPerEdge,self.edgeLengths.keys())          
    
    def ComputeEdgeLogLikelihoodsUsingTwoDirectionalTransitionMatrices(self):
        def ComputeLogLikelihoodsPerEdge(edge):
            p_name, c_name = edge
            seq_p = self.sequences[p_name]
            seq_c = self.sequences[c_name]
            P_pc = GetMaxLikEstimateOfTransitionMatrix(seq_p, seq_c)
            P_cp = GetMaxLikEstimateOfTransitionMatrix(seq_c, seq_p)
            C_pc = GetSubstitutionCountMatrix(seq_p, seq_c)
            C_cp = GetSubstitutionCountMatrix(seq_c, seq_p)
            logLik_pc = 0
            logLik_cp = 0  
            for i in range(4):
                for j in range(4):
                    logLik_pc += C_pc[i,j]*log(P_pc[i,j])
                    logLik_cp += C_cp[i,j]*log(P_cp[i,j]) 
            self.edgeLogLikelihoods[(p_name,c_name)] = logLik_pc
            self.edgeLogLikelihoods[(c_name,p_name)] = logLik_cp
        map(ComputeLogLikelihoodsPerEdge,self.edgeLengths.keys())
    
    def EstimateSequenceAtRootGivenSequencesAtDescendants(self):
        rootSeq = ""
        root = self.root
        rootProb = self.rootProb
        Q = self.Q
        def ComputeExpectedStateForSite(site):
            maxLik = -1
            for state in range(4):
                lik_state =  rootProb[state]
                for c in root.children:
                    seq_c = self.sequences[c.name]
                    t_c = self.GetEdgeLength(root.name, c.name)
                    P_c = ComputeProbabilityMatrixUsingMatrixExponentiation(Q, t_c)
                    lik_state *=  P_c[state,DNA.index(seq_c[site])]
                if lik_state > maxLik:
                    maxLik = maxLik
                    stateWithMaxProb = DNA[state]
             
            rootSeq += stateWithMaxProb
        map(ComputeExpectedStateForSite,range(self.sequenceLength))     
        self.sequences[root.name] = rootSeq

    def OptimizeLengthsOfEdgesIncidentToRoot(self):
        root = self.root
        Q = self.Q
        seq_root = self.sequences[self.root.name]
        for c in root.children:        
            t_0 = self.GetEdgeLength(root.name, c.name)
            seq_c = self.sequences[c.name]
            def ComputeNegLogLikelihood(t):
                P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q, t)
                def ComputeNegLogLikelihoodPerSite(site):
                    log(P[seq_root[site],seq_c[site]])             
    
    def OptimizeQUsingFixedStatesAndComputeLogLik(self, maxIterations = 20, logLikImprovementThreshold = 1.0, verbose=False):                    
        self.InitializeQ(skipRoot=True)
        logLikConvergenceNotReached = True
        iteration = 1
        currentLogLik = self.logLik   
        current_rootSeq = self.sequences[self.root.name]
        while logLikConvergenceNotReached and iteration <= maxIterations:                
            # Estimate sequence at root
            self.EstimateSequenceAtRootGivenSequencesAtDescendants()
            # Initialize rootProb
            self.rootProb = GetBaseFreq(self.sequences[self.root.name])      
            # Given Q, rootProb opt length to root
                    
            # Given rootProb, Q opt length to root
        
        logLikConvergenceNotReached = True
        iteration = 1
        currentLogLik = self.logLik
        while logLikConvergenceNotReached and iteration <= maxIterations:
            iteration += 1
            logLikConvergenceNotReached = abs(self.logLik - currentLogLik) > logLikImprovementThreshold            
            currentLogLik = self.logLik

            
    def PerformEMWithFixedStatesAndEdgeLengthsAndFixedQ(self, maxIterations = 10, logLikImprovementThreshold = 1.0, verbose=False):
        logLikConvergenceNotReached = True
        iteration = 1
        currentLogLik = self.logLik          
        current_rootSeq = self.sequences[self.root.name]
        while logLikConvergenceNotReached and iteration <= maxIterations:
            self.ComputeExpectedAncestralStatesAndUpdateLogLikUsingPatternSaving()
            changeInRootSeq = hamDist(current_rootSeq, self.sequences[self.root.name])
            current_rootSeq = self.sequences[self.root.name]
            if verbose:
                print ("Expected states computed")            
                print (changeInRootSeq)
            # Maximization steps
            # Given root sequence, optimize rootProb            
            self.OptimizeRootProbability()
            if verbose:
                print ("Optimized root probability")           
            # Given rootProb and Q, optimize t
            self.OptimizeEdgeLengths()
            if verbose:
                print ("Optimized edge lengths")
            # Given t and rootProb, optimize Q
            self.OptimizeQ()
            if verbose:
                print ("Optimized Q")
                print ("iteration: "+str(iteration)+". Change in root seq is " + str(changeInRootSeq))
            iteration += 1
            logLikConvergenceNotReached = abs(self.logLik - currentLogLik) > logLikImprovementThreshold and changeInRootSeq > 1/float(self.sequenceLength)            
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

    def ContractVerticesWithOutDegree1(self):
        verticesToRemove = []
        orderedVertices = self.GetPostOrderTraversalWithoutRoot()
        for v in orderedVertices:
            if v.outDegree == 1:
                verticesToRemove.append(v)
        
        for v in verticesToRemove:
            if v.inDegree > 0:
                p = v.parent
                c = v.children[0]
                t_pv = self.GetEdgeLength(p.name, v.name)
                t_vc = self.GetEdgeLength(v.name, c.name)
                self.RemoveDirectedEdge(p.name, v.name)
                self.RemoveDirectedEdge(v.name, c.name)
                self.AddDirectedEdge(p.name, c.name,t_pv+t_vc)
                self.RemoveVertex(v.name)
            
    # use only variable sites
    def ComputeBICForProposedPruning(self,edge):
        F = RootedForest()
        for (u_name, v_name), t in self.edgeLengths.items():
            if v_name != edge[1]:
                pass
            else:
                F.AddDirectedEdge(u_name, v_name, t)
        u = F.GetVertex(u_name)
        if u.outDegree == 1:
            if u.inDegree > 0:
                p = u.parent
                c = u.children[0]
                t_pu = F.GetEdgeLength(p.name, u.name)
                t_uc = F.GetEdgeLength(u.name, c.name)
                F.RemoveDirectedEdge(p.name, u.name)                
                F.AddDirectedEdge(p.name, c.name, t_pu + t_uc)
            F.RemoveDirectedEdge(u.name, c.name)
            F.RemoveVertex(u.name)
    
    def FindEdgeForRooting(self):
        T = Tree()
        T.AddEdges(self.edgeLengths)
        T.AddSequences(self.sequences)
        def ComputeLogLikForRooting(edge):
            logLik = 0
            u_name, v_name = edge
            t = T.GetEdgeLength(u_name, v_name)
            RT = T.GetRootedTree(u_name, v_name, t/2.0, t/2.0)
            preOrderVerticesWithoutRoot = RT.GetPreOrderTraversalWithoutLeaves()[1:]            
            for u in preOrderVerticesWithoutRoot:
                for v in u.children:
                    seq_u = self.sequences[u.name]
                    seq_v = self.sequences[v.name]
                    if (u_name,v_name) in self.QForEdge:
                        Q = self.QForEdge[(u_name,v_name)]
                    else:
                        Q = self.QForEdge[(v_name,u_name)]
                    P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q, t)                    
                    logLik += sum(map(lambda site: log(P[DNA.index(seq_u[site]),DNA.index(seq_v[site])]), range(self.sequenceLength)))
            return logLik
        tup_edge_logLik = [(edge, ComputeLogLikForRooting(edge)) for edge in T.edgeLengths.keys()]
        tup_edge_logLik.sort(key= itemgetter(1),reverse = True)
        return tup_edge_logLik[0][0]
    
    def FindEdgesForRooting(self):
        T = Tree()
        T.AddEdges(self.edgeLengths)
        T.AddSequences(self.sequences)
        def ComputeLogLikForRooting(edge):
            logLik = 0
            u_name, v_name = edge
            t = T.GetEdgeLength(u_name, v_name)
            RT = T.GetRootedTree(u_name, v_name, t/2.0, t/2.0)
            preOrderVerticesWithoutRoot = RT.GetPreOrderTraversalWithoutLeaves()[1:]            
            for u in preOrderVerticesWithoutRoot:
                for v in u.children:
                    seq_u = self.sequences[u.name]
                    seq_v = self.sequences[v.name]
                    if (u_name,v_name) in self.QForEdge:
                        Q = self.QForEdge[(u_name,v_name)]
                    else:
                        Q = self.QForEdge[(v_name,u_name)]
                    P = ComputeProbabilityMatrixUsingMatrixExponentiation(Q, t)                    
                    logLik += sum(map(lambda site: log(P[DNA.index(seq_u[site]),DNA.index(seq_v[site])]), range(self.sequenceLength)))
            return logLik
        tup_edge_logLik = [(edge, ComputeLogLikForRooting(edge)) for edge in T.edgeLengths.keys()]
        tup_edge_logLik.sort(key= itemgetter(1),reverse = True)
        maxLogLik = tup_edge_logLik[0][1]
        edgesToReturn = [tup_edge_logLik[0][0]]
        for edge, logLik in tup_edge_logLik[1:]:
            if logLik == maxLogLik:
                edgesToReturn.append(edge)             
        return edgesToReturn
       
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

    def FindEdgesForRootingForGMMDeprecated(self,pool):
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
                labelAtVertex[p] += "(" + labelAtVertex[c] + ":" + str(format(self.GetEdgeLength(p.name, c.name),'0.20f'))
            else :
                labelAtVertex[p] += "," + labelAtVertex[c] + ":" + str(format(self.GetEdgeLength(p.name, c.name),'0.20f'))
            if p.outDegree - timesVisited[p] == 0:
                labelAtVertex[p] += ")"
        labelAtRoot = labelAtVertex[self.GetRoot()] + ";"
        return labelAtRoot
    
    def GetNewickLabelWithLabeledNodes(self):
        labelAtVertex = {v : "" for v in self.vertices.values() if v.outDegree > 0}
        timesVisited = {v : 0 for v in labelAtVertex.keys() if v.outDegree > 0}
        labelAtVertex.update({v : v.name for v in self.GetLeaves()})    
        postOrderTraversalWithoutRoot = self.GetPostOrderTraversalWithoutRoot()
        for c in postOrderTraversalWithoutRoot:            
            p = c.parent
            timesVisited[p] += 1
            if timesVisited[p] == 1:
                labelAtVertex[p] += "(" + labelAtVertex[c] + ":" + str(format(self.GetEdgeLength(p.name, c.name),'0.20f'))
            else :
                labelAtVertex[p] += "," + labelAtVertex[c] + ":" + str(format(self.GetEdgeLength(p.name, c.name),'0.20f'))
            if p.outDegree - timesVisited[p] == 0:
                labelAtVertex[p] += ")" + p.name
        labelAtRoot = labelAtVertex[self.GetRoot()] + ";"
        return labelAtRoot

    def GetNewickLabelDeprecated(self):
        labelAtVertex = {v : "" for v in self.vertices.values() if v.outDegree > 0}
        timesVisited = {v : 0 for v in labelAtVertex.keys() if v.outDegree > 0}
        labelAtVertex.update({v : v.name for v in self.GetLeaves()})    
        postOrderTraversalWithoutRoot = self.GetPostOrderTraversalWithoutRoot()[:-1]
        for c in postOrderTraversalWithoutRoot:
            p = c.parent
            timesVisited[p] += 1
            bran_len = f'{self.GetEdgeLength(p.name, c.name):.20f}'
            
            if timesVisited[p] == 1:
                labelAtVertex[p] += "("+labelAtVertex[c]+":"+bran_len
            else :
                labelAtVertex[p] += ","+labelAtVertex[c]+":"+bran_len
            if p.outDegree - timesVisited[p] == 0:
                labelAtVertex[p] += ")"
        labelAtRoot = labelAtVertex[self.GetRoot()] + ";"
        return labelAtRoot

    def WriteToFile(self,fileName,fileFormat="edgeList"):
        treeFile = open(fileName,'w')
        if fileFormat == "edgeList":            
            for parent_name, child_name in self.edgeLengths.keys():
                treeFile.write(parent_name+'\t'+child_name+'\t'+str(self.edgeLengths[parent_name,child_name])+'\n')            
        elif fileFormat == "newick":
#             treeFile.write(self.GetNewickLabel()+"\n")
            treeFile.write(self.GetNewickLabelWithLabeledNodes()+"\n")
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

class RootedForest(RootedTree):
    def __init__(self):
        self.vertices = {}
        self.edgeLengths = {}
        self.h_ind = 0 
        self.sequences = {}
        self.QForEdge = {}
    
    def AddVertex(self,u_name):
        self.vertices[u_name] = Vertex(u_name)
    
    def GetVertex(self,u_name):
        return self.vertices[u_name]
    
    def ContainsVertex(self,u_name):
        return u_name in self.vertices.keys()
    
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
        self.edgeLengths[(u_name,v_name)] = w
    
    def AddQForEdge(self,u_name,v_name,Q):
        self.QForEdge[(u_name,v_name)] = Q
    
    def AddEdges(self,edgeLengthsDic):
        for (u_name, v_name), t in edgeLengthsDic.items():
            self.AddEdge(u_name, v_name, t)
    
    def GetEdgeLength(self,u_name,v_name):
        return self.edgeLengths[(u_name,v_name)]               
    
    def AddRootedTree(self,RT):
        for p_name, c_name in RT.edgeLengths.keys():
            if p_name.startswith("h"):
                if int(p_name.split("hiddenVertex_")[1]) > self.h_ind:
                    self.h_ind = int(p_name.split("hiddenVertex_")[1])
            if c_name.startswith("h"):
                if int(c_name.split("hiddenVertex_")[1]) > self.h_ind:
                    self.h_ind = int(c_name.split("hiddenVertex_")[1])
            t_pc = RT.edgeLengths[p_name,c_name]
            self.AddEdge(p_name, c_name, t_pc)
        self.sequences.update(RT.sequences)
    
    def AddLeafSequences(self,sequences):
        if len(sequences) > 0:
            self.sequences.update(sequences)
            self.sequenceLength = len(sequences.values()[0])
            self.IdentifyRepeatedSitePatterns()      
    
    
    def GetSubtreeDescendingFromVertex(self, u):
        RT = RootedTree()
        parents = [u]
        RT.AddSequences(self.sequences[u.name])
        while len(parents) > 0:
            u = parents.pop()
            for c in u.children:
                t = RT.GetEdgeLength(u.name, c.name)
                RT.AddDirectedEdge(u.name, c.name, t)
                parents.append(c)
                RT.AddSequences(self.sequences[c.name])
        RT.SetRoot()
        return RT
        
    def GetComponents(self):
        roots = []
        for v in self.vertices.values():
            if v.inDegree == 0:
                roots.append(v)
        rootedTrees = [self.GetSubtreeDescendingFromVertex(r) for r in roots]
        return rootedTrees
    
    def PerformEMWithEstimatedQ(self, maxIterations = 20, logLikImprovementThreshold = 1.0, verbose=False):
        # Initialize ancestral states using MP
        self.InitializeNonLeafSequencesUsingMPUsingPatternSavings()
#         print "Initialized ancestral states"
        # Initialize rootProb, Q and t
        self.InitializeEdgeLengths()
#         print "Initialized edge lengths"
        self.InitializeQ()
#         print "Initialized Q"
        self.InitializeRootProbUsingPatternSavings()
#         print "Initialized root probability" 
        logLikConvergenceNotReached = True
        iteration = 1
        currentLogLik = -1*pow(10,10)          
        current_rootSeq = self.sequences[self.root.name]
        while logLikConvergenceNotReached and iteration <= maxIterations:        
            # Expectation step
            self.ComputeExpectedAncestralStatesAndUpdateLogLikUsingPatternSaving()
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
            # Given rootProb and Q, optimize t
            self.OptimizeEdgeLengths()
            if verbose:
                print ("Optimized edge lengths")
            # Given t and rootProb, optimize Q
            self.OptimizeQ()
            if verbose:
                print ("Optimized Q")
                print ("iteration: "+str(iteration)+". Change in root seq is " + str(changeInRootSeq))                
            iteration += 1
            logLikConvergenceNotReached = self.logLik - currentLogLik > logLikImprovementThreshold            
            currentLogLik = self.logLik
    
    def WriteToFile(self,fileName):
        edgeListFile = open(fileName,'w')
        for u_name, v_name in self.edgeLengths.keys():            
            edgeListFile.write(u_name+'\t'+v_name)
            edgeListFile.write('\t'+str(self.edgeLengths[(u_name,v_name)]))
            edgeListFile.write('\t'+str(self.QperEdge[(u_name,v_name)]+'\n'))
        edgeListFile.close()


class RootedGraph():
    def __init__(self):
        self.root=[]
        self.vertices={}
        self.vertexCount=0
        self.samplingTime={}
        self.edgeLengths={}
        self.edgeLabels={}
        self.vertexLabels={}
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
        del self.vertices[u_name]
        u.RemoveParent()
        self.vertexCount -= 1
    def AddSamplingTime(self,u_name,t_u):
        self.samplingTime[u_name]=t_u
    def AddVertexLabel(self,u_name,label):
        self.vertexLabels[u_name]=label
    def GetSamplingTime(self,u_name):
        return self.samplingTime[u_name]
    def SetRoot(self,vertexName):
        self.root = self.vertices[vertexName]
    def AddDirectedEdge(self,u_name,v_name):
        u=self.GetVertex(u_name)
        v=self.GetVertex(v_name)
        u.AddChild(v)
        v.AddParent(u)        
    def AddEdgeLength(self,u_name,v_name,weight):
        self.edgeLengths[u_name,v_name] = weight
    def AddEdgeLabel(self,u_name,v_name,label):
        self.edgeLabels[u_name,v_name] = label
    def GetEdgeLabel(self,u_name,v_name):
        return self.edgeLabels[u_name,v_name]
    def RemoveDirectedEdge(self,u_name,v_name):
        u=self.GetVertex(u_name)
        v=self.GetVertex(v_name)
        u.RemoveChild(v)
        v.RemoveParent(u)
    def GetDescendants(self,u_name):
        u=self.GetVertex(u_name)
        descendants=[]
        children=u.children
        while len(children)>0:
            v=children[0]
            descendants.append(v.name)
            del children[0]
            children.extend(v.children)
        return descendants
    def GetDistance(self,u_name,v_name):
        label_u = self.vertexLabels[u_name]
        label_v = self.vertexLabels[v_name]
        d=0
        for i in range(len(label_u)):
            if label_u[i]!=label_v[i]:
                d+=1
        return d
    def SetNumberOfDescendants(self):
        pass
    
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