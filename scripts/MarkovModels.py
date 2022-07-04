# from phylogeneticMethods import GetNJTree
# projectPath= "/local/home/pk/Projects/mstBasedPhylogenetics/results/edgeReliabilityForFixedTrees/"
    
import math as m
from scipy.linalg import expm, logm, eig
from scipy.optimize import minimize
from numpy.linalg import inv, eig, det, lstsq
from numpy.random import uniform, choice
from numpy import diag, array, matrix, append
import numpy as np
from math import exp, log
from decimal import Decimal
import cmath as cm
DNA=["A","C","G","T"]


def GenerateProbabilityDistribution():
    p = map(lambda x: uniform(size=1)[0],range(4))
    p /= sum(p)
    return p

def GenerateProbabilityDistributionWithMinVal(p_min):
    continueTrying = True
    while (continueTrying): 
        p = list(map(lambda x: uniform(size=1)[0],range(4)))
        p /= sum(p)
        if (min(p) > p_min):
            continueTrying = False
    return p

def GetStationaryDistribution(Q):
    b = array([1,0,0,0,0])
    b.shape=(5,1)
    Q_aug = np.c_[np.ones(4),Q]
    Q_aug = Q_aug.transpose()
    pi = lstsq(Q_aug,b,rcond=-1)[0]
    return (pi)
    
def Get11FreeRates(Q):
    rates_11 = []    
    for i in range(4):
        for j in range(4):
            if (i != j):                
                if (len(rates_11) < 11):
                    rates_11.append(Q[i,j]/Q[3,2])
    return (rates_11)

def GetCalibratedRateMatrixFrom11FreeRates(rates_11):
    Q = array([[0.0]*4]*4)    
    a,b,c,d,e,f,g,h,i,j,k = rates_11
    l = 1
    D1 = -(a+b+c)
    D2 = -(d+e+f)
    D3 = -(g+h+i)
    D4 = -(j+k+l)
    Q[0,] = [D1, a, b, c]
    Q[1,] = [d, D2, e, f]
    Q[2,] = [g, h, D3, i]
    Q[3,] = [j, k, l, D4]
    Q_cal = NormalizeQMatrix(Q)
    return (Q_cal)

# GTR rate matrix
def GenerateQ_GTR(stationaryDistribution,rates):
    fa, fc, fg, ft = stationaryDistribution
    alpha, beta, gamma, delta, epsilon, eta = rates
    Q_GTR = array([[-((fc*alpha)+(fg*beta)+(ft*gamma)),(fc*alpha),(fg*beta),(ft*gamma)],
                   [(fa*alpha),-((fa*alpha)+(fg*delta)+(ft*epsilon)),(fg*delta),(ft*epsilon)],
                   [(fa*beta),(fc*delta),-((fa*beta)+(fc*delta)+(ft*eta)),(ft*eta)],
                   [(fa*gamma),(fc*epsilon),(fg*eta),-((fa*gamma)+(fc*epsilon)+(fg*eta))]])
    mu = -1*((fa*Q_GTR[0,0])+(fc*Q_GTR[1,1])+(fg*Q_GTR[2,2])+(ft*Q_GTR[3,3]))
    Q_GTR /= mu
    return(Q_GTR)

# 12 parameter rate matrix 
# [A->C,A->G,A->T,C->A,C->G,C->T,G->A,G->C,G->T,T->A,T->C,T->G]
# [A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G]

def GenerateQ(rates):
    A2C, A2G, A2T, C2A, C2G, C2T, G2A, G2C, G2T, T2A, T2C, T2G = rates
    Q = array([[-(A2C+A2G+A2T),A2C,A2G,A2T],
                   [C2A,-(C2A+C2G+C2T),C2G,C2T],
                   [G2A,G2C,-(G2A+G2C+G2T),G2T],
                   [T2A,T2C,T2G,-(T2A+T2C+T2G)]])
    rootProb = GetStationaryDistribution(Q)
    mu = -1*((rootProb[0]*Q[0,0])+(rootProb[1]*Q[1,1])+(rootProb[2]*Q[2,2])+(rootProb[3]*Q[3,3]))
    Q/=mu
    return(Q)

def GenerateQForStationaryDistribution(stationary_distribution):
    pi = array([0.0]*4)
    for i in range(4):
        pi[i] = stationary_distribution[i]    
    Q = array([[0.0]*4]*4)
    params = [1]*11

    def ConstructRateMatrixFromParams(params):
        a,b,c,d,e,f,g,h,i,j,k = params
        l = 1
        D1 = -(a+b+c)
        D2 = -(d+e+f)
        D3 = -(g+h+i)
        D4 = -(j+k+l)
        Q[0,] = [D1, a, b, c]
        Q[1,] = [d, D2, e, f]
        Q[2,] = [g, h, D3, i]
        Q[3,] = [j, k, l, D4]
        return (Q)

    def penalty(params):
        Q = ConstructRateMatrixFromParams(params)
        product = pi.dot(Q)
        penalty = sum(map(abs,product))
        if (min(params) < 0):
            penalty *= 10
        return (penalty)

    res = minimize(penalty,params,method='Nelder-Mead')
    assert(min(res.x)>0)
    Q_res = ConstructRateMatrixFromParams(res.x)
    pi_s = GetStationaryDistribution(Q_res)
    # print(pi_s)
    Q_norm = NormalizeQMatrix(Q_res)
    return (Q_norm)

def NormalizeQMatrix(Q):
    b = array([1,0,0,0,0])
    b.shape=(5,1)
    Q_aug = np.c_[np.ones(4),Q]
    Q_aug = Q_aug.transpose()
    eq_dist = lstsq(Q_aug,b,rcond=-1)[0]
    mu = 0
    for i in range(4):
        mu -= eq_dist[i]*Q[i,i]
    Q_normalized = Q[:]
    for i in range(4):
        for j in range(4):
            Q_normalized[i,j]/=mu
    return Q_normalized

def ComputeQMatrixFromVectorOfOffDiagonalElements(Q_vector_offDiag):
#     if min(Q_vector_offDiag) < pow(10,-4):
#         for i in range(12):
#             Q_vector_offDiag[i] = pow(10,-4)             
    Q = array([[0.0]*4]*4)
    Q[0,1] = Q_vector_offDiag[0]
    Q[0,2] = Q_vector_offDiag[1]
    Q[0,3] = Q_vector_offDiag[2]
    Q[1,0] = Q_vector_offDiag[3]
    Q[1,2] = Q_vector_offDiag[4]
    Q[1,3] = Q_vector_offDiag[5]
    Q[2,0] = Q_vector_offDiag[6]
    Q[2,1] = Q_vector_offDiag[7]
    Q[2,3] = Q_vector_offDiag[8]
    Q[3,0] = Q_vector_offDiag[9]
    Q[3,1] = Q_vector_offDiag[10]
    Q[3,2] = Q_vector_offDiag[11]
    for i in range(4):
        Q[i,i] = -sum(Q[i,:])
    return Q

def ComputeQVectorOfOffDiagonalElementsFromQMatrix(Q):
    Q_vector_offDiag = [0.0]*12
    Q_vector_offDiag[0] = Q[0,1]
    Q_vector_offDiag[1] = Q[0,2]
    Q_vector_offDiag[2] = Q[0,3]
    Q_vector_offDiag[3] = Q[1,0]
    Q_vector_offDiag[4] = Q[1,2]
    Q_vector_offDiag[5] = Q[1,3]
    Q_vector_offDiag[6] = Q[2,0]
    Q_vector_offDiag[7] = Q[2,1]
    Q_vector_offDiag[8] = Q[2,3]
    Q_vector_offDiag[9] = Q[3,0]
    Q_vector_offDiag[10] = Q[3,1]
    Q_vector_offDiag[11] = Q[3,2]    
    return Q_vector_offDiag


def ComputeProbabilityMatrixFromEigenVectorsAndEigenValuesIncludingImagPart(eigenValuesOfQ, eigenVectorsOfQ, t):   
    scaledEigenValues = eigenValuesOfQ*t
    expD = diag(map(cm.exp,scaledEigenValues))
    P = eigenVectorsOfQ.dot(expD).dot(inv(eigenVectorsOfQ))
    return P   

def ComputeProbabilityMatrixUsingMatrixExponentiation(Q, t): 
    P = expm(Q*t)
    return P   


def ComputeProbabilityMatrixFromEigenVectorsAndEigenValues(eigenValuesOfQ, eigenVectorsOfQ, t):   
    scaledEigenValues = eigenValuesOfQ*t
    expD = diag(map(cm.exp,scaledEigenValues))
    P = eigenVectorsOfQ.dot(expD).dot(inv(eigenVectorsOfQ))
    return P.real   

def GenerateEvolveCharFunction(P):
    
    def EvolveChar(char):
        row = DNA.index(char)
        P_row = P[row,:]
        cumSum = 0
        u = uniform()
        for col in range(4):
            cumSum += P_row[col]
            if cumSum > u:
                break
        return DNA[col]
    
    return EvolveChar


def ComputeProbabilityMatrix(Q,t):
    eigenValuesOfQ, eigenVectorsOfQ = eig(Q)
    scaledEigenValues = eigenValuesOfQ*t
    expD = diag(map(cm.exp,scaledEigenValues))
    P = eigenVectorsOfQ.dot(expD).dot(inv(eigenVectorsOfQ))
    return P

def SampleFromDistribution(rootProb,length):    
    seq=""    
    for _ in range(length):
        u = uniform()
        cumSum=0
        for base in range(4):
            cumSum += rootProb[base]
            if cumSum > u:
                break
        seq+=DNA[base]
    return seq

def GenerateRootSequence(seqLength):
    p = GenerateProbabilityDistribution()
    return SampleFromDistribution(p, seqLength)


def GetRandomlySampledRootSequenceAndRootDist(sequenceLength):
    rootProb = uniform(size=4)
    rootProb /= sum(rootProb)
    seq = SampleFromDistribution(rootProb, sequenceLength)
    return seq, rootProb

def GenerateRandomQ():
    Q = array([[0.0]*4]*4)
    for i in range(4):
        for j in range(4):
            if i != j:
                Q[i,j] = uniform(size=1)
        Q[i,i] = -sum(Q[i,:])    
        
    return NormalizeQMatrix(Q)

def GenerateQForBaseFreq_via_optimization_depr(stationaryDist):
    Q = array([[0.0]*4]*4)
    pi_1 = stationaryDist[0]
    pi_2 = stationaryDist[1]
    pi_3 = stationaryDist[2]
    pi_4 = 1 - pi_1 - pi_2 - pi_3
    a, b, c, d, e, f, g, h = uniform(0,1,8)
    i = (1-(pi_1*(a+b+2*c)+pi_2*(d+e+2*f)+pi_3*(g+h)))/(2*pi_3)
    print ("i", i)
    j = (pi_1*(a+b+c)-pi_2*d-pi_3*g)/pi_4
    print ("j", j)
    k = (pi_2*(d+e+f)-pi_1*a-pi_3*h)/pi_4
    print ("k", k)
    l = (1+pi_3*(g+h)-pi_1*(a+2*c+3*b)-pi_2*(d+2*f+3*e))/(2*pi_4)
    print ("l", l)
    if i < 0:
        print ("i is less than zero")
    if j < 0:
        print ("j is less than zero")
    if k < 0:
        print ("k is less than zero")
    if l < 0:
        print ("l is less than zero")
    D1 = -(a+b+c)
    D2 = -(d+e+f)
    D3 = -(g+h+i)
    D4 = -(j+k+l)
    Q[0,] = [D1, a, b, c]
    Q[1,] = [d, D2, e, f]
    Q[2,] = [g, h, D3, i]
    Q[3,] = [j, k, l, D4]
         
    return Q

def GenerateQForBaseFreq_algebraically_depr(stationaryDist):
    Q = array([[0.0]*4]*4)
    pi_1 = stationaryDist[0]
    pi_2 = stationaryDist[1]
    pi_3 = stationaryDist[2]
    pi_4 = 1 - pi_1 - pi_2 - pi_3
    a, b, c, d, e, f, g, h = uniform(0,1,8)
    i = (1-(pi_1*(a+b+2*c)+pi_2*(d+e+2*f)+pi_3*(g+h)))/(2*pi_3)
    print ("i", i)
    j = (pi_1*(a+b+c)-pi_2*d-pi_3*g)/pi_4
    print ("j", j)
    k = (pi_2*(d+e+f)-pi_1*a-pi_3*h)/pi_4
    print ("k", k)
    l = (1+pi_3*(g+h)-pi_1*(a+2*c+3*b)-pi_2*(d+2*f+3*e))/(2*pi_4)
    print ("l", l)
    if i < 0:
        print ("i is less than zero")
    if j < 0:
        print ("j is less than zero")
    if k < 0:
        print ("k is less than zero")
    if l < 0:
        print ("l is less than zero")
    D1 = -(a+b+c)
    D2 = -(d+e+f)
    D3 = -(g+h+i)
    D4 = -(j+k+l)
    Q[0,] = [D1, a, b, c]
    Q[1,] = [d, D2, e, f]
    Q[2,] = [g, h, D3, i]
    Q[3,] = [j, k, l, D4]
         
    return Q

def GenerateRandomQWithAribitrayPrecisionArithmetic():
    Q = array([[0.0]*4]*4)
    for i in range(4):
        for j in range(4):
            if i != j:
                Q[i,j] = Decimal(uniform(size=1)[0])
        Q[i,i] = -sum(Q[i,:])    
        
    return NormalizeQMatrix(Q)



# Generates a rate matrix by computing the matrix logarithm of a randomly set 
# transition matrix such that the largest entry in each row is on the leading diagonal.  
def GenerateQForTransitionMatrixWithLargeDiagonal():    
    P = GenerateTransitionMatrix()
    Q = logm(P).real    
    return NormalizeQMatrix(Q)

# Generates a transition matrix of three types
# unrestricted: each element is independently sampled from a uniform distribution and scaled such that each row sums to one
# slow: Additionally, the largest element in each row is placed on the leading diagonal. 
# fast: Additionally, the smallest element in each row is placed on the leading diagonal.

def GenerateTransitionMatrix(modelType = "unrestricted"):
    P = array([[0.0]*4]*4)
    for row in range(4):
        p = uniform(size = 4)
        p /= sum(p)        
        if modelType != "unrestricted":
            if modelType == "slow":
                p_diag = max(p)
            elif modelType == "fast":
                p_diag = min(p)
            for ind in range(4):
                if p[ind] == p_diag:
                    swap_index = ind
            p[swap_index] = p[row]
            p[row] = p_diag    
        P[row,:] = p[:]
    return P

# Generates a transition matrix such that
# each diagonal element is sampled from unif [p_min, 1.0)

def GenerateRestrictedTransitionMatrix(p_min = 0.8):
    P = array([[0.0]*4]*4)
    for row in range(4):        
        p = [0.0]*4
        p[row] = uniform(low = p_min, high = 1.0, size = 1)[0]
        p_rest = uniform(low = 0, high = 1.0, size = 3)
        sf = (1 - p[row]) / sum(p_rest)
        p_rest *= sf
        ind = 0
        for i in range(4):
            if i!=row:
                p[i] = p_rest[ind]
                ind += 1                        
        P[row,:] = p[:]
    return P

def GetTimeDurationForTransitionMatrix(P):    
    Q = NormalizeQMatrix(logm(P).real)
    tr_Q = Q[0,0] + Q[1,1] + Q[2,2] + Q[3,3]
    try:
        t = log(det(P))/tr_Q
    except:
        print ('#### P #####')
        print (P)
        t = 1.0
    return t

def GetMaxLikEstimateOfTransitionMatrix(seq_parent,seq_child):
    P = array([[pow(10,-5)]*4]*4)
    for site in range(len(seq_parent)):
        P[DNA.index(seq_parent[site]),DNA.index(seq_child[site])] += 1.0
    for row in range(4):
        rowSum = sum(P[row,:]) 
        for col in range(4):
            P[row,col] /= rowSum
    return P

def GetSubstitutionCountMatrix(seq_parent,seq_child):
    C = array([[0]*4]*4)
    for site in range(len(seq_parent)):
        C[DNA.index(seq_parent[site]),DNA.index(seq_child[site])] += 1
    return C

def EstimateQFromFullyLabeledRootedTree(RT,sequences):
    P = array([[0.0]*4]*4)    
    root = RT.GetRoot()
    unvisitedVertices = root.children
    sequenceLength = len(sequences.values()[0])
    while len(unvisitedVertices) > 0:
        c = unvisitedVertices.pop()
        p = c.parent
        seq_c = sequences[c.name]
        seq_p = sequences[p.name]
        for site in range(sequenceLength):
            P[DNA.index(seq_p[site]),DNA.index(seq_c[site])]+=1
        unvisitedVertices.extend(c.children)
    for i in range(4):
        P[i,:] /= sum(P[i,:])
    Q = NormalizeQMatrix(logm(P).real)    
    return(Q)

# Jukes-Cantor model
a = float(1)/float(3)
Q_JC = array([[-1.0, a, a, a], [a, -1.0, a, a], [a, a, -1.0, a], [a, a, a, -1.0]])

# General time-reversible model with parameters as specified in Waddell and Steel (1997). 
rootProb = [0.3037, 0.3079, 0.1313, 0.2571]
rates = [0.0042, 0.0466, 0.0016, 0.0002, 0.0963, 0.0007]
Q_GTR = GenerateQ_GTR(rootProb, rates)

# Non-reversible models

# Taq polymerase rates as specified in Sanson et al (2002).

Q_taq = NormalizeQMatrix(array([[-0.3595,0.0065,0.2876,0.0654],[0.0065,-0.1046,0.0,0.0980],[0.1242,0.0065,-0.1503,0.0196],[0.0915,0.2810,0.0131,-0.3856]]))

# Cumulative substitutions as specified in Randall et al (2016)
# projectPath = "/home/pk/Projects/MSTBasedForests/data/experimentalPhylogeny/Randall2016/"
# allSequences = ReadAlignment(projectPath+"allDNA.fas")
# RT = ReadRootedTree(projectPath+"RootedExperimentalPhylogeny")
# Q_Randall = EstimateQFromFullyLabeledRootedTree(RT, allSequences)
 
# A2T = float(65)
# A2C = float(16)
# A2G = float(95)
# T2A = float(68)
# T2C = float(96)
# T2G = float(19)
# C2T = float(137)
# C2G = float(15)
# C2A = float(29)
# G2C = float(21)
# G2T = float(39)
# G2A = float(126)
# 
# rates = [A2C,A2G,A2T,C2A,C2G,C2T,G2A,G2C,G2T,T2A,T2C,T2G]
# Q_Randall = GenerateQ(rates)

