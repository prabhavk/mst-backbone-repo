import numpy as np
import math as m
    
# def ComputeGTRGammaDistancesUsingRaxml(alignment,experimentName):
#     devnull=open(os.devnull,'w')
#     alignmentFileName = experimentName+'_sequences.phylip'
#     WriteAlignment(alignment, alignmentFileName, fileFormat='phylip')
#     raxmlCommandForComputingDistances = toolPath+'raxmlHPC-SSE3 -m GTRGAMMA -n\t'+experimentName+'\t-s\t'+alignmentFileName + ' -w ' + tempPath+ '\t-p 1234 -rootProb x'
#     sub.call(raxmlCommandForComputingDistances,stdout=devnull,shell=True)
#     distances = ReadDistances(tempPath+'RAxML_distances.'+experimentName)
#     sub.call('rm '+tempPath+'RAxML*.'+experimentName,stdout=devnull,shell=True)
#     devnull.close()
#     return distances

id2seq = {}

def GetPrior(alignment):
    prior=[10**-7]*4
    for seq in alignment.values():
        for nuc in seq:
            prior[int(nuc)]+=1
    for i in range(4):
        prior[i]/=float(len(alignment.values()[0]))*len(alignment)
    return prior

def GetInsertIndex(dList,d):
    n = len(dList)
        # upper and lower boundaries needed to be remembered
    i=int(m.floor(n/float(2)))
    lower = 0; upper = n-1
    while(True):
        if(dList[i] > d):
            if(i==0): return(i)
            elif(dList[i-1]>d):  upper = i-1; i=int(m.floor((lower+i)/float(2)))
            else: return(i)
        else:
            if(i==n-1): return(i+1)
            elif(dList[i+1]<=d):  lower = i+1; i=int(m.ceil((upper+i)/float(2)))
            else: return(i+1)

# Map nucleotides to indices
n2i = {}
n2i['A'] = [0]
n2i['T'] = [1]
n2i['C'] = [2]
n2i['G'] = [3]
n2i['R'] = [0,3]
n2i['Y'] = [2,1]
n2i['S'] = [3,2]
n2i['W'] = [0,1]
n2i['K'] = [3,1]
n2i['M'] = [0,2]
n2i['B'] = [1,2,3]
n2i['D'] = [0,1,3]
n2i['V'] = [0,2,3]
n2i['H'] = [0,1,2]
n2i['N'] = [0,1,2,3]

DNA = ["A","C","G","T"]

def LogDet_Steel(seq1,seq2):
    F = GetFreqMatrix(seq1,seq2)    
    f1_a, f1_c, f1_g, f1_t = GetBaseFreq(seq1)
    g1 = f1_a*f1_c*f1_g*f1_t
    f2_a, f2_c, f2_g, f2_t = GetBaseFreq(seq2)
    g2 = f2_a*f2_c*f2_g*f2_t 
    ld = -0.25*(m.log(np.linalg.det(F))-0.5*(m.log(g1*g2)))
    return(ld)
        
def GetLogDetDistances(alignment): 
    nameList = alignment.keys()
    nameList.sort()
    pairwiseDistance = {}
    for selectedOTU1 in nameList[:-1]:
        print (float(nameList.index(selectedOTU1))/float(len(nameList)))
        for selectedOTU2 in nameList[nameList.index(selectedOTU1)+1:]:
            dist = LogDet_Steel(alignment[selectedOTU1], alignment[selectedOTU2])
            pairwiseDistance[selectedOTU1,selectedOTU2] = dist
    return(pairwiseDistance)
        
def GetHammingDistances(alignment):
    nameList = alignment.keys()
    nameList.sort()
    pairwiseDistance = {}
    for selectedOTU1 in nameList[:-1]:
        for selectedOTU2 in nameList[nameList.index(selectedOTU1)+1:]:
            dist = hamDist(alignment[selectedOTU1], alignment[selectedOTU2])
            pairwiseDistance[selectedOTU1,selectedOTU2] = dist
    return(pairwiseDistance)      
        
def GetBaseFreq(seq):
    rootProb=[0.0]*4
    for s in seq:
        rootProb[DNA.index(s)]+=1.0
    return map(lambda x: x/float(len(seq)), rootProb)

def GetFreqMatrix(seq1,seq2):
    rootedTree = np.array([[0.0]*4]*4,)
    for s in range(len(seq1)):
        a = DNA.index(seq1[s])
        b = DNA.index(seq2[s])
        rootedTree[a][b]+=1/float(len(seq1))
    return(rootedTree)

def LogDetScaled(seq1,seq2):
    D = Div(seq1,seq2)
    if(np.linalg.det(D)>0):
        ld = -0.25*(m.log(np.linalg.det(D))-0.5*(m.log((D.sum(0)*D.sum(1)).prod())))
        return(round(ld*pow(10,4)))
    else:
        return(pow(10,6))

def LogDet(seq1,seq2):
    D = Div(seq1,seq2)
    logD = np.linalg.slogdet(D)[1]
    ld = -0.25*(logD-0.5*(m.log((D.sum(0)*D.sum(1)).prod())))
    if(np.linalg.det(D)>0):
        return(ld)
    else:
        return(ld)
    
def Div(seq1,seq2):
    D = np.array([[10**-5]*4]*4,np.float32)
    convertNucCode = seq1[0] not in ['0','1','2','3']
    for s in range(0,min(len(seq1),len(seq2))):
        if (seq1[s]!='-' and seq2[s]!='-'):
            if convertNucCode:
                a = n2i[seq1[s]]
                b = n2i[seq2[s]]
                d = 1/float(len(a)*len(b))
                for i in a:
                    for j in b:
                        D[i][j]+=d
            else: 
                a = int(seq1[s])
                b = int(seq2[s])
                D[a][b]+=1
    return(D)

def hamDist(seq1,seq2):
    length = len(seq1)
    d=0
    for i in range(len(seq1)):
        if seq1[i]!=seq2[i]:
            d+=1
    d/=float(length)
    return(d)

def JC69Dist(seq1,seq2):
    seqLength = len(seq1) 
    p=0
    for i in range(seqLength):
        if seq1[i]!=seq2[i]:
            p+=1
    p/=float(seqLength)    
    d = -0.75*m.log(max(10**-7,1-(4*p/float(3))))
    return(d)

def GetVarianceOfJC69Estimate(seq1,seq2):
    seqLength = len(seq1) 
    p=0
    for i in range(seqLength):
        if seq1[i]!=seq2[i]:
            p+=1
    p/=float(seqLength)    
    var = p*(1-p)/float(seqLength)
    var/=(1-4*p/float(3))**2
    return(var)

def F81Dist(seq1,seq2,seqLength,E):
    p=0
    for i in range(seqLength):
        if seq1[i]!=seq2[i]:
            p+=1
    p/=float(seqLength)
    d = -E*m.log(max(10**-7,1-p/E))
    return(d)


def GetDistances(alignment,prior=[0.25]*4,distType='JC69'):
    nameList = alignment.keys()
    nameList.sort()
    pairwiseDistance = {}
    if distType=='JC69':
        seqLength = len(alignment.values()[0])
        for selectedOTU1 in nameList[:-1]:
            for selectedOTU2 in nameList[nameList.index(selectedOTU1)+1:]:
                dist = JC69Dist(alignment[selectedOTU1], alignment[selectedOTU2])
                pairwiseDistance[selectedOTU1,selectedOTU2] = dist
    elif distType=='LogDet':
        for selectedOTU1 in nameList[:-1]:
            for selectedOTU2 in nameList[nameList.index(selectedOTU1)+1:]:
                dist = LogDet(alignment[selectedOTU1], alignment[selectedOTU2])
                pairwiseDistance[selectedOTU1,selectedOTU2] = dist
    elif distType=='F81':
        seqLength = len(alignment.values()[0])
        E = 1-(prior[0]**2+prior[1]**2+prior[2]**2+prior[3]**2)
        for selectedOTU1 in nameList[:-1]:
            for selectedOTU2 in nameList[nameList.index(selectedOTU1)+1:]:
                dist = F81Dist(alignment[selectedOTU1], alignment[selectedOTU2],seqLength,E)
                pairwiseDistance[selectedOTU1,selectedOTU2] = dist
    elif distType=='Hamming':
        for selectedOTU1 in nameList[:-1]:
            for selectedOTU2 in nameList[nameList.index(selectedOTU1)+1:]:
                dist = hamDist(alignment[selectedOTU1], alignment[selectedOTU2])
                pairwiseDistance[selectedOTU1,selectedOTU2] = dist        
    return(pairwiseDistance)

def GetDistanceCorrection(normalizedHammingDistance,E=0.75,distType='JC69'): 
    if distType == 'JC69':
        d = -0.75*m.log(max(10**-7,1-normalizedHammingDistance/0.75))
    elif distType=='F81':
        d = -E*m.log(max(10**-7,1-normalizedHammingDistance/E))
    return(d)

def GetDistanceVectorFromMatrix(distanceMatrix):
    fullVertList = [x[0] for x in distanceMatrix.keys()]
    fullVertList+= [x[1] for x in distanceMatrix.keys()]
    fullVertList = list(set(fullVertList))
    fullVertList.sort()
    vertexNameList = fullVertList[:]
    d=[]
    for i in vertexNameList:
        for j in vertexNameList[vertexNameList.index(i)+1:]:
            d.append(distanceMatrix[(i,j)])
    return d

def GetDistanceVectorFromTree(vertexNameList,T):
    d=[]
    for i in vertexNameList:
        pathsFromi = T.get_all_shortest_paths(T.vs["name"].index(i))
        for j in vertexNameList[vertexNameList.index(i)+1:]:
            path_ij = pathsFromi[T.vs["name"].index(j)]
            pathLength=0
            for pathPos in range(0,len(path_ij)-1):
                pathLength+=T.es[T.get_eid(path_ij[pathPos],path_ij[pathPos+1])]["length"]
            d.append(pathLength)
    return(d)

def GetDistanceMatrixFromTree(vertexNameList,T,idType='name'):
    T.GetRootedTree()
    # keys are index in vertexNameList
    distanceMatrix={}
    if idType =='index':
        for i in vertexNameList:
            pathsFromi = T.get_all_shortest_paths(i)
            for j in vertexNameList[vertexNameList.index(i)+1:]:
                path_ij = pathsFromi[j]
                pathLength=0
                for pathPos in range(0,len(path_ij)-1):
                    pathLength+=T.es[T.get_eid(path_ij[pathPos],path_ij[pathPos+1])]["length"]
                if i<j:
                    distanceMatrix[i,j] = pathLength
                else:
                    distanceMatrix[j,i] = pathLength
    elif idType=='name':   
        for i in vertexNameList:
            pathsFromi = T.get_all_shortest_paths(T.vs["name"].index(i))
            for j in vertexNameList[vertexNameList.index(i)+1:]:
                path_ij = pathsFromi[T.vs["name"].index(j)]
                pathLength=0
                for pathPos in range(0,len(path_ij)-1):
                    pathLength+=T.es[T.get_eid(path_ij[pathPos],path_ij[pathPos+1])]["length"]
                if i<j:
                    distanceMatrix[i,j] = pathLength
                else:
                    distanceMatrix[j,i] = pathLength
    return(distanceMatrix)

def GetDistanceVector(distanceMatrix):
    vertexNameList = [x[0] for x in distanceMatrix.keys()]
    vertexNameList+= [x[1] for x in distanceMatrix.keys()]
    vertexNameList = list(set(vertexNameList))
    vertexNameList.sort()
    d=[]
    for vert_1 in vertexNameList:
        for vert_2 in vertexNameList[vertexNameList.index(vert_1)+1:]:
            d.append(distanceMatrix[vert_1,vert_2])
    return d

# Get a dictionary of key (observed vertex pair) and value (set of indices where a mismatch occurs)
def GetMismatchIndices(alignment):
    vertexNameList = alignment.keys()
    vertexNameList.sort()
    mismatches={}
    seqLength = len(alignment.values()[0])
    for vert1 in vertexNameList:
        seq1=alignment[vert1]
        for vert2 in vertexNameList[vertexNameList.index(vert1)+1:]:
            indexList=[]
            seq2=alignment[vert2]
            for i in range(seqLength):
                if seq1[i]!=seq2[i]:
                    indexList.append(i)
            mismatches[(vert1,vert2)] = set(indexList)
    return (mismatches)
            
def ComputeDistanceFromAlignmentAndWriteToFile(alignment,fileName,distType='JC',symmetric=False):
    rootProb=open(fileName,'w')
    nameList = alignment.keys()
    nameList.sort()
    # whole distance matrix is written 
    if symmetric:
        if distType == 'JC':
            for i in range(len(nameList)):
                for j in range(len(nameList)):
                    dist = JC69Dist(alignment[nameList[i]], alignment[nameList[j]])                
                    rootProb.write(str(dist)+'\t')
                rootProb.write('\n')
        else:
            for i in range(len(nameList)):
                for j in range(i+1,len(nameList)):
                    dist = LogDet(alignment[nameList[i]], alignment[nameList[j]])
                    rootProb.write(nameList[i]+'\t'+nameList[j]+'\t'+ str(dist)+'\n')
    else:
        # variable name pair and distance per line
        if distType == 'JC':
            for i in range(len(nameList)):
                for j in range(i+1,len(nameList)):
                    dist = JC69Dist(alignment[nameList[i]], alignment[nameList[j]])                
                    rootProb.write(nameList[i]+'\t'+nameList[j]+'\t'+ str(dist)+'\n')
        else:
            for i in range(len(nameList)):
                for j in range(i+1,len(nameList)):
                    dist = LogDet(alignment[nameList[i]], alignment[nameList[j]])
                    rootProb.write(nameList[i]+'\t'+nameList[j]+'\t'+ str(dist)+'\n')
            rootProb.close()

def WriteFullDistanceMatrixToFile(distanceMatrix,fileName):
    rootProb=open(fileName,'w')
    N = int((1+m.sqrt(1+8*len(distanceMatrix)))/float(2))
    vertexNameList = ['O' + str(x+1) for x in range(N)]
    vertexNameList.sort()
    # Create full distanceMatrix
    for vert1 in vertexNameList:
        distanceMatrix[vert1,vert1]=0
        for vert2 in vertexNameList[vertexNameList.index(vert1)+1:]:
            distanceMatrix[vert2,vert1] = distanceMatrix[vert1,vert2]
    # Write distances
    for vert1 in vertexNameList:
        for vert2 in vertexNameList:
            rootProb.write(str(distanceMatrix[vert1,vert2])+'\t')
        rootProb.write('\n')
    rootProb.close()     

def GetSeqIndex(N,part,totalNumberOfParts):
    if part==totalNumberOfParts:
        return N
    else:
        k = N-0.5-m.sqrt(4*N*(N-1)*(1-part/float(totalNumberOfParts))+1)/2.0
        return int(m.ceil(k))

# def ComputeDistanceMatrixInParts(alignmentFileName,distanceDirectoryPath,distanceFunction,part,totalNumberOfParts):
#     alignment = ReadAlignment(alignmentFileName)
#     N =  len(alignment)
#     startIndex = GetSeqIndex(N,part-1,totalNumberOfParts)
#     endIndex = GetSeqIndex(N,part,totalNumberOfParts)
#     seqNames = alignment.keys()
#     seqNames.sort()
#     freq = ComputeBaseFrequencies(alignment)
#     distanceFile = open(distanceDirectoryPath + 'distances_'+str(part),'w')
#     for seqIndex1 in range(startIndex,endIndex):
#         for seqIndex2 in range(seqIndex1+1,N):
#             dist = distanceFunction(alignment[seqNames[seqIndex1]], alignment[seqNames[seqIndex2]],freq)
#             distanceFile.write(seqNames[seqIndex1]+'\t'+seqNames[seqIndex2]+'\t'+str(dist)+'\n')
#     distanceFile.close()

# TN93 distance
def ComputeBaseFrequencies(alignment):
    freq = {}
    nucCode=['A','G','T','C']
    for nuc in nucCode:
        freq[nuc]=0
    numberOfUnambiguousBases=0
    numOfSeqs=min(1000,len(alignment))
    for i in range(numOfSeqs):
        for char in alignment.values()[i]:
            numberOfUnambiguousBases+=1
            if char in nucCode:
                freq[char]+=1
    freq['A']/=float(numberOfUnambiguousBases)
    freq['G']/=float(numberOfUnambiguousBases)
    freq['T']/=float(numberOfUnambiguousBases)
    freq['C']/=float(numberOfUnambiguousBases)
    freq['R']=freq['A']+freq['G']
    freq['Y']=freq['T']+freq['C']
    return freq

def Compute_A2G_and_T2C_and_TranversionFreq(seq1,seq2):
    seqLength=float(len(seq1))
    a2gTransitions=0
    t2cTransitions=0
    transversions=0
    for i in range(len(seq1)):
        if seq1[i]!=seq2[i] and seq1[i]!='-' and seq2[i]!='-':
            if (seq1[i]=='A' and seq2[i] =='G') or (seq2[i]=='A' and seq1[i] =='G'):
                a2gTransitions+=1
            elif (seq1[i]=='T' and seq2[i] =='C') or (seq2[i]=='T' and seq1[i] =='C'):
                t2cTransitions+=1
            elif seq1[i] in ['A','G','T','C'] and seq2[i] in ['A','G','T','C']:
                transversions+=1
    a2gTransitions/=seqLength
    t2cTransitions/=seqLength
    transversions/=seqLength
    return (a2gTransitions,t2cTransitions,transversions)

def ComputeTN93GammaDistance(seq1,seq2,freq,gamma=1):
    gamma = float(gamma)
    k1=2*freq['A']*freq['G']/freq['R']
    k2=2*freq['T']*freq['C']/freq['Y']
    k3=2*((freq['R']*freq['Y'])-(freq['A']*freq['G']*freq['Y']/freq['R'])-(freq['T']*freq['C']*freq['R']/freq['Y']))
    k4=2*((freq['A']*freq['G'])+(freq['T']*freq['C'])+(freq['R']*freq['Y']))
    p1,p2,q = Compute_A2G_and_T2C_and_TranversionFreq(seq1,seq2)
    w1 = 1-p1/k1-q/(2*freq['R'])
    w2 = 1-p2/k2-q/(2*freq['Y'])
    w3 = 1-(q/(2*freq['R']*freq['Y']))
    d=gamma*(k1*(w1**(-1/gamma))+k2*(w2**(-1/gamma))+k3*(w3**(-1/gamma))-k4)
    return d


def ComputeTN93Distance(seq1,seq2,freq):
    k1=2*freq['A']*freq['G']/freq['R']
    k2=2*freq['T']*freq['C']/freq['Y']
    k3=2*((freq['R']*freq['Y'])-(freq['A']*freq['G']*freq['Y']/freq['R'])-(freq['T']*freq['C']*freq['R']/freq['Y']))
    p1,p2,q = Compute_A2G_and_T2C_and_TranversionFreq(seq1,seq2)
    w1 = 1-p1/k1-q/(2*freq['R'])
    w2 = 1-p2/k2-q/(2*freq['Y'])
    w3 = 1-(q/(2*freq['R']*freq['Y']))
    d=-1*(k1*m.log(w1)+k2*m.log(w2)+k3*m.log(w3))
    return d

