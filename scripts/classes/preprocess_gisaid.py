# get counts for each segment
from fileIO import ReadAlignment, WriteAlignment
import pandas as pd
lineages = ['H1N1_2009_pandemic','H1N1_seasonal','H3N2']
segmentCount = {'PB2': 0,'PB1': 0,'PA': 0,'HA': 0,'NP': 0,'NA': 0,'MP': 0,'NS': 0,'HE': 0,'P3': 0}

from os import walk
from config import projectPath, toolPath, scriptPath
excelFileNames = []

createBatchFastaFileForAligning = False
concatenateMSAFiles = True
filterSequences = False
stripAllInfoExceptEPIIds = False

def ComputeSequenceCoverageForAllPos(alignment):
    numberOfSequences = len(alignment)
    sequenceLength= len(alignment.values()[0])
    sequenceCoverage = [0.0]*sequenceLength
    for seq in alignment.values():
        for pos in xrange(sequenceLength):
            if seq[pos] in ['A','C','G','T']:
                sequenceCoverage[pos] += 1.0/float(numberOfSequences)    
    return sequenceCoverage

def GetTrimmedSeqAndSeqQuality(seq,posToKeep):
    trimmedSeq = ""
    for pos in posToKeep:
        trimmedSeq += seq[pos]
    numberOfNonAmbiguousChars= sum(map(lambda char: 1 if char in ['A','C','G','T'] else 0,trimmedSeq))
    return trimmedSeq, numberOfNonAmbiguousChars/float(len(trimmedSeq))


# lineage = "H1N1_2009_pandemic"
# reference = "H1N1"

lineage = "H3N2"
reference = "H3N2"
referenceSequence = ReadAlignment(projectPath+'data/gisaid/influenza_'+reference+'_HA_reference.fas')

if createBatchFastaFileForAligning:
#     HA_alignment = {}
#     for (_, _, fileNames) in walk(projectPath+'data/gisaid/'+lineage+'/'):
#         for f in fileNames:
#             isolateInformation={}
#             if f.endswith('.xls'):
#                 excelFileNames.append(projectPath+'data/gisaid/'+lineage+'/'+f)
#     fileNo = 1
#     for excelFileName in excelFileNames:
#         print 'processing file', fileNo, 'of', len(excelFileNames)
#         fileNo += 1
#         fileHeader = excelFileName.split('.xls')[0]
#         fastaFileName = fileHeader.split('isolates')[0]+'sequence'+fileHeader.split('isolates')[1]+'.fasta'
#         xlFile = pd.read_excel(excelFileName, sheetname='Tabelle1')
#         isolate_id_list = xlFile['Isolate_Id']
#         host_list = xlFile['Host']    
#         location_list = xlFile['Location']
#         collection_date_list = xlFile['Collection_Date']
#         for i in xrange(len(host_list)):
#             if location_list.isnull()[i]:
#                 location_list[i] = ""
#             if host_list[i] == 'Human':
#                 isolateInformation[isolate_id_list[i].encode('utf-8')]={'location':location_list[i].encode('utf-8'), 'collection_date':collection_date_list[i].encode('utf-8')}
#         alignment = ReadAlignment(fastaFileName)
#         for seqId in alignment.keys():
#             segment = seqId.split('| ')[1].strip()
#             seq = alignment[seqId]
#             if segment == 'HA':
#                 isolateId = seqId.split(' |')[0].strip()
#                 if isolateId in isolateInformation.keys():
#                     newSeqId = isolateId+';'+lineage+';'+isolateInformation[isolateId]['location']+';'+isolateInformation[isolateId]['collection_date']
#                     HA_alignment[newSeqId] = seq.upper()
    
    HA_alignment = ReadAlignment(projectPath+'data/gisaid/'+lineage+'/HA_sequences_unaligned.fasta')
    # align sequences in batches of 200 sequences
    noOfSeqsInBatch = 0
    smallAlignment = {}
    numberOfParts = 0
    seqsRemaining = len(HA_alignment)
    for seqId, seq in HA_alignment.iteritems():
        smallAlignment[seqId] = seq 
        noOfSeqsInBatch += 1
        if noOfSeqsInBatch == 200:
            numberOfParts+= 1
            noOfSeqsInBatch = 0
            smallAlignment.update(referenceSequence)
            WriteAlignment(smallAlignment, projectPath+'data/gisaid/'+lineage+'/unaligned_sequences/HA_sequences_unaligned_part_'+str(numberOfParts)+'.fasta')
            smallAlignment = {}
            seqsRemaining -= 200
            print 'seqs remaining is', seqsRemaining
    
    numberOfParts += 1
    smallAlignment.update(referenceSequence)
    WriteAlignment(smallAlignment, projectPath+'data/gisaid/'+lineage+'/unaligned_sequences/HA_sequences_unaligned_part_'+str(numberOfParts)+'.fasta')
    
    mafftFile = open(projectPath+'scripts/batch_mafft.sh','w')
    for part in range(1,numberOfParts+1):
        unaligned_sequences_batch_fileName = projectPath+'data/gisaid/'+lineage+'/unaligned_sequences/HA_sequences_unaligned_part_'+str(part)+'.fasta'
        aligned_sequences_batch_fileName = projectPath+'data/gisaid/'+lineage+'/aligned_sequences/HA_sequences_aligned_part_'+str(part)+'.fasta'
        mafftCommand = toolPath+'mafft_bin/mafft --maxiterate 1000 --localpair  --thread 5 --quiet ' + unaligned_sequences_batch_fileName + ' > ' + aligned_sequences_batch_fileName 
        mafftFile.write(mafftCommand+'\n')
    
    
    mafftFile.close()
    
if concatenateMSAFiles:
    HA_alignment = ReadAlignment(projectPath+'data/gisaid/'+lineage+'/HA_sequences_unaligned.fasta')
    # align sequences in batches of 200 sequences
    noOfSeqsInBatch = 0    
    numberOfParts = 0
    seqsRemaining = len(HA_alignment)
    multipleSequenceAlignment = {}
    for seqId, seq in HA_alignment.iteritems(): 
        noOfSeqsInBatch += 1
        if noOfSeqsInBatch == 200:
            numberOfParts += 1
            noOfSeqsInBatch = 0            
            seqsRemaining -= 200
    
    for part in range(1,numberOfParts+1):
        alignedSequences = ReadAlignment(projectPath+'data/gisaid/'+lineage+'/aligned_sequences/HA_sequences_aligned_part_'+str(part)+'.fasta')
        print len(alignedSequences)
        noOfPosInAlignment = len(alignedSequences.values()[0])
        posToKeep = [True]*noOfPosInAlignment
        alignedRefSeq = alignedSequences[referenceSequence.keys()[0]]
        for pos in xrange(noOfPosInAlignment):
            if alignedRefSeq[pos]=='-':
                posToKeep[pos] = False
        for seqid in alignedSequences.keys():
            if seqid != referenceSequence.keys()[0]:
                trimmedSeq = ""
                for pos in xrange(noOfPosInAlignment):
                    if posToKeep[pos]:
                        trimmedSeq += alignedSequences[seqid][pos]
                multipleSequenceAlignment[seqid] = trimmedSeq.upper()
        print 'part', part, 'of', numberOfParts
    WriteAlignment(multipleSequenceAlignment, projectPath+'data/gisaid/'+lineage+'/aligned_seqs.fas', 'fasta')

if filterSequences:
    alignment = ReadAlignment(projectPath+'data/gisaid/'+lineage+'/aligned_seqs.fas')
    sequenceCoverage = ComputeSequenceCoverageForAllPos(alignment)
    posToKeep = []
    for pos in xrange(len(sequenceCoverage)):
        if sequenceCoverage[pos] >= 0.8:
            posToKeep.append(pos)
    
    print "Length of trimmed alignment is", len(posToKeep)        
    numberOfSeqsThatPassQualityCheck = 0
    seqIdAndCollectionTimesTuple = []
    for seqId in alignment.keys():
        seq = alignment[seqId]
        trimmedSeq, seqQuality = GetTrimmedSeqAndSeqQuality(seq, posToKeep)
        if seqQuality == 1.0:
            numberOfSeqsThatPassQualityCheck += 1
            seqIdAndCollectionTimesTuple.append((seqId,seqId.split(';')[3]))
    seqIdAndCollectionTimesTuple.sort(key=lambda x: x[1],reverse=False)
    
    print 'Total number of sequences is', len(alignment)
    print 'Number of sequences that pass quality check are', numberOfSeqsThatPassQualityCheck
    print 'Sorted seqIds of interest wrt collection time'
    uniqueSeqs = set([])
    alignmentToAnalyze = {}
    for seqId, collectionTime in seqIdAndCollectionTimesTuple:
        seq = alignment[seqId]
        trimmedSeq, seqQuality = GetTrimmedSeqAndSeqQuality(seq, posToKeep)
        if trimmedSeq not in uniqueSeqs:
            uniqueSeqs.update(set([trimmedSeq]))
            alignmentToAnalyze[seqId] = trimmedSeq
    print 'Number of distinct sequences that pass quality check are', len(alignmentToAnalyze)
    WriteAlignment(alignmentToAnalyze, projectPath+'data/gisaid/'+lineage+'/H3N2_HA_seqsWithLineageCollectionTimeAndLocationInFastaHeader.fasta')
#     print 'Alignment to use has been written to file'
if stripAllInfoExceptEPIIds:
    allInfoAsATabSeperatedFile = open(projectPath+'data/gisaid/'+lineage+'/'+lineage+'_dataUsed','w')
    allInfoAsATabSeperatedFile.write('EPI_id\tLocation\tTime\tSequence\n')
    alignment = ReadAlignment(projectPath+'data/gisaid/'+lineage+'/'+lineage+'_HA_seqsWithLineageCollectionTimeAndLocationInFastaHeader.fasta')
    newAlignment = {}
    for seqIdWithLocationAndTimeInfo, seq in alignment.iteritems():
        seqId, lineage, location, time = seqIdWithLocationAndTimeInfo.split(';')
        newAlignment[seqId] = seq
        allInfoAsATabSeperatedFile.write(seqId+'\t'+location+'\t'+time+'\t'+seq+'\n')
    allInfoAsATabSeperatedFile.close()
    WriteAlignment(newAlignment, projectPath+'data/gisaid/'+lineage+'/'+lineage+'.fasta')      
    
#     qualityCheckBatchFile = open(scriptPath+'batchQualityCheck.sh','w')
#     alignmentFileName = projectPath+'data/gisaid/H1N1_2009_pandemic/aligned_seqs.fas'
#     for minSeqQuality in range(80,101):
#         minSeqQuality /= float(100)
#         for minSeqCoverage in range(80,101):
#             minSeqCoverage /= float(100)
#             qualityCheckBatchFile.write('python '+scriptPath+'qualityCheck.py\t')
#             qualityCheckBatchFile.write(alignmentFileName+'\t'+str(minSeqCoverage)+'\t'+str(minSeqQuality)+'\n')           
#                  
#     qualityCheckBatchFile.close()
# summarizeQualityCheckFile = open(projectPath+'data/gisaid/H1N1_2009_pandemic/summarizedQualityCheck','w')
# for minSeqQuality in range(80,101):
#     minSeqQuality /= float(100)
#     for minSeqCoverage in range(80,101):
#         minSeqCoverage /= float(100)
#         numberOfUniqueSequencesFile = open(projectPath+'data/gisaid/H1N1_2009_pandemic/qualityCheck/numUniqSeqs_minSeqCoverage_'+str(minSeqCoverage)+'_minSeqQuality_'+str(minSeqQuality),'r')
#         numberOfUniqueSequences = numberOfUniqueSequencesFile.readline().strip()
#         summarizeQualityCheckFile.write(numberOfUniqueSequences+'\t'+str(minSeqCoverage)+'\t'+str(minSeqQuality)+'\n')  
#  
# summarizeQualityCheckFile.close()

# segmentCountFile = open(projectPath+'data/gisaid/segmentCount','w')
# for segment, count in segmentCount.iteritems():
#     segmentCountFile.write(str(segment)+'\t'+str(count)+'\n')
# segmentCountFile.close()
# /TL/euresist_phylodynamics/work/Tools/mafft_bin/mafft --maxiterate 1000 --localpair  --thread 5 --quiet /TL/euresist_phylodynamics/work/Projects/MSTBasedForests/data/gisaid/H1N1_2009_pandemic/unaligned_sequences/HA_sequences_unaligned_part_7.fasta > /TL/euresist_phylodynamics/work/Projects/MSTBasedForests/data/gisaid/H1N1_2009_pandemic/aligned_sequences/HA_sequences_aligned_part_37.fasta

# 96,910 HA sequences. 

# align HA sequences. 