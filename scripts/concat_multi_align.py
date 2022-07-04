from fileIO import ReadAlignment

referenceSequence = ReadAlignment(projectPath+'data/gisaid/influenza_'+reference+'_HA_reference.fas') # modify

if concatenateMSAFiles:
    alignment = ReadAlignment(projectPath+'data/gisaid/'+lineage+'/HA_sequences_unaligned.fasta')  # modify
    # align sequences in batches of 200 sequences
    noOfSeqsInBatch = 0    
    numberOfParts = 0
    seqsRemaining = len(alignment)
    multipleSequenceAlignment = {}
    for seqId, seq in alignment.iteritems(): 
        noOfSeqsInBatch += 1
        if noOfSeqsInBatch == 200:
            numberOfParts += 1
            noOfSeqsInBatch = 0            
            seqsRemaining -= 200
    
    for part in range(1,numberOfParts+1):
        alignedSequences = ReadAlignment(projectPath+'data/gisaid/'+lineage+'/aligned_sequences/HA_sequences_aligned_part_'+str(part)+'.fasta')
        print (len(alignedSequences))
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
        print ('part', part, 'of', numberOfParts)
    WriteAlignment(multipleSequenceAlignment, projectPath+'data/gisaid/'+lineage+'/aligned_seqs.fas', 'fasta')