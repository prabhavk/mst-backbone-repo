#!/usr/bin/env Rscript
library(ape)
args <- commandArgs(TRUE)
treeFileName = args[1]
edgeListFileName = args[2]
tree <- read.tree(treeFileName)
experimentName = paste(sort(tree$tip.label)[1:3],collapse='-')
# Assumes that tree is leaf-labeled
numberOfLatentVertices = Nnode(tree)
latentVertexNames <- paste('h_',experimentName,c(1:numberOfLatentVertices),sep='')
vertexIndexToName <- c(1:(Ntip(tree)+Nnode(tree)))
names(vertexIndexToName)<-c(tree$tip.label,latentVertexNames)
names(vertexIndexToName)[tree$edge[,1]]
names(vertexIndexToName)[tree$edge[,2]]
edgeListToWrite <- paste(names(vertexIndexToName)[tree$edge[,1]],names(vertexIndexToName)[tree$edge[,2]],tree$edge.length,sep='\t')
write.table(x=edgeListToWrite,file=edgeListFileName,quote=F,col.names=F,row.names=F)
quit(save="no")

