#!/usr/bin/env Rscript
if (!require("phytools")) install.packages("phytools")
library(phytools)
if (!require("ape")) install.packages("ape")
library(ape)
if (!require("phangorn")) install.packages("phangorn")
library(phangorn)

argv <- commandArgs(TRUE)


fastRFS_tree<- argv[1]
path_treefile<- argv[2]
#tree <- ape::read.tree("/global/scratch/m.desousaviolante/tmv_mondial/pgSNP/file_fastRFS.single")
tree <- ape::read.tree(fastRFS_tree)


tree$edge.length<-NULL
tree=unroot(tree)
write.tree(tree,paste(path_treefile,'fastRFS_clean',sep=""))

list_files=list.files(path=path_treefile,pattern="*.treefile")
for (i in 1:length(list_files)){
  t=paste0(path_treefile,list_files[i])
  tree<-read.tree(t)
  PatristicDistMatrix<-cophenetic(tree)
  out=paste0(t,"matrix_ok")
  writeDist(PatristicDistMatrix, out, format="phylip")
}
