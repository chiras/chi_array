########################################################################################
######## LIBRARY MANAGEMENT ############################################################
########################################################################################

textlevel = leveltext("Checking R packages","up",textlevel)

if (sum(installed.packages()[,1]=="DBI")==0) {install.packages("DBI", dependencies = TRUE, quiet = TRUE)}
if (sum(installed.packages()[,1]=="RSQLite")==0) {install.packages("RSQLite", dependencies = TRUE, quiet = TRUE)}
if (sum(installed.packages()[,1]=="xtable")==0) {install.packages("xtable", dependencies = TRUE, quiet = TRUE)}
if (sum(installed.packages()[,1]=="XML")==0) {install.packages("XML", dependencies = TRUE, quiet = TRUE)}
if (sum(installed.packages()[,1]=="made4")==0) {install.packages("made4", dependencies=T, quiet = TRUE)}

textlevel = leveltext("Looking for Bioconductor","keep",textlevel)
source("http://bioconductor.org/biocLite.R")

textlevel = leveltext("Checking Bioconductor packages","keep",textlevel)
#biocLite (pkgs=c("DBI","affy","vsn","genefilter","multtest","sigPathway","affyPLM","hgu133plus2.db","KEGG.db","limma","hgu133plus2.db","annotate","RBGL","GOstats","convert","KEGGgraph", "AnnotationDbi","gplots","RSQLite","xtable","XML","SPIA","RColorBrewer","RBGL","grid","fUtilities"),ask=F, quiet = TRUE)

if (sum(installed.packages()[,1]=="affy")==0) {biocLite("affy")}
if (sum(installed.packages()[,1]=="vsn")==0) {biocLite("vsn")}
if (sum(installed.packages()[,1]=="genefilter")==0) {biocLite("genefilter")}
if (sum(installed.packages()[,1]=="multtest")==0) {biocLite("multtest")}
if (sum(installed.packages()[,1]=="sigPathway")==0) {biocLite("sigPathway")}
if (sum(installed.packages()[,1]=="affyPLM")==0) {biocLite("affyPLM")}
if (sum(installed.packages()[,1]=="hgu133plus2.db")==0) {biocLite("hgu133plus2.db")}
if (sum(installed.packages()[,1]=="KEGG.db")==0) {biocLite("KEGG.db")}
if (sum(installed.packages()[,1]=="limma")==0) {biocLite("limma")}
if (sum(installed.packages()[,1]=="annotate")==0) {biocLite("annotate")}
if (sum(installed.packages()[,1]=="RBGL")==0) {biocLite("RBGL")}
if (sum(installed.packages()[,1]=="GOstats")==0) {biocLite("GOstats")}
if (sum(installed.packages()[,1]=="convert")==0) {biocLite("convert")}
if (sum(installed.packages()[,1]=="KEGGgraph")==0) {biocLite("KEGGgraph")}
if (sum(installed.packages()[,1]=="AnnotationDbi")==0) {biocLite("AnnotationDbi")}
if (sum(installed.packages()[,1]=="gplots")==0) {biocLite("gplots")}
if (sum(installed.packages()[,1]=="Rgraphviz")==0) {biocLite("Rgraphviz")}
if (sum(installed.packages()[,1]=="Rgraphviz")==0) {stop("Please make sure you also install Graphviz and rerun this script afterwards. \n Check http://www.graphviz.org/\n\n***Script Aborted***")}
if (sum(installed.packages()[,1]=="SPIA")==0) {biocLite("SPIA")}
if (sum(installed.packages()[,1]=="RColorBrewer")==0) {biocLite("RColorBrewer")}
if (sum(installed.packages()[,1]=="grid")==0) {biocLite("grid")}
#if (sum(installed.packages()[,1]=="fUtilities")==0) {biocLite("fUtilities")}

textlevel = leveltext("Loading packages","keep",textlevel)

library(DBI)
library(AnnotationDbi)
library(annotate)
library(affy)				
library(vsn)
library(genefilter)
library(multtest) 
library(sigPathway)
library(affyPLM)
library(KEGG.db)
library(limma)
library(hgu133plus2.db)
library(RBGL)
library(GOstats)
library(convert)
library(gplots)
#library(Rgraphviz)
#library(KEGGgraph)

library(SPIA)
library(RColorBrewer)
library(grid)
#library(fUtilities)

textlevel = leveltext("Loading chip annotation data","keep",textlevel)
load(paste(PATHress,"/GenesetsU133plus2.RData",sep=""))												# Chip Annotation data


textlevel = leveltext("Preparing annotation data\n","keep",textlevel)
pathinfo <- as.list(hgu133plus2PATH2PROBE)
pregeneinfo <-  hgu133plus2SYMBOL
mapped_probes <- mappedkeys(pregeneinfo)
geneinfo <- as.list(pregeneinfo[mapped_probes])


#end
