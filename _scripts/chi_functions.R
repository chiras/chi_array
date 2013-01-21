########################################################################################
######## LIMMA ANALYSIS ################################################################
########################################################################################
chi_limma <- function(eset,design,contrast.matrix,name){

cat(noquote(paste("Limma Analysis (",name,")\n",sep="")))
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit3 <- fit2

cat(noquote("Limma Analysis: Mapping Annotation Data\n"))
fit2$genes$Symbol=getSYMBOL(fit2$genes$ID, "hgu133plus2")
fit2$genes$GeneName <- unlist(mget(fit$genes$ID, hgu133plus2GENENAME))
fit2$genes$Path <- paste(mget(fit$genes$ID, hgu133plus2PATH))
fit2$genes$EG <- getEG(fit2$genes$ID, "hgu133plus2")

cat(noquote("Limma Analysis: Calculating Statistics\n"))
#F.stat <- FStat(fit3)
#p.value <- pf(F.stat,df1=attr(F.stat,"df1"),df2=attr(F.stat,"df2"),lower.tail=FALSE)
topnumber = length(featureNames(eset))
tab <- topTable(fit2, adjust = "fdr", number=topnumber) # FDR corrected p-values nach dem cut off
#head(tab)



cat(noquote("Limma Analysis: Writing results to disk\n"))
write.table(tab, file=paste('R-Output/results_',name,'_top',topnumber,'-',normalisierung,'.txt', sep=""),sep=";")
analysis_log=c(analysis_log,"Limma Analysis (Single)")

dir_log=c(dir_log,paste('R-Output/results_top',topnumber,'-',normalisierung,'.txt', sep=""))
file_log=c(file_log,"Limma Data")

return(tab)
}

########################################################################################
######## FILTERING DATA ################################################################
########################################################################################



chi_filter <- function(eset,tab,IQR_filter,logFC_filter,AFFX_filter,IQR_intens_above,IQR_intens_probes,IQR_greater_than,logFC_threshold){

myfilter = ""
filtered = tab$ID
number = length(tab$ID)

if (logFC_filter){
	cat(noquote(paste("Filtering: logFC Filter: ", number, "->",sep="")))

	logFC_f1 <- tab[tab$logFC > logFC_threshold | tab$logFC < -logFC_threshold,1]
	logFC_f2 <- logFC_f1  %in% featureNames(eset)
	filterOne = logFC_f1[logFC_f2]
	filtered = intersect(filterOne,filtered)
	myfilter = paste(myfilter,"_logFC",sep ="") 
	analysis_log=c(analysis_log,"Filtering by logFC")
	number = length(filtered)
	cat(noquote(paste(number," probesets\n",sep="")))
}


if (AFFX_filter){
	cat(noquote(paste("Filtering: AFFX Filter: ", number, "->",sep="")))

	affx_probes=grep("AFFX",tab$ID)

	filterTwo <- tab$ID[tab$adj.P.Val<min(tab$adj.P.Val[affx_probes])]
	filtered = intersect(filterTwo,filtered)

	myfilter = paste(myfilter,"_AFFX",sep ="") 
	analysis_log=c(analysis_log,"Filtering by AFFX")
	number = length(filtered)
	cat(noquote(paste(number," probesets\n",sep="")))

}

### OPTION 2:  Filtern nach IQR log ratio
if (IQR_filter){
	cat(noquote(paste("Filtering: IQR Filter: ", number, "->",sep="")))
	IQR_f1 			<- pOverA(IQR_intens_probes, log2(IQR_intens_above)) 				
	IQR_f2 			<- function(x) (IQR(x) > IQR_greater_than) 		
	ff <- filterfun(IQR_f1,IQR_f2) ### Filterfunction
	filterThree <- genefilter(eset, ff) 
	filterThree  = names(filterThree)[filterThree]
	filtered = intersect(filterThree,filtered)

	myfilter = paste(myfilter,"_IQR",sep ="") 
	analysis_log=c(analysis_log,"Filtering by IQR")
	number = length(filtered)
	cat(noquote(paste(number," probesets\n",sep="")))
}

eset <- eset[filtered, ] 
return(eset)
}

########################################################################################
######## Venn Preparation ##############################################################
########################################################################################

chi_eset2sigs <- function (eset, tab.orig){
	cat(noquote(paste("Filtering: Defining significant rows in tabs\n",sep="")))

	inters= !tab.orig$ID %in% featureNames(eset)
	#tab.tmp = 2*(tab.orig$t>0)-1 NOT VERIFIED
	#tab.tmp[inters] = 0
	return(tab.orig[!inters,])
}

chi_eset2bins <- function (eset, tab.orig){
	cat(noquote(paste("Filtering: Defining significant rows in tabs\n",sep="")))

	inters= !tab.orig$ID %in% featureNames(eset)
	tab.tmp = 2*(tab.orig$t>0)-1 
	tab.tmp[inters] = 0
	return(tab.tmp)
}


chi_venn <- function(setlistPLUS,setlistMINUS){
if (length(setlistPLUS)>1){

OLlistPLUS <- overLapper(setlist=setlistPLUS, sep="_", type="vennsets")
OLlistMINUS <- overLapper(setlist=setlistMINUS, sep="_", type="vennsets")

pdf(file=paste('R-Output/MFA-Plots/Venn_',paste(names(setlistPLUS), collapse = "-"),'.pdf', sep=""))
	counts <- list(sapply(OLlistPLUS$Venn_List, length), sapply(OLlistMINUS$Venn_List, length));
	vennPlot(counts=counts, mysub="Top: Up-regulated; Bottom: Down-regulated")
dev.off()

}}

########################################################################################
######## Heat Preparation #(OLD)########################################################
########################################################################################

	
Probesets2Heatmap <- function(eset.norm, VectorOfProbesets){ 
eset2plot <- eset.norm[VectorOfProbesets, ] # Expressionsset für diese Probes generieren

i=1
while (file.exists(paste("R-Output/Paths/path_ind",i,".pdf",sep=""))){
	i=i+1
}
pdf(file=paste("R-Output/Paths/path_ind",i,".pdf",sep=""),  height=length(featureNames(eset2plot))/2+2,width=8)
	
heatmap.2(exprs(eset2plot), labRow=getSYMBOL(featureNames(eset2plot), "hgu133plus2"), Rowv = F, Colv = F,dendrogram="none",lmat=rbind( c(0, 4,0,0), c(2,1,1,0 ), c(0,3,0,0) ), lwid=c(0.5, 2, 4,0.5 ), lhei= c(2,length(featureNames(eset2plot))/2,0.1), density.info="none", scale="row", key = T, symkey = T, trace ="none")
	 # # hier wird die eigentliche Heatmap erstellt!
mtext("           Individual Heatmap Selection", side=2, line=-2, adj=0.0, cex=1, col="black", outer=TRUE)
dev.off()

}

PathSigHeat <- function(eset_sig, pathwaynr){
	
cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Obtaining pathway associated probesets: ", sep="")))
eset.sig.norm.path <- eset.sig.norm[grep(pathwaynr,kegg)]

symbs<- mget(featureNames(eset.sig.norm.path), hgu133plus2SYMBOL, ifnotfound = NA)

symbs<- unlist(symbs)

used_symbs=c()
cat(noquote(paste(length(symbs),"/",length(featureNames(eset.sig.norm.path)), " of ",length(kegg),"\n", sep="")))

cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Getting Gene Symbols\n", sep="")))
for (i in 1:length(symbs)){
	xxnum=0;
	while (sum(used_symbs == paste(symbs[i], xxnum, sep="-"))){
		xxnum= xxnum +1;
	}
	symbs[i]=paste(symbs[i], xxnum, sep="-")	
	used_symbs = c(used_symbs, symbs[i])
}

head(symbs)

cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Defining new Matrix\n", sep="")))
mat <- exprs(eset.sig.norm.path)
cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Renaming Matrix\n",sep="")))
rownames(mat) <- ifelse(!is.na(symbs), as.vector(symbs), featureNames(eset.sig.norm.path))

cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Plotting to PDF\n", sep="")))
pdf(file=paste('R-Output/Paths/path_top',topnumber,'-KEGG', pathwaynr,'-',auswertung,'-',normalisierung,'-',myfilter,'.pdf', sep=""), height=length(rownames(mat))/2+2,width=8)

heatmap.2(mat[order(rownames(mat)), ] , Rowv = F, Colv = F,dendrogram="none", lmat=rbind( c(0, 4,0,0), c(2,1,1,0 ), c(0,3,0,0) ), lwid=c(0.5, 2, 4,0.5 ), lhei= c(2,length(rownames(mat))/2,0.1),density.info="none", scale="row", key = T, symkey = T, trace ="none")
mtext(paste("           ",as.character(mget(pathwaynr, KEGGPATHID2NAME))), side=2, line=-2, adj=0.0, cex=1, col="black", outer=TRUE)

dev.off()

}

InteractKeggPlot <- function(top_table, pathName){ 

	top= top_table
	pathName = pathName
	x <- hgu133plus2ENTREZID
	top$ENTREZ <- unlist(as.list(x[top$ID]))
	top <- top[!is.na(top$ENTREZ),] 
	top <- top[!duplicated(top$ENTREZ),] 
	tg1 <- top #[top$adj.P.Val < 0.05,]
	DE_tg1 <- tg1$logFC
	names(DE_tg1) <- as.vector(tg1$ENTREZ)	
	ALL_tg1 <- top$ENTREZ

	tmp <- paste("../KEGGdata/hsa", pathName,".xml",sep="")
	
	if (!file.exists(tmp)){download.file(paste("http://www.genome.jp/kegg-bin/download?entry=hsa", pathName,"&format=kgml",sep=""), tmp)}
	
	pdf(file=paste('R-Output/Paths/path_inter',topnumber,'-KEGG', pathName,'-',auswertung,'-',normalisierung,'-',myfilter,'.pdf', sep=""), width=20, height=20)
	g <- parseKGML2Graph(tmp)
	
	deKID <- translateGeneID2KEGGID(names(DE_tg1))
	allKID <- translateGeneID2KEGGID(ALL_tg1)
	
	isDiffExp <- nodes(g) %in% deKID
	# sprintf("%2.2f%% genes differentially-expressed", mean(isDiffExp)*100)
	
	g.old = g
	g <- subGraph(nodes(g)[isDiffExp], g)
	
	ar <- 20
	cols <- heatPalette(ar) # colorRampPalette(brewer.pal(6, "RdBu"))(ar)
	colsO <- heatPalette(2) # colorRampPalette(brewer.pal(6, "RdBu"))(ar)
	
	logfcs <- DE_tg1[match(nodes(g), deKID)]
	names(logfcs) <- nodes(g)
	logfcs[is.na(logfcs)] <- 0
#	incol <- round((logfcs+2)*5);
	incol <- round((logfcs+min(logfcs)*-1)*5);

	outcol <-logfcs
	outcol[outcol>0]<-2
	outcol[outcol<0]<-1
	incol[incol>ar] <- ar
	undetected <- !nodes(g) %in% allKID
		
	logcol <- cols[incol]; 
	logcol[logfcs==0] <- "darkgrey"; 
	logcol[undetected] <- "yellow";
	
	logcolO <- colsO[outcol]; 
	
	names(logcol) <- names(logfcs)
	nA <- makeNodeAttrs(g, fillcolor=logcol, color= logcolO, label="", width=10, height=1.2,shape="triangle")
	par(mar=c(3,5,0,5), mgp=c(0,0,0))
	
	gGeneID <- translateKEGGID2GeneID(nodes(g))
	gSymbol <-  sapply(gGeneID, function(x) mget(x, org.Hs.egSYMBOL, ifnotfound=NA)[[1]])
		
	#gnA$width <- makeAttr(g,list("0.8"=toprbccName))
	
	layout(mat=matrix(c(rep(1,8),2), ncol=1, byrow=TRUE))
	gnA <- makeNodeAttrs(g, fillcolor=logcolO, lwd=3, label=gSymbol, fontsize=1,fixedsize=T)
	plot(g, "twopi", nodeAttrs=gnA)
	dev.off()
}

chi_gostats <- function(analysisname, eset.sig.norm,eset.norm){

cat(noquote("GOstats: Obtaining Gene Clusters and Universe\n"))
geneCluster=as.vector(na.omit(getEG(featureNames(eset.sig.norm),"hgu133plus2")))
geneCluster2=as.vector(na.omit(getSYMBOL(featureNames(eset.sig.norm),"hgu133plus2")))
geneCluster3=featureNames(eset.sig.norm)

geneUniverse=as.vector(na.omit(getEG(featureNames(eset.norm),"hgu133plus2"))) 
geneUniverse2=as.vector(na.omit(getSYMBOL(featureNames(eset.norm),"hgu133plus2"))) 
geneUniverse3=featureNames(eset.norm)

gos = c("KEGG","MF","BP","CC")
gosigs=list()

topnumber=length(featureNames(eset.sig.norm))

cat(noquote("GOstats: Creating Outfile: "))
cat(noquote(paste('R-Output/GO_',analysisname,'_','top',topnumber,'-PATH-',auswertung,'-',normalisierung,'.txt\n', sep="")))


outfile = paste('R-Output/GO_',analysisname,'_','top',topnumber,'-PATH-',normalisierung,'.txt', sep="")

dir_log=c(dir_log, outfile)
file_log=c(file_log,"GoStats Results")

sink(outfile)
writeLines(paste("","GOid","Type","Name","Pvalue","OddsRatio","ExpCount","Count","Size","GOGenesNumber","SigGenesNumber","SigGenesLinCorr","GOGenesSymbol","GOGenesEG","SigGenesSymbol","SigGenesEG",sep=";"), con=paste('R-Output/GO_',analysisname,'_','top',topnumber,'-PATH-',normalisierung,'.txt', sep=""))	
sink()


for (id in 1:4){
	cat(noquote(paste("GOstats: Using Over Data for: ",gos[id],"\n", sep="")))
	cat(noquote("GOstats: Calculation\n"))
	if (gos[id] != "KEGG"){
		paramsGOunder <- new("GOHyperGParams", geneIds = geneCluster, universeGeneIds = geneUniverse, annotation = "hgu133plus2", ontology = gos[id], pvalueCutoff = 1, conditional = FALSE, testDirection = "under")
		paramsGOover <- new("GOHyperGParams", geneIds = geneCluster, universeGeneIds = geneUniverse, annotation = "hgu133plus2", ontology = gos[id], pvalueCutoff = 1, conditional = FALSE, testDirection = "over")
		tryCatch(hgOverGO <- hyperGTest(paramsGOover),error = function(e) {print('error GO')})
		#tryCatch(hgUnderGO <- hyperGTest(paramsGOunder),error = function(e) {print('error GO')})
		#htmlReport(hgUnderGO,file=paste('GOunder_',gos[id],'_top',topnumber,'-PATH-',auswertung,'-',normalisierung,'.html', sep=""))
		#htmlReport(hgOverGO,file=paste('GOover_',gos[id],'_top',topnumber,'-PATH-',auswertung,'-',normalisierung,'.html', sep=""))
	}else{
		paramsGOunder <- new("KEGGHyperGParams", geneIds = geneCluster, universeGeneIds = geneUniverse, annotation = "hgu133plus2",pvalueCutoff = 1, testDirection = "under")
		paramsGOover <- new("KEGGHyperGParams", geneIds = geneCluster, universeGeneIds = geneUniverse, annotation = "hgu133plus2",pvalueCutoff = 1, testDirection = "over")
		tryCatch(hgOverGO <- hyperGTest(paramsGOover),error = function(e) {print('error KEGG')})
		#tryCatch(hgUnderKEGG <- hyperGTest(paramsKEGGunder),error = function(e) {print('error KEGG')})
		#htmlReport(hgOverKEGG,file=paste('KEGGover_top',topnumber,'-PATH-',auswertung,'-',normalisierung,'-',myfilter,'.html', sep=""))
		#htmlReport(hgUnderKEGG,file=paste('KEGGunder_top',topnumber,'-PATH-',auswertung,'-',normalisierung,'-',myfilter,'.html', sep=""))
	}
	cat(noquote("GOstats: Summarizing\n"))
	gosigs[[gos[id]]]=summary(hgOverGO)[summary(hgOverGO)$Pvalue<0.05,][[1]]
	output = gos[id]

	listofGOs = gosigs[[output]]

if (length(listofGOs)>0){

for (i in 1:length(listofGOs)){ 
	
analyzeGO=F;
if (gos[id] == "KEGG" & length(unlist(mget(listofGOs[i],hgu133plus2PATH2PROBE, ifnotfound=NA))[is.na(unlist(mget(listofGOs[i],hgu133plus2PATH2PROBE, ifnotfound=NA)))==F])>0){analyzeGO<-T}
if (gos[id] != "KEGG" & length(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE  , ifnotfound=NA))[is.na(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE  , ifnotfound=NA)))==F])>0){analyzeGO<-T}

if(analyzeGO){
genefunction=list();

cat(noquote(paste("GOstats: Working on GO-Set: ",gos[id] ," - ", i, " (",length(listofGOs),") ", sep="")))
cat(noquote("[1] Collecting Data"))
myinfo = summary(hgOverGO)[summary(hgOverGO)[,1]==listofGOs[i],]

cat(noquote(" [2] Parsing Data"))
genefunction$GOid					= myinfo[,1]	
genefunction$Type					= names(myinfo)[1]	
genefunction$Name 					= myinfo$Term
genefunction$Pvalue 				= myinfo$Pvalue
genefunction$OddsRatio 				= myinfo$OddsRatio
genefunction$ExpCount 				= myinfo$ExpCount
genefunction$Count 					= myinfo$Count
genefunction$Size 					= myinfo$Size

if(gos[id] != "KEGG"){
	genefunction$GOGenesNumber 			= length(unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F))))
	tmpGOGenesSymbol 					= unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F)))
	tmpGOGenesEG 						= unique(as.vector(unlist(getEG(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F)))
	#genefunction$AllGenesChipNumber 	= length(unique(intersect(tmpGOGenesSymbol, geneUniverse2)))
	genefunction$SigGenesNumber 		= length(unique(intersect(tmpGOGenesSymbol, geneCluster2)))
	genefunction$SigGenesCorrected 		= (genefunction$Size/genefunction$GOGenesNumber)*genefunction$SigGenesNumber
	genefunction$GOGenesSymbol 			= paste(unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F))), collapse=", ")
	genefunction$GOGenesEG 				= paste(unique(as.vector(unlist(getEG(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F))), collapse=", ")
	#genefunction$AllGenesChipSymbol 	= paste(unique(intersect(tmpGOGenesSymbol, geneUniverse2)), collapse=", ")
	#genefunction$AllGenesChipEG 		= paste(unique(intersect(tmpGOGenesEG, geneUniverse)), collapse=", ")
	genefunction$SigGenesSymbol 		= paste(unique(intersect(tmpGOGenesSymbol, geneCluster2)), collapse=", ")
	genefunction$SigGenesEG				= paste(unique(intersect(tmpGOGenesEG, geneCluster)), collapse=", ")
}else{
	genefunction$GOGenesNumber 			= length(unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F))))
	tmpGOGenesSymbol 					= unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F)))
	tmpGOGenesEG 						= unique(as.vector(unlist(getEG(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F)))
	#genefunction$AllGenesChipNumber 	= length(unique(intersect(tmpGOGenesSymbol, geneUniverse2)))
	genefunction$SigGenesNumber 		= length(unique(intersect(tmpGOGenesSymbol, geneCluster2)))
	genefunction$SigGenesCorrected 		= (genefunction$Size/genefunction$GOGenesNumber)*genefunction$SigGenesNumber
	genefunction$GOGenesSymbol 			= paste(unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F))), collapse=", ")
	genefunction$GOGenesEG 				= paste(unique(as.vector(unlist(getEG(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F))), collapse=", ")
	#genefunction$AllGenesChipSymbol 	= paste(unique(intersect(tmpGOGenesSymbol, geneUniverse2)), collapse=", ")
	#genefunction$AllGenesChipEG 		= paste(unique(intersect(tmpGOGenesEG, geneUniverse)), collapse=", ")
	genefunction$SigGenesSymbol 		= paste(unique(intersect(tmpGOGenesSymbol, geneCluster2)), collapse=", ")
	genefunction$SigGenesEG				= paste(unique(intersect(tmpGOGenesEG, geneCluster)), collapse=", ")	
}

cat(noquote(" [3] Writing Data\n"))
cat(noquote(paste(c("",t(as.data.frame(genefunction))),collapse=";")),file=outfile, append=T, sep="\n")

}else{
	cat(noquote(paste("GOstats: GO-Set excluded: ",gos[id] ," - ", i, " (",length(listofGOs),") \n", sep="")))
}# end if (analyzeGO)
}# end for length(listofGOs))
}# end if (length(listofGOs)>0)

cat(noquote(" Making Barplot\n"))
cat(noquote(paste(c("",t(as.data.frame(genefunction))),collapse=";")),file=outfile, append=T, sep="\n")

sig.hgOverGO=summary(hgOverGO)[summary(hgOverGO)$Pvalue<0.1,]

sig.hgOverGO =sig.hgOverGO[order(sig.hgOverGO$Pvalue, decreasing=T),]

pdf(file=paste('R-Output/GO_',analysisname,'_',gos[id],'_top',topnumber,'-PATH-',normalisierung,'.pdf', sep=""), height=4+length(sig.hgOverGO$Count)/10,width=10)
 	par(mar=c(3,20,0.1,0.2))
	barplot(rbind(sig.hgOverGO$ExpCount,sig.hgOverGO$Count), col=c("gray80","gray20"), beside=T, horiz=T, names.arg=paste(sig.hgOverGO$KEGGID,sig.hgOverGO$Term,"- p =", round(sig.hgOverGO$Pvalue, digits=4)), las=2, cex.names=0.5)


 	for (i in 1:(max(sig.hgOverGO$Count)/5)){
 		abline(v=i*5, lty=3)
 	}
 	par(xpd=NA)
 	legend(-15,0,legend=c("Count", "Expected"), fill=c("gray20","gray80"))
dev.off()


}# end for MF/BP/CC/KEGG
}# end function analyzeGOstats
