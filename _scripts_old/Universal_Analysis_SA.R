
########################################################################################
######## LIMMA ANALYSIS ################################################################
########################################################################################
chi_limma <- function(eset,design,contrast.matrix){

cat(noquote("Limma Analysis\n"))
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
write.table(tab, file=paste('R-Output/results_top',topnumber,'-',normalisierung,'.txt', sep=""),sep=";")
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
######## LIMMA AFTER FILTER ############################################################
########################################################################################

cat(noquote("Limma Analysis: Recalculation after Filtering\n"))
fit <- lmFit(eset.sig.norm, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

cat(noquote("Limma Analysis: Recalculation after Filtering: Mapping Annotation Data\n"))
fit2$genes$Symbol=getSYMBOL(fit2$genes$ID, "hgu133plus2")
fit2$genes$GeneName <- unlist(mget(fit$genes$ID, hgu133plus2GENENAME))
fit2$genes$Path <- paste(mget(fit$genes$ID, hgu133plus2PATH))
fit2$genes$EG <- getEG(fit2$genes$ID, "hgu133plus2")

cat(noquote("Limma Analysis: Recalculation after Filtering: Getting Top Hits\n"))
topnumber = length(featureNames(eset.sig.norm))
tab <- topTable(fit2, adjust = "fdr", number=topnumber) # FDR corrected p-values nach dem cut off
#head(tab)

cat(noquote("Limma Analysis: Writing results to disk\n"))
write.table(tab, file=paste('R-Output/results_top',topnumber,'-',normalisierung,myfilter,'.txt', sep=""),sep=";")

analysis_log=c(analysis_log,"Limma Analysis II")
dir_log=c(dir_log,paste('R-Output/results_top',topnumber,'-',normalisierung,myfilter,'.txt', sep=""))
file_log=c(file_log,"Limma Data II")

########################################################################################
######## PATHWAY ANALYSIS ##############################################################
########################################################################################

if(parseGenes2Pathways){
analysis_log=c(analysis_log,"Parsing Pathways")

cat(noquote("Pathway Parsing\n"))
 #Pathway-Zuordnung
 tab2 = tab[1,]
 tab2[1,] = c(rep("NA",length(tab2[1,])))
 rownames(tab2) = "Start"

 for (i in 1:(length(tab$Path))){
 	print(i);
 	if (is.na(tab$Path[i]) | tab$Path[i]  == "NA"){
 		tab2 =rbind(tab2,tab[i,])
 	}else{
 		pathxways= gregexpr("[0-9]+", tab$Path[i])
 		for (j in 1:length(attr(pathxways[[1]],"match.length"))){
 			bla=substr(tab$Path[i],pathxways[[1]][j],pathxways[[1]][j]+attr(pathxways[[1]],"match.length")[j]-1)
 			tab2 =rbind(tab2,tab[i,])
 			tab2$Path[length(tab2$Path)]=paste(bla,"-",unlist(mget(bla, KEGGPATHID2NAME)))
 			}
 			
 	}#endelse
 	}#endif

cat(noquote("Pathway Parsing: Writing results to disk\n"))
write.table(tab2, file=paste('R-Output/pathways_top',topnumber,'-',normalisierung,myfilter,'.txt', sep=""),sep=";")

	dir_log=c(dir_log,paste('R-Output/pathways_top',topnumber,'-',normalisierung,myfilter,'.txt', sep=""))
	file_log=c(file_log,"Pathway Mapping")

cat(noquote("Pathway Parsing: Abundance Diagram to PDF\n"))
tab3=table(tab2$Path[tab2$Path!="NA"])[order(table(tab2$Path[tab2$Path!="NA"]),decreasing = F)]
pdf(file=paste('R-Output/pathways_top',topnumber,'-',normalisierung,myfilter,'.pdf', sep=""), height=2+length(tab3[tab3>2])/10,width=10)
 	par(mar=c(3,20,0.1,0.2))
 	barplot(tab3[tab3>2], horiz=T, las=2, cex.names=0.75)
 	for (i in 1:(max(tab3)/5)){
 		abline(v=i*5, lty=3)
 	}
dev.off()

tab3=table(tab2$Path[tab2$Path!="NA"])[order(table(tab2$Path[tab2$Path!="NA"]),decreasing = F)]
pdf(file=paste('R-Output/pathways_top',topnumber,'-',normalisierung,myfilter,'_10.pdf', sep=""), height=2+length(tab3[tab3>9])/10,width=10)
 	par(mar=c(3,20,0.1,0.2))
 	barplot(tab3[tab3>9], horiz=T, las=2, cex.names=0.75)
 	for (i in 1:(max(tab3)/5)){
 		abline(v=i*5, lty=3)
 	}
dev.off()

dir_log=c(dir_log,paste('R-Output/pathways_top',topnumber,'-',normalisierung,myfilter,'.pdf', sep=""))
file_log=c(file_log,"Pathway Abundance Diagram")
}

########################################################################################
######## HEATMAP #######################################################################
########################################################################################

cat(noquote("Loading Advanced Functions\n"))

adv_log =c("Probesets2Heatmap(eset.norm, VectorOfProbesets)")

adv_log =c(adv_log, "PathSigHeat(eset.sig, PathwayKeggID)")



if(plotPathHeatmaps){
analysis_log=c(analysis_log,"Plotting Pathway Heatmaps")


dir.create("R-Output/Paths")
kegg<- mget(featureNames(eset.sig.norm), hgu133plus2PATH, ifnotfound = NA)
all.sig.paths = names(table(as.character(unlist(kegg)))[table(as.character(unlist(kegg)))>1])
for (i in 1:length(all.sig.paths)){
	PathSigHeat(eset.sig.norm,all.sig.paths[i])
}
}
dir_log=c(dir_log,"R-Output/Paths/")
file_log=c(file_log,"Pathway Heatmaps")


if(plotPathInteractions){
kegg<- mget(featureNames(eset.sig.norm), hgu133plus2PATH, ifnotfound = NA)
all.sig.paths = names(table(as.character(unlist(kegg)))[table(as.character(unlist(kegg)))>2])
for (i in 1:length(all.sig.paths)){
	cat(noquote(paste("Interaction Plot: KEGG", all.sig.paths[i],"\n", sep="")))
	InteractKeggPlot(tab,all.sig.paths[i])
}
}
dir_log=c(dir_log,"R-Output/Paths/")
file_log=c(file_log,"Pathway Interaction Maps")

#pathName= "04010"
#InteractKeggPlot(top, pathName)

########################################################################################
######## VENN DIAGRAMM #################################################################
########################################################################################

if(plotVenn){
cat(noquote("Venn Diagramm\n"))
analysis_log=c(analysis_log,"Plotting Venn Diagramm")

results <- classifyTestsF(fit_orig, p.value=min(p.value[irs]))
a <- vennCounts(results)

pdf(file=paste('R-Output/venn_top',topnumber,'-',normalisierung,myfilter,'.pdf', sep=""), height=10+topnumber/400,width=10)
vennDiagram(results,include=c("up","down"),counts.col=c("red","green"))
dev.off()
}
dir_log=c(dir_log,paste('R-Output/venn_top',topnumber,'-',normalisierung,myfilter,'.pdf', sep=""))
file_log=c(file_log,"Venn Diagramm")

########################################################################################
######## GO ENRICHMENT #################################################################
########################################################################################

if(analyzeGOstats){
analysis_log=c(analysis_log,"GoStat Analysis")

chi_gostats <- function(eset.sig.norm,eset.norm){

cat(noquote("GOstats: Obtaining Gene Clusters and Universe\n"))
geneCluster=as.vector(na.omit(getEG(featureNames(eset.sig.norm),"hgu133plus2")))
geneCluster2=as.vector(na.omit(getSYMBOL(featureNames(eset.sig.norm),"hgu133plus2")))
geneCluster3=featureNames(eset.sig.norm)

geneUniverse=as.vector(na.omit(getEG(featureNames(eset.norm),"hgu133plus2"))) 
geneUniverse2=as.vector(na.omit(getSYMBOL(featureNames(eset.norm),"hgu133plus2"))) 
geneUniverse3=featureNames(eset.norm)

gos = c("KEGG","MF","BP","CC")
gosigs=list()

cat(noquote("GOstats: Creating Outfile: "))
cat(noquote(paste('R-Output/GOgenes_','top',topnumber,'-PATH-',auswertung,'-',normalisierung,'-',myfilter,'.txt\n', sep="")))


outfile = paste('R-Output/goGenes_','top',topnumber,'-PATH-',normalisierung,'-',myfilter,'.txt', sep="")

dir_log=c(dir_log, outfile)
file_log=c(file_log,"GoStats Results")

sink(outfile)
writeLines(paste("","GOid","Type","Name","Pvalue","OddsRatio","ExpCount","Count","Size","GOGenesNumber","SigGenesNumber","SigGenesLinCorr","GOGenesSymbol","GOGenesEG","SigGenesSymbol","SigGenesEG",sep=";"), con=paste('R-Output/goGenes_','top',topnumber,'-PATH-',normalisierung,'-',myfilter,'.txt', sep=""))	
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
	cat(noquote(paste("GOstats: GO-Set excluded: ",gos[id] ," - ", i, " (",length(listofGOs),") ", sep="")))
}# end if (analyzeGO)
}# end for length(listofGOs))
}# end if (length(listofGOs)>0)

cat(noquote(" Making Barplot\n"))
cat(noquote(paste(c("",t(as.data.frame(genefunction))),collapse=";")),file=outfile, append=T, sep="\n")

sig.hgOverGO=summary(hgOverGO)[summary(hgOverGO)$Pvalue<0.1,]

sig.hgOverGO =sig.hgOverGO[order(sig.hgOverGO$Pvalue, decreasing=T),]

pdf(file=paste('R-Output/goGenes_',gos[id],'_top',topnumber,'-PATH-',normalisierung,'-',myfilter,'.pdf', sep=""), height=4+length(sig.hgOverGO$Count)/10,width=10)
 	par(mar=c(3,20,0.1,0.2))
	barplot(rbind(sig.hgOverGO$ExpCount,sig.hgOverGO$Count), col=c("gray80","gray20"), beside=T, horiz=T, names.arg=paste(sig.hgOverGO$KEGGID,sig.hgOverGO$Term,"- p =", round(sig.hgOverGO$Pvalue, digits=4)), las=2, cex.names=0.5)


 	for (i in 1:(max(sig.hgOverGO$Count)/5)){
 		abline(v=i*5, lty=3)
 	}
 	par(xpd=NA)
 	legend(-15,0,legend=c("Count", "Expected"), fill=c("gray20","gray80"))
dev.off()


}# end for MF/BP/CC/KEGG
}# end if analyzeGOstats
}
###### FROM HERE ON OLD/BAD STUFF, NOT USE, KEPT FOR REFERENCE ONLY
	
# 	##########################################################################################
	# ######## HEATMAPS ######################################################################
	# ########################################################################################
	
	# KEGG2heatmap("04810", sig.eset, "hgu133plus2")
	
	
	# if(length(heatmaps4GO)){
	
	# geneSub2plot <- c("204493_at",
#"213373_s_at",
#"202531_at",
#"203275_at",
#"204211_x_at",
#"207196_s_at",
#"1552648_a_at",
#"210405_x_at",
#"214329_x_at",
#"209546_s_at",
#"233110_s_at",
#"218696_at",
#"225164_s_at",
#"209808_x_at",
#"241985_at",
#"1555897_at",
#"204285_s_at",
#"218849_s_at",
#"203625_x_at",
#"206513_at",
#"219716_at",
#"223349_s_at",
#"211367_s_at",
#"226032_at",
#"216598_s_at",
#"218943_s_at",
#"224336_s_at",
#"219209_at",
#"208173_at",
#"208436_s_at",
#)
	# # einfach mit Probeset
	
	# eset2plot <- eset.norm[geneSub2plot, ] # Expressionsset für diese Probes generieren
	
	# heatmap.2(exprs(eset2plot), labRow=paste(getSYMBOL(featureNames(eset2plot), "hgu133plus2")," (",featureNames(eset2plot),")", sep=""), Rowv = F, Colv = F,dendrogram="none", lmat=rbind( c(0, 4,3), c(2,1,0 ), c(0,0,0) ), lwid=c(0.1, 3, 1 ), lhei= c(0.85,4,0.25),density.info="none", scale="row", key = T, symkey = T, trace ="none")
	 # # hier wird die eigentliche Heatmap erstellt!
	
	# ########################################################################################
	# ######## PATHWAYS ######################################################################
	# ########################################################################################
	

	# x <- layoutGraph(g.old)
	# fill <- c(rep(logcolO[1], length(logcol[logfcs>0])),rep(logcolO[2], length(logcol[logfcs<0])))
	# names(fill) <- c(names(logcol[logfcs>0]),names(logcol[logfcs<0]))
	# nodeRenderInfo(x) <- list(fill = fill)
	# edgeRenderInfo(x) <- list(fill = fill)
	
	# layout(mat=matrix(c(1,1,2,3), ncol=2, byrow=T))
	# plotKEGGgraph(g.old, y="twopi",defAttrs = list(rankdir="LR"))
	
	# KEGGgraphLegend()
	
	
	
	# image(as.matrix(seq(1,ar)), col=cols, yaxt="n", xaxt="n")
	# mtext("down-regulation", side=1,  at=0.03, line=2, cex=2)
	# mtext("up-regulation", side=1,  at=0.97, line=2, cex=2)
	# #dev.off()

