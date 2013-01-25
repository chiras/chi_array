#vennlist
setlistPLUS <- list()#all=featureNames(eset.norm))
setlistMINUS <- list()#all=featureNames(eset.norm))

for (v in 2:length(analysesList)){
		
	analysesList[[v]]$tab 		= chi_limma(	eset.norm[,analysesList[[v]]$use_subset],
												analysesList[[v]]$design,
												analysesList[[v]]$contrast.matrix,
												analysesList[[v]]$name)
												
	analysesList[[v]]$eset.sig 	= chi_filter(	eset.norm[,analysesList[[v]]$use_subset],
												analysesList[[v]]$tab,
												analysesList[[v]]$toDo$IQR_filter,
												analysesList[[v]]$toDo$logFC_filter,
												analysesList[[v]]$toDo$AFFX_filter,
												analysesList[[v]]$Filtering$IQR_intens_above,
												analysesList[[v]]$Filtering$IQR_intens_probes,
												analysesList[[v]]$Filtering$IQR_greater_than,
												analysesList[[v]]$Filtering$logFC_threshold)


	textlevel = leveltext("Writing filtered data to file","up",textlevel)
												
	analysesList[[v]]$tab.hit 	= chi_eset2sigs(analysesList[[v]]$eset.sig,analysesList[[v]]$tab)
	analysesList[[v]]$tab.sig 	= chi_eset2bins(analysesList[[v]]$eset.sig,analysesList[[v]]$tab)
	tx = length(analysesList[[v]]$tab.hit[,1])
	
	write.table(analysesList[[v]]$tab.hit, file=paste(PATHoutp,"/",analysis,'/results_',analysesList[[v]]$name,'_top',tx,'-',normalisierung,'.txt', sep=""),sep=";")

	if(analysesList[[v]]$toDo$includeMFA){
		setlistPLUS[[analysesList[[v]]$name]] = analysesList[[v]]$tab$ID[analysesList[[v]]$tab.sig==1]
		setlistMINUS[[analysesList[[v]]$name]] = analysesList[[v]]$tab$ID[analysesList[[v]]$tab.sig==-1]
	}

	textlevel = leveltext("Analyses completed\n","keep",textlevel)
	textlevel = leveltext("","down",textlevel)
		
}

#### Multifaktorial Steps
# Venn
textlevel = leveltext("Multifactorial overlap","keep",textlevel)
textlevel = leveltext("Plotting Venn diagrams","up",textlevel)
chi_venn(setlistPLUS,setlistMINUS)

textlevel = leveltext("Writing overlaps to file\n","keep",textlevel)
overlap = chi_overlap(setlistPLUS,setlistMINUS)
write.list(overlap,filename=paste(PATHoutp,"/",analysis,"/overlap_probesets.txt",sep=""))


# Heatmaps
textlevel = leveltext("Plotting heatmaps","down",textlevel)
textlevel = leveltext("Converting Genes to Probesets","up",textlevel)

if (length(mygroups)>0){
for(z in 1:length(mygroups)){
	y = unique(mygroups[[z]])
	myprobes[[names(mygroups[z])]]=c()
	for (r in 1:length(y)){
		myprobes[[paste("gl",names(mygroups[z]),sep="_")]] <- c(myprobes[[names(mygroups[z])]],names(geneinfo[geneinfo==y[r]]))
	}
}	
}

textlevel = leveltext("Converting Pathways to Probesets","keep",textlevel)

if (length(mypaths)>0){
for(z in 1:length(mypaths)){
	for (r in 1:length(mypaths[[z]])){
		 y <- pathinfo[[mypaths[[z]][r]]]
	}
	myprobes[[paste("pw",names(mypaths[z]),sep="_")]] = unique(y)

}
}
	
############# MIT PROBESETS DANN GEMEINSAM WEITER
textlevel = leveltext("Preparing Heatmaps","keep",textlevel)

for(z in 1:length(myprobes)){
	y = unique(myprobes[[z]])
	nameplot=names(myprobes)[[z]]
	tmp_lst=list()

	for (i in 2:length(analysesList)){
		tmp_lst[[analysesList[[i]]$name]]=analysesList[[i]]$tab[which(analysesList[[i]]$tab$ID==""),c(1,2,3,6,10)]

		for (counterA in 1:length(y)){
			tmp_lst[[analysesList[[i]]$name]] = rbind(tmp_lst[[analysesList[[i]]$name]],analysesList[[i]]$tab[which(analysesList[[i]]$tab$ID==y[counterA]),c(1,2,3,6,10)])
		}
	tmp_lst[[analysesList[[i]]$name]] = tmp_lst[[analysesList[[i]]$name]][with(tmp_lst[[analysesList[[i]]$name]], order(Symbol,ID)), ]
					# sorting the lists to make them directly comparable
	}
	
	logfcvalues = data.frame(scale=tmp_lst[[analysesList[[i]]$name]]$ID)
	for (i in 1:length(tmp_lst)){
			logfcvalues = cbind(logfcvalues,tmp_lst[[i]]$logFC)
			names(logfcvalues)[length(logfcvalues)]<-names(tmp_lst)[i]
			}	
			
	logfcvalues$scale = seq(from=max(logfcvalues[2:length(logfcvalues)]), to= -max(logfcvalues[2:length(logfcvalues)]), by= -(max(logfcvalues[2:length(logfcvalues)])*2)/(length(tmp_lst[[1]]$ID)-1))

		heat.matrix = as.matrix(logfcvalues)

	row.names(heat.matrix) <- paste(tmp_lst[[1]]$Symbol," (",tmp_lst[[1]]$ID,")", sep="")
	hmcol<-brewer.pal(11,"RdBu")

  textlevel = leveltext(paste("Plotting Heatmaps (",nameplot,")", sep=""),"keep",textlevel)
    
pdf(file=paste(PATHoutp,"/",analysis,'/MFA-Plots/Heat_',nameplot,'.pdf', sep=""), height=2+length(tmp_lst[[1]]$ID)/5,width=(length(tmp_lst)+1)*2)
heatmap.chi(heat.matrix, key=F, dendrogram="none", trace="none",col= rev(hmcol),   margins = c(12, 12), cellnote=round(heat.matrix, digits=3), notecol="black", Rowv=NA, Colv=NA, scale="none",
 add.expr = abline(h = c(0.5,cumsum(rev(table(tmp_lst[[1]]$Symbol)))+0.5), v=c(0.5,1.5,5.5), lwd = 3), xlab= nameplot, cexRow = 1.05,
 cexCol=1,labCol = c(paste("Scale (",analysesList$general$scalename,")",sep=""),colnames(heat.matrix)[2:length(colnames(heat.matrix))]),
 lmat=rbind( c(0,0,0),  c(0,1,0), c(0,0,0) ), lwid=c(0.005,0.99,0.005), lhei=c( 0.005,0.99, 0.005 ), keysize=0)


#heatmap.2(heat.matrix, key=F, dendrogram="none", trace="none",lwid=c(5,5), col= hmcol, lhei=c(5,5),   margins = c(12, 12), cellnote=round(heat.matrix, digits=3), notecol="black", Rowv=NA, Colv=NA, scale="none",add.expr = abline(h = c(0.5,cumsum(rev(table(tab_wisp1_x_sort_msc$Symbol)))+0.5), v=c(0.5,1.5,5.5), lwd = 3), xlab= nameplot, labCol = c("Scale (sh-scr)","WISP1 MSCs","WISP3 MSCs","WISP1 Cho","WISP1 Cho"))
dev.off()
    
}
textlevel = leveltext("Done with Heatmaps\n","keep",textlevel)


#missing: SA !!!!

########################################################################################
######## LIMMA AFTER FILTER ############################################################
########################################################################################

########################################################################################
######## PATHWAY ANALYSIS ##############################################################
########################################################################################

########################################################################################
######## HEATMAP #######################################################################
########################################################################################

########################################################################################
######## VENN DIAGRAMM #################################################################
########################################################################################

########################################################################################
######## GO ENRICHMENT #################################################################
########################################################################################
textlevel = leveltext("GO/KEGG enrichment analysis","keep",0)
for (i in 2:length(analysesList)){
	if (analysesList[[i]]$toDo$analyzeGOstats){
		textlevel = leveltext(analysesList[[i]]$name,"up",textlevel)
		chi_gostats (analysesList[[i]]$name, analysesList[[i]]$eset.sig,eset.norm)
		textlevel = leveltext("","down",textlevel)
	}
}

# 	##########################################################################################
	# ######## HEATMAPS ######################################################################
	# ########################################################################################
	
	# ########################################################################################
	# ######## PATHWAYS ######################################################################
	# ########################################################################################
	
textlevel = leveltext("Pathway analysis","keep",0)
for (i in 2:length(analysesList)){
	if (analysesList[[i]]$toDo$parseGenes2Pathways){
		textlevel = leveltext(analysesList[[i]]$name,"up",textlevel)
		chi_genes2Pathways (analysesList[[i]]$name, analysesList[[i]]$tab.hit)
		textlevel = leveltext("","down",textlevel)
	}
}