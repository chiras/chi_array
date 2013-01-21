
#vennlist
setlistPLUS <- list()#all=featureNames(eset.norm))
setlistMINUS <- list()#all=featureNames(eset.norm))

for (i in 2:length(analysesList)){
		
	analysesList[[i]]$tab 		= chi_limma(	eset.norm[,analysesList[[i]]$use_subset],
												analysesList[[i]]$design,
												analysesList[[i]]$contrast.matrix,
												analysesList[[i]]$name)
												
	analysesList[[i]]$eset.sig 	= chi_filter(	eset.norm[,analysesList[[i]]$use_subset],
												analysesList[[i]]$tab,
												analysesList[[i]]$toDo$IQR_filter,
												analysesList[[i]]$toDo$logFC_filter,
												analysesList[[i]]$toDo$AFFX_filter,
												analysesList[[i]]$Filtering$IQR_intens_above,
												analysesList[[i]]$Filtering$IQR_intens_probes,
												analysesList[[i]]$Filtering$IQR_greater_than,
												analysesList[[i]]$Filtering$logFC_threshold)
												
	analysesList[[i]]$tab.hit 	= chi_eset2sigs(analysesList[[i]]$eset.sig,analysesList[[i]]$tab)
	analysesList[[i]]$tab.sig 	= chi_eset2bins(analysesList[[i]]$eset.sig,analysesList[[i]]$tab)
	
	if(analysesList[[i]]$toDo$includeMFA){
		setlistPLUS[[analysesList[[i]]$name]] = analysesList[[i]]$tab$ID[analysesList[[i]]$tab.sig==1]
		setlistMINUS[[analysesList[[i]]$name]] = analysesList[[i]]$tab$ID[analysesList[[i]]$tab.sig==-1]
	}
	cat("\n")		
}

#### Multifaktorial Steps
# Venn
cat("Plotting Venn Diagrams\n")		

chi_venn(setlistPLUS,setlistMINUS)

# Heatmaps

cat(noquote("Converting Genes to Probesets\n"))
analysis_log=c(analysis_log,"Converting Genes to Probesets")

if (length(mygroups)>0){
for(z in 1:length(mygroups)){
	y = unique(mygroups[[z]])
	myprobes[[names(mygroups[z])]]=c()
	for (r in 1:length(y)){
		myprobes[[names(mygroups[z])]] <- c(myprobes[[names(mygroups[z])]],names(geneinfo[geneinfo==y[r]]))
	}
}	
}

cat(noquote("Converting Pathways to Probesets\n"))
analysis_log=c(analysis_log,"Converting Pathways to Probesets")

if (length(mypaths)>0){
for(z in 1:length(mypaths)){
	for (r in 1:length(mypaths[[z]])){
		 y <- pathinfo[[mypaths[[z]][r]]]
	}
	myprobes[[names(mypaths[z])]] = unique(y)

}
}
	
############# MIT PROBESETS DANN GEMEINSAM WEITER
cat(noquote("Preparing Heatmaps\n"))
analysis_log=c(analysis_log,"Preparing Heatmaps")


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

  cat(noquote(paste("Plotting Heatmaps (",nameplot,")\n", sep="")))
  analysis_log=c(analysis_log,paste("Plotting Heatmaps (",nameplot,")", sep=""))
    
pdf(file=paste('R-Output/MFA-Plots/Heat_',nameplot,'.pdf', sep=""), height=2+length(tmp_lst[[1]]$ID)/5,width=(length(tmp_lst)+1)*2)
heatmap.chi(heat.matrix, key=F, dendrogram="none", trace="none",col= rev(hmcol),   margins = c(12, 12), cellnote=round(heat.matrix, digits=3), notecol="black", Rowv=NA, Colv=NA, scale="none",
 add.expr = abline(h = c(0.5,cumsum(rev(table(tmp_lst[[1]]$Symbol)))+0.5), v=c(0.5,1.5,5.5), lwd = 3), xlab= nameplot, cexRow = 1.05,
 #labCol = c(paste("Scale (",scalename,")",sep=""),complist$wisp1_msc$name,complist$wisp3_msc$name,complist$wisp1_cho$name ,complist$wisp3_cho$name),
 lmat=rbind( c(0,0,0),  c(0,1,0), c(0,0,0) ), lwid=c(0.005,0.99,0.005), lhei=c( 0.005,0.99, 0.005 ), keysize=0)
#heatmap.2(heat.matrix, key=F, dendrogram="none", trace="none",lwid=c(5,5), col= hmcol, lhei=c(5,5),   margins = c(12, 12), cellnote=round(heat.matrix, digits=3), notecol="black", Rowv=NA, Colv=NA, scale="none",add.expr = abline(h = c(0.5,cumsum(rev(table(tab_wisp1_x_sort_msc$Symbol)))+0.5), v=c(0.5,1.5,5.5), lwd = 3), xlab= nameplot, labCol = c("Scale (sh-scr)","WISP1 MSCs","WISP3 MSCs","WISP1 Cho","WISP1 Cho"))
dev.off()
    
    dir_log=c(dir_log,paste('R-Output/MFA-Plots/MFA_',nameplot,'.pdf', sep=""))
    file_log=c(file_log,paste("MFA Heatmap (",nameplot,")", sep=""))

}


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

for (i in 2:length(analysesList)){
	chi_gostats (analysesList[[i]]$name, analysesList[[i]]$eset.sig,eset.norm)
}

# 	##########################################################################################
	# ######## HEATMAPS ######################################################################
	# ########################################################################################
	
	# ########################################################################################
	# ######## PATHWAYS ######################################################################
	# ########################################################################################
	
