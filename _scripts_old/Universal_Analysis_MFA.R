

########################################################################################################################################
########################################################################################################################################

##### from here on automatically

fitx=list()
tabx=list()
vennx=list()

for (r in 1:length(complist)){
	complist[[r]]$contrast.matrix <- makeContrasts(complist[[r]]$comparison, levels=complist[[r]]$design)  	
	complist[[r]]$fit <- lmFit(eset.norm[,complist[[r]]$subset], complist[[r]]$design)	

fit2 <- contrasts.fit(complist[[r]]$fit, complist[[r]]$contrast.matrix)
fit2 <- eBayes(fit2)
fit_orig <- fit2

cat(noquote(paste("Limma Analysis: Mapping Annotation Data (",complist[[r]]$name,")\n", sep="")))
fit2$genes$Symbol=getSYMBOL(fit_orig$genes$ID, "hgu133plus2")
fit2$genes$GeneName <- unlist(mget(complist[[r]]$fit_orig$genes$ID, hgu133plus2GENENAME))
fit2$genes$Path <- paste(mget(complist[[r]]$fit_orig$genes$ID, hgu133plus2PATH))
fit2$genes$EG <- getEG(fit_orig$genes$ID, "hgu133plus2")

pathinfo <- as.list(hgu133plus2PATH2PROBE)
pregeneinfo <-  hgu133plus2SYMBOL
mapped_probes <- mappedkeys(pregeneinfo)
geneinfo <- as.list(pregeneinfo[mapped_probes])

cat(noquote("   --> Calculating Statistics\n"))
F.stat <- FStat(fit2)
p.value <- pf(F.stat,df1=attr(F.stat,"df1"),df2=attr(F.stat,"df2"),lower.tail=FALSE)
topnumber = length(featureNames(eset.norm))
tab <- topTable(fit2, adjust = "fdr", number=topnumber) # FDR corrected p-values nach dem cut off
head(tab)
    analysis_log=c(analysis_log,paste("Limma Analysis (",complist[[r]]$name,")", sep=""))

tabx[[names(complist)[[r]]]]<-tab
}
names(tabx)


#####

cat(noquote("Converting Genes to Probesets\n"))
analysis_log=c(analysis_log,"Converting Genes to Probesets")

if (length(mygroups)>0){
for(z in 1:length(mygroups)){
	y = unique(mygroups[[z]])
	myprobes[[names(mygroups[z])]]=c()
	for (r in 1:length(y)){
		myprobes[[names(mygroups[z])]] <- c(myprobes[[names(mygroups[z])]],names(geneinfo[geneinfo==y[r]]))
	}
}}	

cat(noquote("Converting Pathways to Probesets\n"))
analysis_log=c(analysis_log,"Converting Pathways to Probesets")

if (length(mypaths)>0){
for(z in 1:length(mypaths)){
	for (r in 1:length(mypaths[[z]])){
		 y <- pathinfo[[mypaths[[z]][r]]]
	}
	myprobes[[names(mypaths[z])]] = unique(y)

}}	

cat(noquote("Preparing Heatmaps\n"))
analysis_log=c(analysis_log,"Preparing Heatmaps")



for(z in 1:length(myprobes)){
	y = unique(myprobes[[z]])
	nameplot=names(myprobes)[[z]]
	


####
tab_wisp1 = tabx$wisp1_msc
tab_wisp3 = tabx$wisp3_msc
tab_wisp1_x=tabx$wisp1_msc[which(tabx$wisp1_msc$ID==""),c(1,2,3,6,10)]

for (i in 1:length(y)){
	tab_wisp1_x = rbind(tab_wisp1_x,tab_wisp1[which(tab_wisp1$ID==y[i]),c(1,2,3,6,10)])
}

tab_wisp3_x=tab_wisp3[which(tab_wisp3$ID==""),c(1,2,3,6,10)]

for (i in 1:length(y)){
	tab_wisp3_x = rbind(tab_wisp3_x,tab_wisp3[which(tab_wisp3$ID==y[i]),c(1,2,3,6,10)])
}

tab_wisp1_x_sort_msc = tab_wisp1_x[with(tab_wisp1_x, order(Symbol,ID)), ]
tab_wisp3_x_sort_msc = tab_wisp3_x[with(tab_wisp3_x, order(Symbol,ID)), ]

####
tab_wisp1 = tabx$wisp1_cho
tab_wisp3 = tabx$wisp3_cho
tab_wisp1_x=tab_wisp1[which(tab_wisp1$ID==""),c(1,2,3,6,10)]

for (i in 1:length(y)){
	tab_wisp1_x = rbind(tab_wisp1_x,tab_wisp1[which(tab_wisp1$ID==y[i]),c(1,2,3,6,10)])
}

tab_wisp3_x=tab_wisp3[which(tab_wisp3$ID==""),c(1,2,3,6,10)]

for (i in 1:length(y)){
	tab_wisp3_x = rbind(tab_wisp3_x,tab_wisp3[which(tab_wisp3$ID==y[i]),c(1,2,3,6,10)])
}

tab_wisp1_x_sort_cho = tab_wisp1_x[with(tab_wisp1_x, order(Symbol,ID)), ]
tab_wisp3_x_sort_cho = tab_wisp3_x[with(tab_wisp3_x, order(Symbol,ID)), ]
####


if(sum(tab_wisp1_x_sort_msc$ID == tab_wisp3_x_sort_msc$ID)!=length(tab_wisp3_x_sort_msc$ID)){print("******************WRONG***************")}
if(sum(tab_wisp1_x_sort_cho$ID == tab_wisp3_x_sort_cho$ID)!=length(tab_wisp1_x_sort_cho$ID)){print("******************WRONG***************")}
if(sum(tab_wisp1_x_sort_msc$ID == tab_wisp3_x_sort_cho$ID)!=length(tab_wisp3_x_sort_msc$ID)){print("******************WRONG***************")}




pseudomax = max(c(tab_wisp1_x_sort_msc$logFC,tab_wisp3_x_sort_msc$logFC,tab_wisp1_x_sort_cho$logFC,tab_wisp3_x_sort_cho$logFC,3))
pseudo = seq(from=pseudomax, to= -pseudomax, by= -(pseudomax*2)/(length(tab_wisp1_x_sort_msc$ID)-1))


heat.matrix = as.matrix(data.frame(scale=pseudo, WISP1_msc=tab_wisp1_x_sort_msc$logFC,WISP3_msc=tab_wisp3_x_sort_msc$logFC,WISP1_cho=tab_wisp1_x_sort_cho$logFC,WISP3_cho=tab_wisp3_x_sort_cho$logFC))

row.names(heat.matrix) <- paste(tab_wisp1_x_sort_msc$Symbol," (",tab_wisp1_x_sort_msc$ID,")", sep="")

hmcol<-brewer.pal(11,"RdBu")


  cat(noquote(paste("Plotting Heatmaps (",nameplot,")\n", sep="")))
    analysis_log=c(analysis_log,paste("Plotting Heatmaps (",nameplot,")", sep=""))
    
pdf(file=paste('R-Output/MFA-Plots/MFA_',nameplot,'.pdf', sep=""), height=2+length(tab_wisp1_x$ID)/5,width=10)
heatmap.chi(heat.matrix, key=F, dendrogram="none", trace="none",col= rev(hmcol),   margins = c(12, 12), cellnote=round(heat.matrix, digits=3), notecol="black", Rowv=NA, Colv=NA, scale="none",
 add.expr = abline(h = c(0.5,cumsum(rev(table(tab_wisp1_x_sort_msc$Symbol)))+0.5), v=c(0.5,1.5,5.5), lwd = 3), xlab= nameplot, cexRow = 1.05,labCol = c(paste("Scale (",scalename,")",sep=""),complist$wisp1_msc$name,complist$wisp3_msc$name,complist$wisp1_cho$name ,complist$wisp3_cho$name),lmat=rbind( c(0,0,0),  c(0,1,0), c(0,0,0) ), lwid=c(0.005,0.99,0.005), lhei=c( 0.005,0.99, 0.005 ), keysize=0)
#heatmap.2(heat.matrix, key=F, dendrogram="none", trace="none",lwid=c(5,5), col= hmcol, lhei=c(5,5),   margins = c(12, 12), cellnote=round(heat.matrix, digits=3), notecol="black", Rowv=NA, Colv=NA, scale="none",add.expr = abline(h = c(0.5,cumsum(rev(table(tab_wisp1_x_sort_msc$Symbol)))+0.5), v=c(0.5,1.5,5.5), lwd = 3), xlab= nameplot, labCol = c("Scale (sh-scr)","WISP1 MSCs","WISP3 MSCs","WISP1 Cho","WISP1 Cho"))
dev.off()
    
    dir_log=c(dir_log,paste('R-Output/MFA-Plots/MFA_',nameplot,'.pdf', sep=""))
    file_log=c(file_log,paste("MFA Heatmap (",nameplot,")", sep=""))

}

