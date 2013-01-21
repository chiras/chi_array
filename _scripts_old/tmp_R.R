if(plotVenn){
cat(noquote("Venn Diagramm\n"))
analysis_log=c(analysis_log,"Plotting Venn Diagramm")
irs <- grep("AFFX",featureNames(eset.norm)) ##### rausfinden, welchen p-value unregulierte gene haben "AFFX"


source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")
chi_eset2sigs <- function (eset, tab.orig){
	inters= !tab$ID %in% featureNames(eset)
	tab.tmp = 2*(tab$t>0)-1
	tab.tmp[inters] = 0
	return(tab.tmp)
}


setlistPLUS <- list()#all=featureNames(eset.norm))
setlistMINUS <- list()#all=featureNames(eset.norm))

for (i in 2:length(analysesList)){
	results = chi_eset2sigs(xx, tab)

	setlistPLUS[[analysesList[[i]]$name]] = w1plus=tab$ID[results==1]
	setlistMINUS[[analysesList[[i]]$name]] = w1minus=tab$ID[results==-1]

}

OLlistPLUS <- overLapper(setlist=setlistPLUS, sep="_", type="vennsets")
OLlistMINUS <- overLapper(setlist=setlistMINUS, sep="_", type="vennsets")

counts <- list(sapply(OLlistPLUS$Venn_List, length), sapply(OLlistMINUS$Venn_List, length));

vennPlot(counts=counts, mysub="Top: Up-regulated; Bottom: Down-regulated", yoffset=c(0.3, -0.2))



#results <- classifyTestsF(xxx)

#p_values = tabx$wisp1_msc$adj.P.Val
#affx_probes=grep("AFFX",tabx$wisp1_msc$ID)

#all= w1plus=tabx$wisp1_msc$ID
#p.all = 0.01

#p.w1 = p.all1 #min(tabx$wisp1_msc$adj.P.Val[affx_probes])

w1plus=tab$ID[results==1]
w1minus=tab$ID[results==-1]

#table(classifyTestsF(tabx$wisp1_msc$t, p.value=p.w1))

p.w2 = p.all1 #min(tabx$wisp1_msc$adj.P.Val[irs])
w2plus=tabx$wisp3_msc$ID[classifyTestsF(tabx$wisp3_msc$t, p.value=p.w2)==1]
w2minus=tabx$wisp3_msc$ID[classifyTestsF(tabx$wisp3_msc$t, p.value=p.w2)==-1]





#pdf(file=paste('R-Output/venn_top',topnumber,'-',normalisierung,myfilter,'.pdf', sep=""), height=10+topnumber/400,width=10)
#vennDiagram(results,include=c("up","down"),counts.col=c("red","green"))
#dev.off()
dir_log=c(dir_log,paste('R-Output/venn_top',topnumber,'-',normalisierung,myfilter,'.pdf', sep=""))
file_log=c(file_log,"Venn Diagramm")

}



source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")

setlist <- list(A=sample(letters, 18), B=sample(letters, 16), C=sample(letters, 20), D=sample(letters, 22), E=sample(letters, 18), F=sample(letters, 22, replace=T))

OLlist <- overLapper(setlist=setlist, sep="", type="vennsets"); OLlist; names(OLlist)

setlist4 <- setlist[1:4]
OLlist4 <- overLapper(setlist=setlist4, sep="_", type="vennsets")

counts <- list(sapply(OLlist4$Venn_List, length), sapply(OLlist4$Venn_List, length));

vennPlot(counts=counts, mysub="Top: var1; Bottom: var2", yoffset=c(0.3, -0.2))