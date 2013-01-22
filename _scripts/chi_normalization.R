##################################################################################################
#todo: 
#- Venn up & Down + AFFX cut + Filter
#- heatmaps
#- interaktion


########################################################################################
######## READ & SETTING ENVIRONMENT ####################################################
########################################################################################
auswertung 		= analysis
normalisierung 	= normalization

textlevel = leveltext("Loading Chip Annotation Data","up",textlevel)

load(paste(PATHress,"/GenesetsU133plus2.RData",sep=""))												# Chip Annotation data


if (use_norm){ 																			# when new normalization is used
	
########################################################################################
######## READ RAW DATA #################################################################
########################################################################################
	textlevel = leveltext("Normalization","keep",textlevel)

	textlevel = leveltext("Loading Raw Data","up",textlevel)
	pd <- read.AnnotatedDataFrame(phenodataFile, header = TRUE, row.names = 1, sep="") 	# read Phenotype data
	pData(pd) 																			# Map Phenotype Data
	rawAffyData <- ReadAffy(filenames= rownames(pData(pd)), phenoData=pd) 				# read raw data
	
	# look at genes, probes, probesets
	length(featureNames(rawAffyData))	# gene (probesets)
	length(probeNames(rawAffyData))		# probes -> 11 probes pro probeset

########################################################################################
######## QUALITY CONTROL ###############################################################
########################################################################################

	if (qualityPlots){	

	textlevel = leveltext("Quality Control","keep",textlevel)

	# plot boxplots
	textlevel = leveltext("Boxplots","up",textlevel)
	pdf(file=paste(PATHoutp,"/",analysis,"/Q-Plots/box_nonorm.pdf", sep=""))
		boxplot(rawAffyData)
	dev.off()

	# plot densities
	textlevel = leveltext("Density Plots","keep",textlevel)
	pdf(file=paste(PATHoutp,"/",analysis,"/Q-Plots/dens_nonorm.pdf", sep=""))	
	plot(density(log2(intensity(rawAffyData[, 1]))), col =0, main="Densities", ylim=c(0,1))
	for(i in 1:dim(intensity(rawAffyData))[2])
	{
		lines(density(log2(intensity(rawAffyData[,i]))), col = i)
	}
	dev.off()

	textlevel = leveltext("MA Plots","keep",textlevel)
	cat(noquote("Normalization: Quality Control: MA-Plots\n"))
	pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/MA_nonorm%03d.pdf', sep=""), onefile=FALSE) 	
	MAplot(rawAffyData)
	dev.off()

	# look at RNA degradation
	textlevel = leveltext("RNA Degradation Plots","keep",textlevel)
	pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/rna-deg_nonorm.pdf', sep=""))	
	deg <- AffyRNAdeg(rawAffyData, log.it=FALSE)
	plotAffyRNAdeg(deg)
	dev.off()

	}

########################################################################################
######## Normalization #################################################################
########################################################################################

	# use expresso and summarization for normalization
	textlevel = leveltext("Calculating","down",textlevel)

	if (normalisierung == "quantiles"){
		textlevel = leveltext("Quantiles","up",textlevel)
		eset.norm <- expresso(rawAffyData, bg.correct=FALSE, normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		write.exprs(eset.norm,file=paste(PATHdata,"/",analysis,'/eset_norm-',normalisierung,'.txt', sep=""))
	}else if (normalisierung == "rma"){
		textlevel = leveltext("RMA","up",textlevel)
		eset.norm <- expresso(rawAffyData,bgcorrect.method="rma", normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		write.exprs(eset.norm,file=paste(PATHdata,"/",analysis,'/eset_norm-',normalisierung,'.txt', sep=""))
	}else if (normalisierung == "vsn"){
		textlevel = leveltext("VSN (currently disabled)","up",textlevel)
#		eset.norm <- expresso(rawAffyData, bg.correct=FALSE, normalize.method="vsn", pmcorrect.method="pmonly", summary.method="medianpolish")
#		write.exprs(eset.norm,file=paste(PATHdata,"/",analysis,'/eset_norm-',normalisierung,'.txt', sep=""))
	}else{
		textlevel = leveltext(" Quantiles/RMA/VSN (VSN currently disabled)","up",textlevel)
		eset.quantiles <- expresso(rawAffyData, bg.correct=FALSE, normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		eset.rma <- expresso(rawAffyData,bgcorrect.method="rma", normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
#		eset.vsn <- expresso(rawAffyData, bg.correct=FALSE, normalize.method="vsn", pmcorrect.method="pmonly", summary.method="medianpolish")

		# look at normalise data
		textlevel = leveltext("Normalized boxplots","keep",textlevel)

		par(mfrow=c(1,3))
		pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/box_norm%03d.pdf', sep=""), onefile=FALSE) 	
			boxplot(data.frame(exprs(eset.rma)))
			boxplot(data.frame(exprs(eset.quantiles)))
	#		boxplot(data.frame(exprs(eset.vsn)))
		dev.off()

		# look at density plots
		textlevel = leveltext("Normalized Density Plots","keep",textlevel)
		pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/dens_norm%03d.pdf', sep=""), onefile=FALSE) 	
		plot(density(log2(exprs(eset.rma[, 1]))), col=0, main="RMA")
		for(i in 1:dim(exprs(eset.rma))[2])
		{
			lines(density(log2(exprs(eset.rma[,i]))), col=i)

		}

		plot(density(log2(exprs(eset.quantiles[, 1]))), col=0, main="Quantiles")
		for(i in 1:dim(exprs(eset.quantiles))[2])
		{
			lines(density(log2(exprs(eset.quantiles[,i]))), col=i)
		}

	#	plot(density(log2(exprs(eset.vsn[, 1]))), col =0, main="VSN")
	#	for(i in 1:dim(exprs(eset.vsn))[2])
	#	{
	#		lines(density(log2(exprs(eset.vsn[,i]))), col = i)
	#	}
		dev.off()

		#write.exprs(eset.vsn,file=paste('R-Output/eset_norm-vsn.txt', sep=""))

		textlevel = leveltext("Writing normalized data","down",textlevel)
		write.exprs(eset.quantiles,file=paste(PATHdata,"/",analysis,'/eset_norm-quantiles.txt', sep=""))
		write.exprs(eset.rma,file=paste(PATHdata,"/",analysis,'/eset_norm-rma.txt', sep=""))

		textlevel = leveltext("No norm method specified, using RMA. Alternatives are saved to disk","keep",textlevel)
		eset.norm=eset.rma
	}

	if (qualityPlots){	
	textlevel = leveltext("Quality Control after Normalization","keep",textlevel)

	# plot boxplots
	textlevel = leveltext("Boxplots","up",textlevel)

	pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/box_norm.pdf', sep=""))
		boxplot(data.frame(exprs(eset.norm)))
	dev.off()

	# plot densities
	textlevel = leveltext("Density Plots","keep",textlevel)
	pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/dens_norm.pdf', sep=""))	
		plot(density(log2(exprs(eset.norm[, 1]))), col=0, main="Normalization")
		for(i in 1:dim(exprs(eset.norm))[2])
		{
			lines(density(log2(exprs(eset.norm[,i]))), col=i)

		}
	dev.off()
	}
	textlevel = leveltext("","down",textlevel)
	textlevel = leveltext("","down",textlevel)

}else{
	textlevel = leveltext("Reading pre-normalized data from file","keep",textlevel)

	eset.norm=readExpressionSet(paste('eset_norm-',auswertung,'-',normalisierung,'.txt', sep="")) # reading already normalized data
}

	textlevel = leveltext("(Checking Subsets) - Currently disabled!","keep",textlevel)

#if (length(use_subset)==length(sampleNames(eset.norm))){
	
#}else{
#	cat(noquote("Limiting Data to subset, if specified\n"))
#	eset.norm = eset.norm[,use_subset] 
#	analysis_log=c(analysis_log,"Restricting to subset")
#}

textlevel = leveltext("Renaming Samples","keep",textlevel)
	textlevel = leveltext("","up",textlevel)

for (zz in 1:length(sampleNames(eset.norm))){
	textlevel = leveltext(paste (sampleNames(eset.norm)[zz], "-->", analysesList$general$sampleName[zz]),"keep",textlevel)
	}
# renaming datasets	
sampleNames(eset.norm) <- analysesList$general$sampleNames 	
textlevel = leveltext("Done reading data!\n","down",textlevel)



