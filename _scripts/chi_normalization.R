########################################################################################
######## READ & SETTING ENVIRONMENT ####################################################
########################################################################################
auswertung 		= analysis
normalisierung 	= normalization

textlevel = leveltext("","up",textlevel)

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

	if (analysesList$general$qualityPlots){	

	textlevel = leveltext("Quality Control","keep",textlevel)

	# plot boxplots
	textlevel = leveltext("Boxplots","up",textlevel)
	pdf(file=paste(PATHoutp,"/",analysis,"/Q-Plots/BOX_nonorm.pdf", sep=""))
		boxplot(rawAffyData)
	dev.off()

	# plot densities
	textlevel = leveltext("Density Plots","keep",textlevel)
	pdf(file=paste(PATHoutp,"/",analysis,"/Q-Plots/DENS_nonorm.pdf", sep=""))	
	plot(density(log2(intensity(rawAffyData[, 1]))), col =0, main="Densities", ylim=c(0,1))
	for(i in 1:dim(intensity(rawAffyData))[2])
	{
		lines(density(log2(intensity(rawAffyData[,i]))), col = i)
	}
	dev.off()

	textlevel = leveltext("MA Plots","keep",textlevel)
	pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/MA_nonorm%03d.pdf', sep=""), onefile=FALSE) 	
	MAplot(rawAffyData)
	dev.off()

	# look at RNA degradation
	textlevel = leveltext("RNA Degradation Plots","keep",textlevel)
	pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/RNAdeg_nonorm.pdf', sep=""))	
	deg <- AffyRNAdeg(rawAffyData, log.it=FALSE)
	plotAffyRNAdeg(deg)
	dev.off()
	textlevel = leveltext("","down",textlevel)

	}

########################################################################################
######## Normalization #################################################################
########################################################################################

	# use expresso and summarization for normalization
	textlevel = leveltext("Calculating","keep",textlevel)

	if (normalisierung == "quantiles"){
		textlevel = leveltext("Quantiles","up",textlevel)
		eset.quantiles <- chi_expresso(rawAffyData, bg.correct=FALSE, normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		write.exprs(eset.quantiles,file=paste(PATHdata,"/",analysis,'/eset_norm-',normalisierung,'.txt', sep=""))
	}else if (normalisierung == "rma"){
		textlevel = leveltext("RMA background-corrected quantiles","up",textlevel)
		eset.rma <- chi_expresso(rawAffyData,bgcorrect.method="rma", normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		write.exprs(eset.rma,file=paste(PATHdata,"/",analysis,'/eset_norm-',normalisierung,'.txt', sep=""))
	}else if (normalisierung == "vsn"){
		textlevel = leveltext("VSN (currently disabled)","up",textlevel)
		stop("Please choose a different normalization method")
#		eset.vsn <- chi_expresso(rawAffyData, bg.correct=FALSE, normalize.method="vsn", pmcorrect.method="pmonly", summary.method="medianpolish")
#		write.exprs(eset.vsn,file=paste(PATHdata,"/",analysis,'/eset_norm-',normalisierung,'.txt', sep=""))
	}else{
		textlevel = leveltext("Quantiles","up",textlevel)
		eset.quantiles <- chi_expresso(rawAffyData, bg.correct=FALSE, normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		textlevel = leveltext("","down",textlevel)

		textlevel = leveltext("RMA background-corrected quantiles","up",textlevel)
		eset.rma <- chi_expresso(rawAffyData,bgcorrect.method="rma", normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		textlevel = leveltext("","down",textlevel)
		
		textlevel = leveltext("VSN currently disabled!","up",textlevel)
#		eset.vsn <- chi_expresso(rawAffyData, bg.correct=FALSE, normalize.method="vsn", pmcorrect.method="pmonly", summary.method="medianpolish")
		textlevel = leveltext("","down",textlevel)

		textlevel = leveltext("Writing normalized data","down",textlevel)
		write.exprs(eset.quantiles,file=paste(PATHdata,"/",analysis,'/eset_norm-quantiles.txt', sep=""))
		write.exprs(eset.rma,file=paste(PATHdata,"/",analysis,'/eset_norm-rma.txt', sep=""))
#		write.exprs(eset.vsn,file=paste(PATHdata,"/",analysis,'/eset_norm-vsn.txt', sep=""))

	}

	if (analysesList$general$qualityPlots){	
		textlevel = leveltext("Quality control plots after Normalization","keep",textlevel)

		if ((normalisierung == "quantiles") | (normalisierung == "")){
			pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/BOX_quantiles.pdf', sep=""))
			boxplot(data.frame(exprs(eset.quantiles)))
			dev.off()
			pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/DENS_quantiles.pdf', sep=""))	
				plot(density(log2(exprs(eset.quantiles[, 1]))), col=0, main="Quantiles Normalization")
				for(i in 1:dim(exprs(eset.quantiles))[2])
					{
					lines(density(log2(exprs(eset.quantiles[,i]))), col=i)
					}
			dev.off()
		}
		
		if (normalisierung == "rma" | normalisierung == ""){
			pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/BOX_rma.pdf', sep=""))
			boxplot(data.frame(exprs(eset.rma)))
			dev.off()
			pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/DENS_rma.pdf', sep=""))	
				plot(density(log2(exprs(eset.rma[, 1]))), col=0, main="Quantiles Normalization")
				for(i in 1:dim(exprs(eset.rma))[2])
					{
					lines(density(log2(exprs(eset.rma[,i]))), col=i)
					}
			dev.off()
		}
		
	#	if ((normalisierung == "vsn" || (normalisierung == ""))){
	#		pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/BOX_vsn.pdf', sep=""))
	#		boxplot(data.frame(exprs(eset.vsn)))
	#		dev.off()
	#		pdf(file=paste(PATHoutp,"/",analysis,'/Q-Plots/DENS_rma.pdf', sep=""))	
	#			plot(density(log2(exprs(eset.vsn[, 1]))), col=0, main="Quantiles Normalization")
	#			for(i in 1:dim(exprs(eset.vsn))[2])
	#				{
	#				lines(density(log2(exprs(eset.vsn[,i]))), col=i)
	#				}
	#		dev.off()
	#	}

	}

	textlevel = leveltext("No norm method specified, using RMA corrected quantiles. Alternatives are saved to disk","keep",textlevel)
	eset.norm=eset.rma

}else{
	textlevel = leveltext("Reading pre-normalized data from file","keep",textlevel)
	eset.norm=readExpressionSet(paste('eset_norm-',auswertung,'-',normalisierung,'.txt', sep="")) # reading already normalized data
}

	textlevel = leveltext("(Checking Subsets) - Currently disabled!","keep",1)

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



