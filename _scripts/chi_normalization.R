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

cat("Loading Chip Annotation Data\n")
load("../_scripts/GenesetsU133plus2.RData")													# Chip Annotation data


if (use_norm){ 																			# when new normalization is used
	
########################################################################################
######## READ RAW DATA #################################################################
########################################################################################

	cat("Normalization: Loading Raw Data\n")
	pd <- read.AnnotatedDataFrame(phenodataFile, header = TRUE, row.names = 1, sep="") 	# read Phenotype data
	pData(pd) 																			# Map Phenotype Data
	rawAffyData <- ReadAffy(filenames= rownames(pData(pd)), phenoData=pd) 				# read raw data
	
	analysis_log=c(analysis_log,"Reading Raw Data")
	
	# look at genes, probes, probesets
	length(featureNames(rawAffyData))	# gene (probesets)
	length(probeNames(rawAffyData))		# probes -> 11 probes pro probeset

########################################################################################
######## QUALITY CONTROL ###############################################################
########################################################################################

	if (qualityPlots){	

	cat("Normalization: Quality Control\n")
	analysis_log=c(analysis_log,"Quality Control Plots")
	dir_log=c(dir_log,"R-Output/Plots/")
	file_log=c(file_log,"Quality Plots")

	# plot boxplots
	cat(noquote("Normalization: Quality Control: Generating Boxplots\n"))
	pdf(file=paste('R-Output/Plots/box_nonorm.pdf', sep=""))
		boxplot(rawAffyData)
	dev.off()

	# plot densities
	cat(noquote("Normalization: Quality Control: Generating Density Plots\n"))
	pdf(file=paste('R-Output/Plots/dens_nonorm.pdf', sep=""))	
	plot(density(log2(intensity(rawAffyData[, 1]))), col =0, main="Densities", ylim=c(0,1))
	for(i in 1:dim(intensity(rawAffyData))[2])
	{
		lines(density(log2(intensity(rawAffyData[,i]))), col = i)
	}
	dev.off()

	cat(noquote("Normalization: Quality Control: MA-Plots\n"))
	pdf(file=paste('R-Output/Plots/MA_nonorm%03d.pdf', sep=""), onefile=FALSE) 	
	MAplot(rawAffyData)
	dev.off()

	# look at RNA degradation
	cat(noquote("Normalization: Quality Control: RNA Degradation Plots\n"))
	pdf(file=paste('R-Output/Plots/rna-deg_nonorm.pdf', sep=""))	
	deg <- AffyRNAdeg(rawAffyData, log.it=FALSE)
	plotAffyRNAdeg(deg)
	dev.off()

	}

########################################################################################
######## Normalization #################################################################
########################################################################################

	# use expresso and summarization for normalization
	cat(noquote("Normalization: Process\n"))
	if (normalisierung == "quantiles"){
	cat(noquote("Normalization: Process: Quantiles\n"))
		eset.norm <- expresso(rawAffyData, bg.correct=FALSE, normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		write.exprs(eset.norm,file=paste('R-Output/eset_norm-',normalisierung,'.txt', sep=""))
		analysis_log=c(analysis_log,"Normalization: Quantiles")
		dir_log=c(dir_log,paste('R-Output/eset_norm-',normalisierung,'.txt', sep=""))
		file_log=c(file_log,"Quantiles Normalized Data")
	}else if (normalisierung == "rma"){
	cat(noquote("Normalization: Process: RMA\n"))
		eset.norm <- expresso(rawAffyData,bgcorrect.method="rma", normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		write.exprs(eset.norm,file=paste('R-Output/eset_norm-',normalisierung,'.txt', sep=""))
		analysis_log=c(analysis_log,"Normalization: RMA")
		dir_log=c(dir_log,paste('R-Output/eset_norm-',normalisierung,'.txt', sep=""))
		file_log=c(file_log,"RMA Normalized Data")
	}else if (normalisierung == "vsn"){
#	cat(noquote("Normalization: Process: VSN\n"))
#		eset.norm <- expresso(rawAffyData, bg.correct=FALSE, normalize.method="vsn", pmcorrect.method="pmonly", summary.method="medianpolish")
#		write.exprs(eset.norm,file=paste('R-Output/eset_norm-',normalisierung,'.txt', sep=""))
#		analysis_log=c(analysis_log,"Normalization: VSN")
#		dir_log=c(dir_log,paste('R-Output/eset_norm-',normalisierung,'.txt', sep=""))
#		file_log=c(file_log,"VSN Normalized Data")
	}else{
	cat(noquote("Normalization: Process: Quantiles/RMA/VSN\n"))
		eset.quantiles <- expresso(rawAffyData, bg.correct=FALSE, normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
		eset.rma <- expresso(rawAffyData,bgcorrect.method="rma", normalize.method="quantiles", pmcorrect.method="pmonly", summary.method="medianpolish")
#		eset.vsn <- expresso(rawAffyData, bg.correct=FALSE, normalize.method="vsn", pmcorrect.method="pmonly", summary.method="medianpolish")
		analysis_log=c(analysis_log,"Normalization: Quantiles/RMA/VSN")

		# look at normalise data
		par(mfrow=c(1,3))
		pdf(file=paste('R-Output/Plots/box_norm%03d.pdf', sep=""), onefile=FALSE) 	
			boxplot(data.frame(exprs(eset.rma)))
			boxplot(data.frame(exprs(eset.quantiles)))
	#		boxplot(data.frame(exprs(eset.vsn)))
		dev.off()

		# look at density plots

		pdf(file=paste('R-Output/Plots/dens_norm%03d.pdf', sep=""), onefile=FALSE) 	
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
		write.exprs(eset.quantiles,file=paste('R-Output/eset_norm-quantiles.txt', sep=""))
		write.exprs(eset.rma,file=paste('R-Output/eset_norm-rma.txt', sep=""))

		dir_log=c(dir_log,paste('R-Output/eset_norm-quantiles.txt', sep=""))
		file_log=c(file_log,"Quantiles Normalized Data")
		dir_log=c(dir_log,paste('R-Output/eset_norm-rma.txt', sep=""))
		file_log=c(file_log,"RMA Normalized Data")
		dir_log=c(dir_log,paste('R-Output/eset_norm-vsn.txt', sep=""))
		file_log=c(file_log,"VSN Normalized Data")


		cat(noquote("Normalization: No method specified, using RMA. Alternates are saved to disk\n"))
		eset.norm=eset.rma
	}

	if (qualityPlots){	
	cat(noquote("Normalization: Quality Control after Normalization\n"))

	# plot boxplots
	cat(noquote("Normalization: Quality Control after Normalization: Generating Boxplots\n"))
	pdf(file=paste('R-Output/Plots/box_norm.pdf', sep=""))
		boxplot(data.frame(exprs(eset.norm)))
	dev.off()

	# plot densities
	cat(noquote("Normalization: Quality Control after Normalization: Generating Density Plots\n"))
	pdf(file=paste('R-Output/Plots/dens_norm.pdf', sep=""))	
		plot(density(log2(exprs(eset.norm[, 1]))), col=0, main="Normalization")
		for(i in 1:dim(exprs(eset.norm))[2])
		{
			lines(density(log2(exprs(eset.norm[,i]))), col=i)

		}
	dev.off()
	}
	
}else{
	cat(noquote("Reading Normalized Data from file\n"))
	eset.norm=readExpressionSet(paste('eset_norm-',auswertung,'-',normalisierung,'.txt', sep="")) # reading already normalized data
	analysis_log=c(analysis_log,"Reading Normalized Data")
}

	cat(noquote("(Checking Subsets) - Currently disabled!\n"))

#if (length(use_subset)==length(sampleNames(eset.norm))){
	
#}else{
#	cat(noquote("Limiting Data to subset, if specified\n"))
#	eset.norm = eset.norm[,use_subset] 
#	analysis_log=c(analysis_log,"Restricting to subset")
#}

cat(noquote("Renaming Samples\n"))
for (zz in 1:length(sampleNames(eset.norm))){
	 cat(paste ("-",sampleNames(eset.norm)[zz], "-->", analysesList$general$sampleName[zz],"\n"))
	}
# renaming datasets	
sampleNames(eset.norm) <- analysesList$general$sampleNames 	
analysis_log=c(analysis_log,"Renaming Samples")



