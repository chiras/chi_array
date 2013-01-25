########################################################################################
# This is the parameter file for the chi-array script. 
# Duplicate this file for each separate analysis and rename it to "YOUR_ANALYSIS_NAME.R". Leave it in the ./Parameters folder
#
# Abbreviations: 
# 	- Y_A_NAME: Name of your analysis, no spaces or OS-shell invalid characters allowed
# 	- SA: Single Analysis (e.g. SampleTypeA vs. SampleTypeB)
# 	- MFA: Multifactorial Analysis == SA vs. SA (e.g. (SampleTypeA vs. SampleTypeB) vs. (SampleTypeC vs. SampleTypeD))
#
# This file is structured into blocks, which have to be edited to match your needs:
#	- General information of the analysis
#	- SA-analysis informtion block, has to be copied and pasted for each SA. Needs at least two for MFA
#	- Heatmap block: Heatmaps that you want: may be list of probesets, list of genes or pathways
# Stop editing below the Heatmap-block. There are some initial steps performed before the analytical scripts can be loaded. They need to be in this file.

########################################################################################
### GENERAL INFORMATION FOR THIS ANALYSIS

analysesList=list()
analysesList$general$name					=	"julia_koku_MSCvsOstdif"						# Analysis name: must match the subfolder name in ./Data/Y_A_NAME/
analysesList$general$phenodataFile			=	"phenoData_julia_koku_MSCvsOstdif.txt"			# How your phenodata file is called in ./Data/Y_A_NAME/
analysesList$general$sampleNames		=c("MSC905_ctrl","MSC918_ctrl","MSC920_ctrl","MSC932_ctrl","MSC933_ctrl","MSC905_koku","MSC918_koku","MSC920_koku","MSC932_koku","MSC933_koku","MSC918_ostdif_ctrl","MSC941_ostdif_ctrl","MSC946_ostdif_ctrl","MSC947_ostdif_ctrl","MSC949_ostdif_ctrl","MSC918_ostdif_koku","MSC941_ostdif_koku","MSC946_ostdif_koku","MSC947_ostdif_koku","MSC949_ostdif_koku")
																								# Names of the samples in the same order as in your phenodata file
analysesList$general$MFAplotVenn		= TRUE													# MFA only: Do you want a Venn Plot between the SAs (TRUE/FALSE)
analysesList$general$MFAplotHeatmap		= TRUE													# MFA only: Do you want a comparative Venn Plot between the SAs (TRUE/FALSE)
analysesList$general$scalename 			= "Treat-Ctrl"											#		-> Name of the HeatMap Scale (usual design of differential expression)
analysesList$general$use_norm 			= TRUE													# TRUE when you want to normalize, FALSE if you have already normalized data (needs eset file in folder ./Data/Y_A_NAME/)
analysesList$general$normalization 		= "rma"													# Which normalization should be applied or has been used in prenormalized data (one of: "rma","quantiles","vsn" or "" for all three in comparison)
analysesList$general$qualityPlots 		= TRUE	 												# 		-> Only with new normalization: You want to see Quality Plots?
### END: GENERAL INFORMATION BLOCK
########################################################################################

########################################################################################
### BEGIN: SA-BLOCK 1 (copy and paste this block for each separate SA, at least 2 blocks needed for MFA)
# SA information
tmpAnalysis=list()
tmpAnalysis$name 						= "KoKu_MSC"											# Name of this differential expression SA
tmpAnalysis$use_subset 					= c(1:5,6:10)											# Which samples in the order of the phenodata file are used in this SA?
tmpAnalysis$design 						= model.matrix(~ 0+factor(c(1,1,1,1,1,2,2,2,2,2)))		# Adapt the vector c(...) to match different test design categories (e.g. 1 for KnockDown, 2 for Control) in the order of the phenodata file
colnames(tmpAnalysis$design) 			= c("ctr", "KoKu")										# Names of these categories
tmpAnalysis$comparison					= "KoKu-ctr"											# Test design: how these categories should be compared (e.g. KnockDown-Control or Control-KnockDown). Defines the direction of diff expression.
tmpAnalysis$normalization 				= "rma" 												# name of the normalization algorithm (one of: "rma","quantiles","vsn")

# What steps to perform
tmpAnalysis$toDo$includeMFA	 			= TRUE													# Shall this SA be used in the overall MFA							
tmpAnalysis$toDo$includeSA	 			= FALSE													# Shall advanced SA analysis steps be performed, see below or restrict to steps necessary for MFA	
tmpAnalysis$toDo$analyzeGOstats 		= F														# 		- GO/KEGG enrichment analysis (TRUE/FALSE)
tmpAnalysis$toDo$plotPathHeatmaps 		= FALSE													# 		- (currently disabled) plot heatmaps of significant enriched pathways (TRUE/FALSE)
tmpAnalysis$toDo$parseGenes2Pathways	= T														# 		- yields a table with probsets expanded according to pathway designation -> multiple occurrences of sets likely (TRUE/FALSE)
tmpAnalysis$toDo$plotVenn				= FALSE													# 		- (currently disabled) plot a Venn diagram for this comparision only
tmpAnalysis$toDo$plotPathInteractions	= FALSE													# 		- (currently disabled) plot an interaction plot of pathways... errors due to bugs in the actual version of a necessary package
tmpAnalysis$toDo$Normalize 				= FALSE	 												# Currently disabled option, leave with FALSE

# Filtering options
tmpAnalysis$toDo$AFFX_filter 			= TRUE 													# Cut-Off at the maximum p-value matching AFFX-Probesets (TRUE/FALSE)
tmpAnalysis$toDo$IQR_filter				= FALSE 												# Cut-Off under special IQR conditions, see below (TRUE/FALSE)
tmpAnalysis$Filtering$IQR_intens_above 	= 100													#		- the intensity of a probe should be above XXX (e.g. 100) ...
tmpAnalysis$Filtering$IQR_intens_probes = 0.75 													#		  ...in at least XX (e.g. 75) percent of the samples
tmpAnalysis$Filtering$IQR_greater_than	= 1.0 													#		-  the interquartile range of log2–intensities should be at least X (1.00)
tmpAnalysis$toDo$logFC_filter 			= TRUE	 												# Cut-Off at the maximum p-value matching AFFX-Probesets (TRUE/FALSE)
tmpAnalysis$Filtering$logFC_threshold	= 0.5													#		- cut-off value

analysesList[[tmpAnalysis$name]]=tmpAnalysis
### END: SA-BLOCK 1
########################################################################################

########################################################################################
### BEGIN: SA-BLOCK 2 (copy and paste this block for each separate SA, at least 2 blocks needed for MFA)
# SA information
tmpAnalysis=list()
tmpAnalysis$name 						= "KoKu_Ostdif"											# Name of this differential expression SA
tmpAnalysis$use_subset 					= c(11:15,16:20)											# Which samples in the order of the phenodata file are used in this SA?
tmpAnalysis$design 						= model.matrix(~ 0+factor(c(1,1,1,1,1,2,2,2,2,2)))		# Adapt the vector c(...) to match different test design categories (e.g. 1 for KnockDown, 2 for Control) in the order of the phenodata file
colnames(tmpAnalysis$design) 			= c("ctr", "KoKu")										# Names of these categories
tmpAnalysis$comparison					= c("KoKu-ctr")												# Test design: how these categories should be compared (e.g. KnockDown-Control or Control-KnockDown). Defines the direction of diff expression.
tmpAnalysis$normalization 				= "rma" 												# name of the normalization algorithm (one of: "rma","quantiles","vsn")

# What steps to perform
tmpAnalysis$toDo$includeMFA	 			= TRUE													# Shall this SA be used in the overall MFA							
tmpAnalysis$toDo$includeSA	 			= FALSE													# Shall advanced SA analysis steps be performed, see below or restrict to steps necessary for MFA	
tmpAnalysis$toDo$analyzeGOstats 		= F														# 		- GO/KEGG enrichment analysis (TRUE/FALSE)
tmpAnalysis$toDo$plotPathHeatmaps 		= FALSE													# 		- (currently disabled) plot heatmaps of significant enriched pathways (TRUE/FALSE)
tmpAnalysis$toDo$parseGenes2Pathways	= T														# 		- yields a table with probsets expanded according to pathway designation -> multiple occurrences of sets likely (TRUE/FALSE)
tmpAnalysis$toDo$plotVenn				= FALSE													# 		- (currently disabled) plot a Venn diagram for this comparision only
tmpAnalysis$toDo$plotPathInteractions	= FALSE													# 		- (currently disabled) plot an interaction plot of pathways... errors due to bugs in the actual version of a necessary package
tmpAnalysis$toDo$Normalize 				= FALSE	 												# Currently disabled option, leave with FALSE

# Filtering options
tmpAnalysis$toDo$AFFX_filter 			= TRUE 													# Cut-Off at the maximum p-value matching AFFX-Probesets (TRUE/FALSE)
tmpAnalysis$toDo$IQR_filter				= FALSE 												# Cut-Off under special IQR conditions, see below (TRUE/FALSE)
tmpAnalysis$Filtering$IQR_intens_above 	= 100													#		- the intensity of a probe should be above XXX (e.g. 100) ...
tmpAnalysis$Filtering$IQR_intens_probes = 0.75 													#		  ...in at least XX (e.g. 75) percent of the samples
tmpAnalysis$Filtering$IQR_greater_than	= 1.0 													#		-  the interquartile range of log2–intensities should be at least X (1.00)
tmpAnalysis$toDo$logFC_filter 			= TRUE	 												# Cut-Off at the maximum p-value matching AFFX-Probesets (TRUE/FALSE)
tmpAnalysis$Filtering$logFC_threshold	= 0.5													#		- cut-off value

analysesList[[tmpAnalysis$name]]=tmpAnalysis
### END: SA-BLOCK 2
########################################################################################


########################################################################################
### BEGIN: HEATMAP BLOCK (copy and paste this Block for each analysis)
complist = mypaths = mygroups = myprobes = list()

# each list you compile here will be put in a separate heatmap

# PATHWAYS to look at (important: use the original KEGG 5 digits ID e.g. "04210" not "4210")
# you may delete or copy these two example pathways to any number of pathways, each will be put in a seperate heatmap
mypaths$Apoptosis_04210 					= c("04210")											# rename Apoptosis04210 to a unique name suiting the pathway, replace 04062 by KEGG pathway number of interest				
mypaths$ChemokineSignal_04062 				= c("04062")											# rename ChemokineSig_04062 to a unique name suiting the pathway, replace 04062 by KEGG pathway number of interest	

# GENES to look at
# you may delete or copy these two example lists to any number of lists, each will be put in a seperate heatmap
#groups$ossification 					= c("CDH11","FGF2","LEF1","SP7","VEGFA","WNT9A")		# rename ossification to a unique name suiting the list, replace vector c(...) by Gene-Symbols of interest
mygroups$others 						= c("CHI3L1","CILP","CILP2","CRTAC1","MSX1")			# rename others to a unique name suiting the list, replace vector c(...) by Gene-Symbols of interest

# PROBESETS to look at
# you may delete or copy these two example lists to any number of lists, each will be put in a seperate heatmap
myprobes$ADAMs_MMPs						= c("204475_at","205680_at","203877_at","203878_s_at")	# rename ADAMs_MMPs to a unique name suiting the list, replace vector c(...) by probeset-IDs of interest
myprobes$others 						= c("1552726_at","1553234_at","242823_at","230040_at")	# rename others to a unique name suiting the list, replace vector c(...) by probeset-IDs of interest

### END: HEATMAP BLOCK
########################################################################################



				



########################################################################################
########################################################################################
######### DANGER ZONE: STOP EDITING BELOW ##############################################
########################################################################################
########################################################################################
########################################################################################
devel_mode=F		# shall this script be run in developmental mode? 
					#		TRUE will result in truncated results but faster runtime
########################################################################################
######### READ & SETTING ENVIRONMENT ###################################################
########################################################################################

tmplog=c()
tmplog=append(tmplog,"\n\n Setting global parameters\n")

tmplog=append(tmplog,"Setting up working directories")
frame_files <- lapply(sys.frames(), function(x) x$ofile)
if(!length(frame_files)>0){stop('\n\n\nERROR\n!!!Script musst be called WITHIN R with >source("ABSOLUTE_PATH_TO_PARAMETER_FILE")!!!\n\n')}

frame_files <- Filter(Negate(is.null), frame_files)
PATH <- dirname(frame_files[[length(frame_files)]])

PATHbase = paste(PATH,"/../",sep="")
PATHscpt = paste(PATH,"/../_scripts",sep="")
PATHress = paste(PATH,"/../_ressources",sep="")
PATHdata = paste(PATH,"/../data",sep="")
PATHoutp = paste(PATH,"/../output",sep="")


if (file.exists(paste(PATHoutp,"/",analysesList$general$name,sep=""))){
	tmplog=append(tmplog,"Results directory already exists, overwriting files!")
} else {
	tmplog=append(tmplog,"Setting up results directory")
    dir.create(paste(PATHoutp,"/",analysesList$general$name,sep=""))
    dir.create(paste(PATHoutp,"/",analysesList$general$name,"/Q-Plots",sep=""))
    dir.create(paste(PATHoutp,"/",analysesList$general$name,"/SA-Plots",sep=""))
   # dir.create("R-Output/KEGGdata")
    dir.create(paste(PATHoutp,"/",analysesList$general$name,"/MFA-Plots",sep=""))
}

tmplog=append(tmplog,"Checking local Environment")
direrr = c()
direrr$general = "Directory/file check:\n"
direrr$script = ifelse(file.exists(PATHscpt),"OK","[ERR] BASEDIR/_scripts missing\n")
direrr$ress = ifelse(file.exists(PATHress),"OK","[ERR] BASEDIR/_ressources folder missing\n")
direrr$data = ifelse(file.exists(PATHdata),"OK","[ERR] BASEDIR/data folder missing\n")

if (direrr$data == "OK"){
	tmplog=append(tmplog,"Setting local directory")
	setwd(paste(PATHdata,"/", analysesList$general$name, sep="")) #setting working directory
	
	filename_norm_data= paste("eset_norm-", analysesList$general$name,"-", analysesList$general$normalization,".txt", sep="")
	
	if(!analysesList$general$use_norm){
	    direrr$eset_norm = ifelse(file.exists(paste(PATHdata,"/",analysesList$general$name,"/",filename_norm_data,sep="")),"OK",paste("[ERR] Normalized data BASEDIR/data/",analysesList$general$name,"/",filename_norm_data," missing",sep=""))
	}
    direrr$phenodata = ifelse(file.exists(paste(PATHdata,"/",analysesList$general$name,"/",analysesList$general$phenodataFile,sep="")),"OK",paste("[ERR] Phenodata BASEDIR/data/",analysesList$general$name,"/",analysesList$general$phenodataFile," missing",sep=""))
}

tmplog=append(tmplog,paste("Scripts folder: ", direrr$script, sep=""))
tmplog=append(tmplog,paste("Ressources folder: ", direrr$ress, sep=""))
tmplog=append(tmplog,paste("Data folder: ", direrr$data, sep=""))
if(!analysesList$general$use_norm){
	tmplog=append(tmplog,paste("Normalized data: ",direrr$eset_norm, sep=""))
}
tmplog=append(tmplog,paste("Phenodata: ",direrr$phenodata, sep=""))

direrr$results =  ifelse(direrr$script == direrr$ress && direrr$ress == direrr$data && direrr$ress == direrr$phenodata && direrr$ress == "OK", "everything fine!", "-> please address the errors above!\n")
if(!analysesList$general$use_norm && direrr$results == "everything fine!"){
	direrr$results =  ifelse(direrr$eset_norm == "OK", "everything fine!", "-> please address the errors above!\n")
}

tmplog=append(tmplog,paste(direrr$result,"\n",sep=""))

if (direrr$results == "everything fine!"){
	source(paste(PATHscpt,'/chi_workflow.R', sep=""))
}else{
	cat (tmplog, sep="\n")
}
