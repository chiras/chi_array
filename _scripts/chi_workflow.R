textlevel = leveltext("Loading Libraries","keep",0)
source(paste(PATHscpt,"/chi_libraries.R",sep=""))

#cat("Looking up script version\n")
#source("../_scripts/chi_versioning.R")					
		
textlevel = leveltext("Loading analysis specific functions","keep",0)
source(paste(PATHscpt,"/chi_functions.R",sep=""))

		
textlevel = leveltext("Checking input parameters","keep",0)
source(paste(PATHscpt,"/chi_settings.R",sep=""))

textlevel = leveltext("Retrieving data","keep",0)
source(paste(PATHscpt,"/chi_normalization.R",sep=""))

textlevel = leveltext("Starting analysis","keep",0)
source(paste(PATHscpt,"/chi_analysis.R",sep=""))

##source("../_scripts/Universal_Analysis_MFA.R")
##source("../_scripts/Universal_Analysis_Results.R")

textlevel = leveltext("\n ******************************************","keep",0)
textlevel = leveltext("* Analysis completed","keep",textlevel)
textlevel = leveltext(paste("* You can find your results in:\n * ",PATHoutp,"/",analysis,"/",sep=""),"keep",textlevel)
textlevel = leveltext("******************************************","keep",textlevel)

warnings(file=errfile, append=T)
sink()