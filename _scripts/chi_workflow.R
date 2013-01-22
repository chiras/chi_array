textlevel = leveltext("Loading Libraries","keep",textlevel)
source(paste(PATHscpt,"/chi_libraries.R",sep=""))

#cat("Looking up script version\n")
#source("../_scripts/chi_versioning.R")					
		
textlevel = leveltext("Loading analysis specific functions","down",textlevel)
source(paste(PATHscpt,"/chi_functions.R",sep=""))

		
textlevel = leveltext("Checking input parameters","down",textlevel)
source(paste(PATHscpt,"/chi_settings.R",sep=""))

textlevel = leveltext("Retrieving data","down",textlevel)
source(paste(PATHscpt,"/chi_normalization.R",sep=""))

textlevel = leveltext("Starting analysis","down",textlevel)
source(paste(PATHscpt,"/chi_analysis.R",sep=""))

##source("../_scripts/Universal_Analysis_MFA.R")
##source("../_scripts/Universal_Analysis_Results.R")