logfile = paste(PATHoutp,"/",analysis,"/log_",analysis,".log",sep="")
parfile = paste(PATHoutp,"/",analysis,"/log_",analysis,".par",sep="")
if (!devel_mode){errfile = paste(PATHoutp,"/",analysis,"/log_",analysis,".err",sep="")}

leveltext <- function (text,levelchange="keep",currlevel){
	currlevel = ifelse(levelchange=="up", currlevel+1, currlevel)
	currlevel = ifelse(levelchange=="down", currlevel-1, currlevel)
	if (text != ""){
		if (currlevel != 0){cat("  ")}
		cat(paste(c(rep("-",currlevel)," ",text,"\n"), sep="", collapse=""))
		cat(paste(c(rep("-",currlevel)," ",text,"\n"), sep="", collapse=""), file = logfile, append = T)
	}
	return(currlevel)
}
tmplog=append(tmplog,"Setting log files")


time <- Sys.time()

if (!devel_mode){
	cat(paste("# Errors for run at: ",time,"\n",sep=""), file = errfile, append = F)
	errconn = file(errfile, open="wt")
	sink(errconn, type="message")
}

if (sum(installed.packages()[,1]=="marray")==0) {install.packages("marray", dependencies = TRUE, quiet = TRUE)}
library(marray)

cat(paste("# Chi analysis started at ",time,"\n",sep=""), file = logfile, append = F)
if (devel_mode){cat("# !!!DEVELOPMENTAL MODE!!!\n# !!!Results are truncated to speed up runs.\n# !!!Do NOT use these results!\n\n" , file = logfile, append = T)}

cat(paste("# Parameters for run at: ",time,"\n",sep=""), file = parfile, append = F)
write.list(analysesList,filename=parfile)

for (i in 1:length(tmplog)){
	level = "keep"
	if (i==3){level = "up"}
	if (i==length(tmplog)){level = "down"}
	textlevel = leveltext(tmplog[i],level,textlevel)
}
