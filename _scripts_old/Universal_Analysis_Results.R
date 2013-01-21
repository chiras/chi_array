
	cat("------------------------------------------------------------------------------------------------\n")	
	cat("End of Analyses, the following procedures have been conducted: \n")
	for (i in 1:length(analysis_log)){cat(paste("[",i,"] ", analysis_log[i], "\n",sep =""))}
	cat("------------------------------------------------------------------------------------------------\n")	
	cat("Your Files are in: \n")
	for (i in 1:length(dir_log)){cat(paste("[",i,"] ", file_log[i],": ",dir_log[i],"\n", sep =""))}
	cat("------------------------------------------------------------------------------------------------\n")	
#	cat("Further Analytical Tools are: \n")
#	for (i in 1:length(adv_log)){cat(paste("[",i,"] ", adv_log[i],"\n", sep =""))}
#	cat("------------------------------------------------------------------------------------------------\n")	
#
