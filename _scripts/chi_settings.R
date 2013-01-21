
########################################################################################
######## READ & SETTING ENVIRONMENT ####################################################
########################################################################################


analysis_log=c()
dir_log=c()
file_log=c()
adv_log=c()


if (file.exists("R-Output")){
} else {
    cat("Setting up results directory\n")
    analysis_log=c(analysis_log,"Setting up results directory")
    
    dir.create("R-Output")
    dir.create("R-Output/Plots")
    dir.create("R-Output/KEGGdata")
    dir.create("R-Output/MFA-Plots")
}


cat("(Checking Analysis Parameters) - Currently disabled!\n")

# checking whether all parameters have been correctly defined
#continue = c();

#if (!file.exists(phenodataFile)) 									{continue = c(continue, '- No/Bad Phenodata File specified\n')}
#if (use_IQR_filter){
#	if (!is.numeric(IQR_intens_above))								{continue = c(continue, '- the intensity of a probe should be# specified for IQR filter\n')}
#	if (!is.numeric(IQR_intens_probes))								{continue = c(continue, '- the percentage of a probes should be specified for IQR filter\n')}
#	if (IQR_intens_probes<0 | IQR_intens_probes>1)					{continue = c(continue, '- the percentage of a probes should be within 0 and 1\n')}
#	if (!is.numeric(IQR_greater_than))								{continue = c(continue, '- the minimum IQR should be specified for IQR filter\n')}
#}


#if (!use_norm){
#	if ((normalization != "rma" & normalization != "vsn" & normalization != "quantiles" ) & use_norm == F) {continue = c(continue, '- Please refer to used Normalization method\n')}
#	filename_norm_data= paste("eset_norm-", analysis,"-", normalization,".txt", sep="")
	
#	if (!file.exists(filename_norm_data)) 									{continue = c(continue, paste('- File with normalized data not found, should be: ', filename_norm_data,"\n", sep=""))}
#}