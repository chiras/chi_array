chi_array
=========

This script was created for a special purpose, yet has been adapted over the last years to serve in very different situtations with various datasets. The features are constantly accumulating; it has become a growing bioinformatical pipeline to analyze microarray chip data.

Please be aware that it has been designed for specific purposes needed in our research group. It may or may not work, but will likely yield errors, when used under different conditions as ours (see "How to use it").

What it does:
--------------

Analyses:
- Various normalization methods implemented in Bioconductor
- Various quality and data constitution plots 
- Linear model statistics
- Filtering via AFFX cut, IQR and/or logFC cutoff
- Calculating overlaps and Venn plots in comparison between differential expressions
- Heatmaps for specific probesets, pathways or genes
- GO/KEGG enrichment analysis

Further features:
- Automatically detecting/checking/setting the global enviroment for the analysis
- Logging analytical steps, used parameters and errors/warnings into files
- Automatically downloading and loading required packages for R and Bioconductor
- Loading custom analytical functions
- Validating input parameters


How to use it:
--------------

1. Valid chip data: Affymetrix HGU133plus2.0 only. The script may work with other chips, provided that you adapt annotation data accordingly
2. Testing design: 
  * Single analysis:
   SampleA vs. SampleB  
   -> Differential gene expression analysis  
   Example: KnockDown Cells vs. Control Cells  
  * Multiple analysis:
   (SampleA vs. SampleB) vs. (SampleC vs. SampleD)  
   -> Comparison between differential expressions  
   Example: KO/Control-differential expression in Chondrocytes vs. KO/Control-differential expression in Stem Cells
3. Valid folder structure: maintain the folders and files as downloaded
4. Download the Gene sets U133plus2 .Rdata file at http://compbio.med.harvard.edu/Supplements/PNAS05.html
5. use the template.R file in ./Parameters: rename according to your analysis and edit the Parameters inside according to your needs.
6. Store your .CEL files in the folder ./data within a subfolder named according to the setting "analysesList$general$name" in the parameter file.
7. Be sure to have a working internet connection to able to download needed libraries
8. Call the script WITHIN R using 'source("PATH_TO_YOUR_PARAMETER_FILE.R")' or by using the "load script" option in the GUI menu
  * Mac only option: you may simply drag and drop your parameter.R file onto the Mac-Tool. It will start your analysis and check for script-updates at the git-repos.
9. Please acknowledge my and others work by reading the CONTRIB.txt and citing the relevant articles.


Tested on:
--------------
so far only on a couple of Macs

To do:
--------------
many things, soon coming: input parameter checking and providing test data

--------------
All scripts written by Alexander Keller, DNA Analytics Core Facility, Biocenter, University of Würzburg