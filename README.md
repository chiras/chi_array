chi_array
=========

This script was created for special analysis, yet has been adapted over the last years to serve in very differnet questions with various features. It has become a growning bioinformatical pipeline to analyze microarray chip data.

Please be aware that it has been designed for specific purposes needed in our research area. It may or may not, but will likely yield errors, when used under different conditions as ours (see "How to use it").

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

(1) Valid chip data: Affymetrix HGU133plus2.0 only. The script may work with other chips, provided that you adapt annotation data accordingly

(2) Testing design: 

    - Single analysis: SampleA vs. SampleB
    -> Differential gene expression analysis
    Example: KnockDown Cells vs. Control Cells

    - Multiple analysis: (SampleA vs. SampleB) vs. (SampleC vs. SampleD)
    -> Comparison between differential expressions
    Example: KO/Control-differential expression in Chondrocytes vs. KO/Control-differential expression in Stem Cells

(3) valid folder structure: maintain the folders and files as downloaded

(4) use the Parameters.R.template file in ./Parameters, remove the .template extension and edit the Parameters according to your needs.

(5) Be sure to have a working internet connection to able to download needed libraries

(6) Call the script WITHIN R using source("PATH_TO_YOUR_PARAMETER_FILE") or by using the "load script" option in the GUI menu

(7) Please acknowledge my work in resulting publications by mentioning the github ressource URL along with my author name "Alexander Keller" (e.g. as an electronic manual or online ressource).
