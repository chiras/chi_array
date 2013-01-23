chi_array
=========

This script was created for special analysis, yet has been adapted over the years to serve in multiple analyses with various features. It has become a growning bioinformatical pipeline to analyze microarray chip data.

Please be aware that it has been designed for specific purposes needed in our research area. It may or may not, but will likely yield errors, when used under different conditions as ours.

(1) Valid chip data: Affymetrix HGU133plus2. The script may work with other chips, provided that you adapt annotation data accordingly

(2) Testing design: 
  - Single analysis: SampleA vs. SampleB
    -> Differential gene expression analysis
    Example: KnockDown Cells vs. Control Cells

  - Multiple analysis: (SampleA vs. SampleB) vs. (SampleC vs. SampleD)
    -> Comparison between differential expressions
    Example: KO/Control-differential expression in Chondrocytes vs. KO/Control-differential expression in Stem Cells

(3) valid folder structure: maintain the folders and files as downloaded

(4) use the Parameters.R.template file in ./Parameters, remove the .template extension and edit the Parameters according to your needs.

