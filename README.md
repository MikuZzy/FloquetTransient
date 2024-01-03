# FloquetTransient

This is the code repository for the research paper *"Modeling Transient Changes in Circadian Rhythms"* currently on arXiv *(https://doi.org/10.48550/arXiv.2304.07412)*. To generate analysis results and figures, run codes in the following orders:

 - 00_Preprocessing/\*.r
  
 - 00_Reproducibility/\*.m, 00_Reproducibility/\*.r
  
 - 00_SyntheticData/\*.r
  
 - 01novogene_job.r
  
 - 02_SimulatedFNs/\*.r
  
 - 03_Residuals/\*.r
  
 - 04_Phase/\*.r
  
 - 05_GRN/\*.r
  
 - 06_MultiLambda/\*.r
  
 - 07_Additional/\*.r

Most R scripts include commands of setting working directory. One may change the directory *"~/R/novogene_trimmed/"* to the directory where this repo is cloned. Requirements for external libraries are specified in each script. In general, the `quantreg` and the `parallel` package are required for running the main algorithm, and plots are created with `ggplot2` and `igraph` packages. 