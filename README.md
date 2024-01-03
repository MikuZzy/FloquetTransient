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

Most R scripts include commands of setting working directory. One may change the directory *"~/R/novogene_trimmed/"* to the directory where this repo is cloned.

Requirements for external libraries are specified in each script. In general, the `quantreg` v5.97 and the `parallel` package are required for running the main algorithm, and plots are created with `ggplot2` v3.4.2 and `igraph` v1.4.2 packages. Source codes can be executed without any further installations. Codes has been tested and successfully executed with `R` v4.1.1 on Northwestern University's Quest High-Performance Computing Cluster (with multiple CPUs) and `Matlab` R2022b on Win10 system.

The experimental data are included as .csv files in the repo in addition to the GEO database, thus the expected outputs can be viewed in the preprint. The algorithm is computation-intensive - a typical runtime for experimental and synthetic data is ~3min per gene with a single CPU - and thus we recommend running it on an HPC cluster. 
