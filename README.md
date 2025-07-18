# Introduction
This package will be a collection of tools that are routinely used by Texas A&M University Rose Breeding Program for analysis of genetic data from the WagRhSNP68k SNP array
    
Tools in this package were used in this paper:       
[Rose Rosette Disease Resistance Loci Detected in Two Interconnected Tetraploid Garden Rose Populations](https://www.frontiersin.org/articles/10.3389/fpls.2022.916231/full)     
    
    
**IMPORTANT:** Must have internet connection for some of the functions as some functions need annotation files that are hosted online.    
   
To download this package in R or Rstudio, use the following lines of code:          
     
		 
```
install.packages('devtools')
devtools::install_github("jeekinlau/RoseArrayTools")
```
      
     
For R-package to run correctly, you must also install the R-packages: "data.table", "splitstackshape"     	 
The two dependencies below should install automatically using the two lines of code above but in case it does not install them manually.    
    
```
install.packages('data.table')
install.packages('splitstackshape')
```       
    
Use the link to the Vignette below for the functions included in this package.      
*[Link to Vignette](https://jeekinlau.github.io/RoseArrayTools/RoseArrayTools_Vignette.html)

# How to cite
Lau, Jeekin, Ellen L. Young, Sara Collins, Mark T. Windham, Patricia E. Klein, David H. Byrne, and Oscar Riera-Lizarazu. “Rose Rosette Disease Resistance Loci Detected in Two Interconnected Tetraploid Garden Rose Populations.” Frontiers in Plant Science 13 (2022). https://doi.org/10.3389/fpls.2022.916231.


# Acknowledgements
This work was supported by grants from the United States Department of Agriculture, Specialty Crops Research Initiative grant “Combating Rose Rosette Disease: Short and Long Term Approaches” Grant No. 2014-51181-22644 and “Tools for Genomic-Assisted Breeding in Polyploids” Grant No. 2020-51181-32156, the Basye Endowment in Rose Genetics and the American Rose Society.
