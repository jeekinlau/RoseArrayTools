# Introduction
This package will be a collection of tools that are routinely used by Texas A&M University Rose Breeding Program for analysis of genetic data from the WagRhSNP68k SNP array
     
To download this package in R or Rstudio, use the following lines of code:          
 
The R-package "devtools" is used to install "RoseArrayTools" from Github        
		 
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
