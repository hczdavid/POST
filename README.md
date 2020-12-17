# POST
Phylogeny-Guided Microbiome OTU-Specific Association Test

### Overview

The POST package is the implementation of our recent method Phylogeny-Guided Microbiome OTU-Specific Association Test. POST boosts the testing power by adaptively borrowing information from phylogenetically close OTUs of the target OTU. Whether or not borrowing information or the amount of information from the neighboring OTUs is data adaptive and supervised by phylogenetic distance and the outcome variable. POST is built on a kernel machine regression framework and inherited the advantages including flexibly model complex microbiome effects (e.g., effects from opposite direction), easily adjust for covariates, and accommodate both continuous and binary outcomes. POST extends the current global kernel tests \citep{zhao2015testing, wu2016adaptive, koh2017powerful} to OTU-level test using a OTU-specific kernel.


### Dependencies

The following packages are required for functions and examples in the POST package: 

GUniFrac, CompQuadForm, ACAT. 

We adapt two functions (CKAT.bin and CKAT.con) from CKAT pacakge calculate p-value (Debashis Ghosh and Xiang Zhan. "CKAT software").


### Installation

```{r}
devtools::install_github("hczdavid/POST")
```











