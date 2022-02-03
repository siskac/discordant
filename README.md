discordant
======

R package for determining differential correlation of molecular 
feature pairs from -omics data using mixture models.

### Information

Method paper published in 
[Bioinformatics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5006287/)

>Siska C, Bowler R, Kechris K. The discordant method: a novel approach for 
>differential correlation [published correction appears in Bioinformatics. 2017 
>Jan 1;33(1):150]. Bioinformatics. 2016;32(5):690-696. 
>doi:10.1093/bioinformatics/btv633

Software paper published in 
[BMC Research Notes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5244536/)

>Siska C, Kechris K. Differential correlation for sequencing data. BMC Res 
>Notes. 2017;10(1):54. Published 2017 Jan 19. doi:10.1186/s13104-016-2331-9

### Installation

Install via Bioconductor:

    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("discordant")

Install via Github:

    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("siskac/discordant")

### Bug Reports

Report bugs as issues on the [GitHub repository new
issue](https://github.com/siskac/discordant/issues/new)