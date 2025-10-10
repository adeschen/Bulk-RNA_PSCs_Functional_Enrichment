# Functional Enrichment on PCSs Differential Expression from Mucciolo G. et al 2024

This capsule contains the functional enrichment analysis conducted on the 
bulk RNA-seq PSCs from 
[Mucciolo et al. 2024](https://doi.org/10.1016/j.ccell.2023.12.002).
The capsule also includes the code to generate an enrichment map with the 
selected terms of interest related to axon guidance obtained during the 
functional enrichment step.

More specifically, the data used in this analysis are from [Mucciolo et al. 2024](https://doi.org/10.1016/j.ccell.2023.12.002):

- *Table S1.* RNA-sequencing analysis of control-treated or TGF-beta-treated PSCs - differential expression analysis (related to Figure 1).
- *Table S2.* RNA-sequencing analysis of control-treated or CM-treated PSCs with or without ERBBi - differential expression analysis (related to Figures 2 and 3).

## Authors

[J&eacute;r&eacute;my Nigri](https://orcid.org/0000-0003-1358-1863), [Astrid Desch&ecirc;nes](https://orcid.org/0000-0001-7846-6749) and [David A Tuveson](https://tuvesonlab.labsites.cshl.edu/)

### Code, Data, Results

All associated code is provided in the *code* directory.

The *data* directory contains all necessary data.

Output files and plot will be saved in the *results* directory. The *results* 
directory must be created before running the code.

### How to replicate analysis

This code has been used to create the Code Ocean capsule.

The R software v4.5.1 was used for the analysis.

The specific R package versions used for the analysis are:

* CRAN stringr v1.5.2 
* CRAN gprofiler2 v0.2.3 
* CRAN ggplot2 v4.0.0
* CRAN scatterpie v0.2.6
* CRAN ggtangle v0.0.7
* CRAN ggrepel v0.9.6
* CRAN ggnetwork v0.5.14
* CRAN igraph v2.1.4
* CRAN devtools v2.4.5
* Bioconductor enrichViewNet v1.7.3


### License ###

This package and the underlying  code are distributed under 
the **GNU GENERAL PUBLIC LICENSE**.

### Maintainer

[Astrid Desch&ecirc;nes](http://ca.linkedin.com/in/astriddeschenes "Astrid Desch&ecirc;nes")

### Questions

[Please contact us](https://github.com/adeschen/Bulk-RNA_PSCs_Functional_Enrichment/issues) for related questions.

Thanks!