# clustlasso

This repository hosts the **clustlasso** R package supporting the paper [Interpreting k-mer based signatures for antibiotic resistance prediction](https://academic.oup.com/gigascience) (Jaillard *et al.*, 2020).

## Description

The **clustlasso** package is a generic package implementing both the lasso and cluster-lasso strategies described in the paper [Interpreting k-mer based signatures for antibiotic resistance prediction](https://academic.oup.com/gigascience).

While this paper focuses on the task of building predictive models of anti-microbial resistance (AMR) from whole-genome sequences and k-mers, the **clustlasso** package is generic and can be applied to any other applications, not necessarily involving k-mers nor AMR phenotypes.

 A step-by-step procedure describing its integration with [DBGWAS](https://gitlab.com/leoisl/dbgwas) is available on another dedicated GitLab repository: [clustlasso-dbgwas-integration](https://gitlab.com/biomerieux-data-science/clustlasso-dbgwas-integration).

## Installation 
The **clustlasso** package can be installed as any standard R package (as described for instance in the [R documentation](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages)).
This will typically involve typing the following command from your GitLab workspace :  
 ` install.package("clustlasso", type = "source", repos = NULL)`
Once the package is installed, it can be loaded using : 
`library(clustlasso)`

## Usage

A [vignette](https://gitlab.com/biomerieux-data-science/clustlasso/-/blob/master/vignettes/vignette.pdf) included in the package describes in detail a typical usage of the package on a real dataset related to *Neisseria gonorrhoeae* resistance to cefixime [[Eyre, 2017](https://pubmed.ncbi.nlm.nih.gov/28333355/)]. 

## Support

Feel free to contact [Magali Jaillard-Dancette](mailto:magali.jaillard-dancette@biomerieux.com) or [Pierre Mahé](mailto:pierre.mahe@biomerieux.com) for further information or to report any issue.


 

