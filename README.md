
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mrQC

<!-- badges: start -->

<!-- badges: end -->

mrQC is an R package for checking the accuracy of meta- and summary-data
for outcome traits prior to their use in Mendelian randomisation
analyses. For example, the package provides tools for checking that the
reported effect allele and effect allele frequency columns are correct.
It also checks for possible errors in the summary data that might
introduce bias into downstream Mendelian randomisation analyses.

## Installation

To install the latest version of mrQC, perform as normal:

``` r
install.packages("devtools")
devtools::install_github("MRCIEU/mrQC")
```

## General overview

This package exploits three groups of SNPs: 1) **the MAF 1KG reference
set**. This a set of 2297 SNPs that have the same minor allele across
the 1000 genomes super populations and that have a minor allele
frequency between 0.1 and 0.3. These SNPs are used to check for allele
frequency errors; 2) **the GWAS catalog top hits**. These are SNPs in
the GWAS catalog that are strongly associated with the outcome trait of
interest; and 3) the **exposure SNPs** or the genetic instruments for
the exposure of interest.

The package performs 4 types of checks:

1)  confirm the identity of the effect allele column. This step relies
    on comparison of GWAS catalog top hits.  
2)  confirm the identity of the effect allele frequency column. This
    step relies on comparison of the MAF reference set.  
3)  assess the ancestral origins of a study. This step relies on
    comparison of the MAF reference set.  
4)  for case-control studies, confirm that the provided effect sizes are
    log odds ratios or identify potential summary data errors. This step
    can use any of the SNP sets.

## Example. Extract and format the target data

The first step is to make a list of SNP rsids. Let’s assume out outcome
trait of interest is glioma. Let’s identify the top GWAS hits for glioma
in the GWAS catalog. We will search on both the reported trait name as
well as the trait experimental factor ontology (EFO). Note that not all
studies in the GWAS catalog have up-to date EFO annotations. It is
therefore advisable to search on both efo and reported trait, to
maximise the number of retrieved hits.

``` r
library(mrQC)
snplist1<-make_snplist(efo="glioma",ref1000G_superops=FALSE) #search on EFO
snplist2<-make_snplist(trait="glioma",ref1000G_superops=FALSE) #search on EFO
snplist<-unique(c(snplist1,snplist2))
head(snplist)
length(snplist)
```

We see that searching on the glioma EFO retrieves 42 SNPs, searching on
glioma as the reported trait retrieves 37 SNPs, with 54 non duplicated
SNPs across the two lists. We could also have searched on the glioma EFO
ID, which retrieves the same SNPs as searching on the “glioma” EFO.
e.g.,

``` r
library(mrQC)
snplist3<-make_snplist(efo_id="EFO_0005543",ref1000G_superops=FALSE) 
length(snplist)
```

In the above examples we set ref1000G\_superops to FALSE. We now set
this to TRUE, which will result in a SNP list that includes the “MAF
reference set”. This is a set of 2297 SNPs that have the same minor
allele across all 1000 genomes super populations and a minor allele
frequency of 0.1-0.3.
