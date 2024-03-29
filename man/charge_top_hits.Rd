% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/charge_top_hits.R
\docType{data}
\name{charge_top_hits}
\alias{charge_top_hits}
\title{GWAS top hits for arachidonic acid in the CHARGE consortium}
\format{
A charactor vector of length 210:
}
\source{
\url{http://www.chargeconsortium.com/main/results}
}
\usage{
charge_top_hits
}
\description{
The dataset contains rsids for single nucleotide polymorphisms extracted from a genome-wide association study of arachidonic acid in the CHARGE consortium. The list was generated by 1) extracting all SNPs with P values <5e-8 (1063 SNPs in total); and then 2) performing LD clumping on the 1063 extracted SNPs (clump_r2 = 0.01, clump_kb=10000 ) using European participants from UK Biobank as a reference dataset. Clumping was performed using ieugwasr::ld_clump. No post-GWAS filtering on allele frequency, imputation info score or number of studies was performed on the GWAS summary statistics prior to the extraction of the SNPs.
}
\keyword{datasets}
