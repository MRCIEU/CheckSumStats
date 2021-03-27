#' A example dataset of genetic summary data
#'
#' The dataset contains summary association statistics for 98 SNPs, generated in logistic regression models, from a genome-wide association study of glioma conducted by the GliomaScan consortium. 
#'
#' @format A data frame with 98 rows and 20 variables:
#' \describe{
#'   \item{Locus}{SNP rsid}
#'   \item{Allele1}{non-effect allele}
#'   \item{Allele2}{effect allele}
#'   \item{MAF}{SNP minor allele frequency in controls|cases}
#'   \item{Geno_Counts}{genotype counts in controls/cases}
#'   \item{Subjects}{Number of participants in study }
#'   \item{p }{p value statistic describing the association between the SNP and glioma}
#'   \item{OR }{odds ratio for glioma}
#'   \item{  OR_95._CI_l }{lower 95% confidence interval}
#'   \item{OR_95._CI_u }{upper 95% confidence interval}
#'   \item{CHROMOSOME}{chromosome number}
#'   \item{LOCATION}{genomic coordinates in base pairs}
#'   \item{controls}{number of controls}
#'   \item{cases}{number of cases}
#'   \item{eaf.controls}{effect allele frequency in controls}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/22886559/}
"glioma_test_dat"



  