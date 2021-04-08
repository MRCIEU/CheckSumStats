#' A example dataset of genetic summary data for arachidonic acid
#'
#' The dataset contains summary association statistics for 436 SNPs, generated in linear regression models, from a genome-wide association study of arachidonic acid conducted by the CHARGE consortium. No post-GWAS filtering on allele frequency, imputation info score or number of studies has been performed. The selected SNPs correspond to three groups: 1) A MAF 1KG reference set, 2) GWAS catalog top hits for arachidonic acid and 3) GWAS top hits for arachidonic acid in the CHARGE study  
#'
#' @format A data frame with 436 rows and 9 variables:
#' \describe{
#'   \item{snp}{SNP rsid}
#'   \item{effect_allele}{effect allele}
#'   \item{other_allele}{non-effect allele}
#'   \item{effect_allele_freq}{effect allele frequency}
#'   \item{beta}{change in arachidonic acid per copy of the effect allele}
#'   \item{se}{standard error for beta}
#'   \item{p }{p value statistic describing the association between the SNP and arachidonic acid}
#'   \item{n }{number of study participants}
#'   \item{path_to_target_file}{name of file used to generate example dataset}
#' }
#' @source \url{http://www.chargeconsortium.com/main/results}
"ara_test_dat"



  