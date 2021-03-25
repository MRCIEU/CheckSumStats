#' A dataset of reference allele frequencies from 1000 genomes superpopulations
#'
#' The dataset contains minor allele frequency for 2297 SNPs that have minor allele frequency 0.1-0.3 across each superpopulation in the 1000 genomes project. 
#'
#' @format A data frame with 13782 rows and 8 variables:
#' \describe{
#'   \item{CHR}{chromosome number}
#'   \item{SNP}{SNP rsid}
#'   \item{minor_allele}{SNP minor allele}
#'   \item{major_allele}{SNP major allele}
#'   \item{MAF}{SNP minor allele frequency}
#'   \item{NCHROBS}{number of observed chromosomes}
#'   \item{population}{1000 genomes superpopulation: AFR=African; ALL=all individuals;  AMR = Ad Mixed American; EAS=East Asian; EUR=European;  SAS=South Asian}
#' }
#' @source \url{https://www.internationalgenome.org/home}
"refdat_1000G_superpops"