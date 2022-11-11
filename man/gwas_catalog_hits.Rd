% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_snplist.R
\name{gwas_catalog_hits}
\alias{gwas_catalog_hits}
\title{GWAS top hits}
\usage{
gwas_catalog_hits(
  trait = NULL,
  efo = NULL,
  efo_id = NULL,
  map_association_to_study = TRUE,
  gwas_catalog = NULL
)
}
\arguments{
\item{trait}{the trait of interest as reported in the GWAS catalog}

\item{efo}{trait of intersest in the experimental factor ontology}

\item{efo_id}{ID for trait of interest in the experimental factor ontology}

\item{map_association_to_study}{map associations to study in GWAS catalog. This supports matching of results on PMID and study ancestry, which increases accuracy of comparisons, but is slow when there are large numbers of associations. Default = TRUE}
}
\value{
data frame
}
\description{
Extract results for top hits for the trait of interest from the NHGRI-EBI GWAS catalog
}