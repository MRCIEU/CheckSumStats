% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_snps.R
\name{transform_betas}
\alias{transform_betas}
\title{Transform betas}
\usage{
transform_betas(dat = NULL, effect = "lnor", effect.se = "se")
}
\arguments{
\item{dat}{the target dataset
rsids}

\item{effect}{the column containing the beta. We wish to transform this to a log odds ratio scale}

\item{effect.se}{standard error for the beta}
}
\value{
data frame
}
\description{
Transform betas from a linear model to a log odds ratio scale. Assumes betas have been derived from a linear model of case-control status regressed on SNP genotype (additively coded).
}
