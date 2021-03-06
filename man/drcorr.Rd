% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dr.corr.R
\name{drcorr}
\alias{drcorr}
\title{Correlation function
this function will calculate all correlations between data matrix 1 and 2 by row}
\usage{
drcorr(datafiles, job_dir, outformat = "json", method = "spearman",
  alternative = "two.sided", min_samples = 5, pvalcorrection = "fdr")
}
\arguments{
\item{datafiles}{contain file1, data matrix used for analysis, each column is a data vector
and file2, data matrix used for analysis, each column is a data vector
if file2 not given, calculate correlations within file1}

\item{method, }{pearson, spearman, kendall (default spearman)}

\item{alternative, }{two.sided, less, greater (default two.sided)}

\item{min_samples, }{minimum samples for correlation calculation}

\item{pvalcorrection, }{fdr, bonferroni}
}
\value{
results
}
\description{
Correlation function
this function will calculate all correlations between data matrix 1 and 2 by row
}
