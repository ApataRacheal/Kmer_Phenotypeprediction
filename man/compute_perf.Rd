% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{compute_perf}
\alias{compute_perf}
\title{compute_perf}
\usage{
compute_perf(preds, probs, y)
}
\arguments{
\item{preds}{[n_sample x n_models] matrix containing binary (-1/+1) predictions}

\item{probs}{[n_sample x n_models] matrix containing posterior probability p(y=1|x) (more generally : score associated to positive class)}

\item{y}{reference labels}
}
\description{
This function computes several performance indicators for binary classification
}
\examples{
compute_perf()
}
