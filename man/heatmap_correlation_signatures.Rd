% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visu.R
\name{heatmap_correlation_signatures}
\alias{heatmap_correlation_signatures}
\title{heatmap_correlation_signatures}
\usage{
heatmap_correlation_signatures(
  X,
  model.lasso,
  model.cluster,
  clust.min = 10,
  plot.title = ""
)
}
\arguments{
\item{model.lasso}{lasso model (returned by "extract_best_model" function)}

\item{model.cluster}{lasso model (returned by "extract_best_model" function)}

\item{clust.min}{minimum size of clusters to consider in the representation (default : 10)}

\item{plot.title}{title for the plot (default = "" - empty)}
}
\description{
This function generates a heatmap figure comparing the signatures of a lasso and a cluster-lasso model (returned by the "extract_best_model" function)
}
\examples{
heatmap_correlation_signatures()
}
