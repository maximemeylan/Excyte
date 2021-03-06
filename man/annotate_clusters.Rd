% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_excyte.R
\name{annotate_clusters}
\alias{annotate_clusters}
\title{Provide annotations for each cluster, based on intensity distribution of all events}
\usage{
annotate_clusters(
  phenograph_obj,
  channels = "all",
  cluster_to_use = "all",
  threshold = c("median", "tertile", "quartile")[1],
  positivity_threshold = 0.5,
  channel_names = c("channel_only", "marker_only", "both")[3]
)
}
\arguments{
\item{phenograph_obj}{list containing result of phenograph clustering and processed fcs}

\item{channels}{vector of channels to use, default uses all channels}

\item{cluster_to_use}{vector of cluster to use, default uses all clusters}

\item{threshold}{character defining if threshold should be caracterized as the median, tertile or quartile}

\item{positivity_threshold}{numeric value between 0 and 1 defining the percentage of cells needed to call positivity to a threshold}

\item{channel_names}{character, edit channels names accordingly}
}
\description{
Provide annotations for each cluster, based on intensity distribution of all events
}
