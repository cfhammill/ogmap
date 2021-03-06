% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/concepts.R
\name{nodes-concept}
\alias{nodes-concept}
\title{Nodes - Concept}
\description{
Nodes are the data that drive both fitting and inference from a probability field (see \link{probField-concept}).
Nodes are a convenient container used to describe a location in covariate space (guides) and contain associated
information. The additional information contained in a node can include an observation of the random variable,
a probability distribution, and extra information used in calculating the probability distribution. The node
represents a point in covariate space (it need not be unique), the coordinates of this point are given
by the set of covariates defining its location (guides for short). Value for the guides are the only
necessary peice of information for a node, the contents and utility of a node is determined by the
type of node it is. The ogmap package recognizes two informal classes of nodes, survey nodes and inference nodes.
\itemize{
\item{Survey Nodes:}{ These nodes are the observations used to create the probability field. A survey node is
a member of a collection of like nodes used for constructing the probability field mapping function. The kernel
weighted average of these nodes is used to give distributions to the inference nodes}
\item{Inference Nodes:}{ These nodes are used to learn about the random variable of interest. After construction
these nodes are given distribution functions via the probability field, so that information about likely (or unlikely)
values of the random variable at each node can be obtained.}
}
}

