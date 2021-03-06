% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/probField.R
\name{getSurvey}
\alias{getSurvey}
\title{Get Survey Nodes}
\usage{
getSurvey(surveySource)
}
\arguments{
\item{surveySource}{The source of survey nodes}
}
\value{
A \code{probFieldNodeList} or NULL if the argument is neither a \code{probField} or
\code{probFieldNodeList}
}
\description{
Return a \code{probFieldNodeList} either identical to its argument
if a \code{probFieldNodeList} is passed, or the contained list of survey nodes is a
\code{probField} is passed
}

