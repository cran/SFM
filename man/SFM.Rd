% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SFM.R
\name{SFM}
\alias{SFM}
\title{The SFM function is to generate Skew Factor Models data.}
\usage{
SFM(n, p, m, xi, omega, alpha, distribution_type)
}
\arguments{
\item{n}{Sample size.}

\item{p}{Sample dimensionality.}

\item{m}{Number of factors.}

\item{xi}{A numerical parameter used exclusively in the "Skew-t" distribution, representing the distribution's xi parameter.}

\item{omega}{A numerical parameter representing the omega parameter of the distribution, which affects the degree of skewness in the distribution.}

\item{alpha}{A numerical parameter representing the alpha parameter of the distribution, which influences the shape of the distribution.}

\item{distribution_type}{The type of distribution.}
}
\value{
A list containing:
\item{data}{A matrix of generated data.}
\item{A}{A matrix representing the factor loadings.}
\item{D}{A diagonal matrix representing the unique variances.}
}
\description{
The function supports various distribution types for generating the data,
including: Skew-Normal Distribution, Skew-Cauchy Distribution, Skew-t Distribution.
}
\examples{
library(MASS)
library(SOPC)
library(sn)
library(matrixcalc)
library(psych)
n <- 100
p <- 10
m <- 5
xi <- 5
omega <- 2
alpha <- 5
distribution_type <- "Skew-Normal Distribution"
X <- SFM(n, p, m, xi, omega, alpha, distribution_type)

}
