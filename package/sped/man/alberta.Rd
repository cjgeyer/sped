\name{alberta}
\docType{data}
\alias{alberta}
\title{Pedigree Data on Equus przewalskii in Alberta}
\description{
  Pedigree of Asian wild horse (also called Mongolian wild horse)
  \emph{Equus przewalskii} individuals living in Alberta, Canada in 1988.
}
\usage{alberta}
\format{
A matrix with 69 rows giving (individual, father, mother) triplets.
The column labels are \code{ind}, \code{pa}, and \code{ma}.
The names (which are numbers) are studbook numbers.
There are 82 individuals in this pedigree in all.
There are \eqn{82 - 69 = 13} founders.
}
\examples{
library(sped)
data(alberta)
}
\keyword{datasets}
