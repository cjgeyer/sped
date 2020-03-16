\name{thompson}
\docType{data}
\alias{thompson}
\title{Example Pedigree Data}
\description{
  Pedigree of Example data from Thompson (1986).
}
\usage{data(thompson)}
\format{
A matrix with 15 rows giving (individual, father, mother) triplets.
The column labels are \code{ind}, \code{pa}, and \code{ma}.
The names are letters of the alphabet.
There are 21 individuals in this pedigree in all.
There are \eqn{21 - 15 = 6} founders.
}
\references{
Thompson, E. A. (1986)
Ancestry of alleles and extinction of genes in populations with defined
    pedigrees.
\emph{Zoo Biology}, \bold{5}, 161--170.
\doi{doi:10.1002/zoo.1430050210}.
}
\examples{
library(sped)
data(thompson)
head(thompson)
}
\keyword{datasets}
