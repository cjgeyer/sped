\name{descent}
\alias{descent}
\title{Multigene Descent Probabilities}
\concept{genetics}
\description{
  Calculate probabilities that genes randomly chosen from some set
  of individuals are descended from some specified set of genes.
}
\usage{
descent(individuals, pedigree, geneset, check.sex=FALSE, debug=FALSE)
}
\arguments{
  \item{individuals}{an atomic vector of names of individuals.  Must be
      of type integer or character.}
  \item{pedigree}{a matrix with three columns. Each row contains the name
      of an individual and its father and mother, in that order.  Must be
      the same type as the preceding argument.}
  \item{geneset}{an integer vector having names that are contained in the
      preceding argument.  Allowed values are 0, 1, or 2, but zero values
      may be omitted.}
  \item{check.sex}{If \code{TRUE} check that no individual appears in both
      the father column and the mother column of argument \code{pedigree}.}
  \item{debug}{If \code{TRUE} provide individuals and value for each
      recursive call.}
}
\details{
We work relative to the pedigree defined by argument \code{pedigree} in
which every individual has either two parents or none specified. Those with
none specified are called \emph{founders}.  Any ancestors of founders
are assumed to not be individuals in the pedigree, that is, we are assuming
all unknown individuals are different from all known individuals.

Thompson (1983) defines multigene descent probabilities
\eqn{g_S(B_1, \ldots, B_n)}{g_S(B_1, ..., B_n)}
to be the probability that genes at one autosomal locus randomly chosen
from each of the individuals \eqn{B_1}, \ldots, \eqn{B_n} are all descended
from genes (not necessarily the same gene) in some set \eqn{S} of genes in
individuals in the pedigree. The individuals need not be distinct (there
can be repeats). The individuals are specified by argument \code{individuals}.
The gene set \eqn{S} is specified by argument \code{geneset}.
}
\value{
If \code{debug == FALSE}, a scalar: the calculated probability.
If \code{debug == TRUE}, a list, having components

\item{value}{what would have been returned if \code{debug == FALSE}.}
\item{individuals}{the argument \code{individuals} to this call.}
\item{type}{a letter from "a" to "f" corresponding to the
    equation numbers in Theorem 1 in the design document, which is
    found in the DesignDoc directory of the installed package.}
\item{calls}{a list of calls (this function calls itself recursively).
    Each item in this list is itself a list having components \code{value},
    \code{individuals}, and (perhaps) \code{calls}.  Depending on the type
    there may be from zero to three calls.}
}
\references{
Geyer, C.~J. (1988)
Software for calculating gene survival and multigene descent
    probabilities and for pedigree manipulation and drawing.
Technical Report No. 153, Department of Statistics,
    University of Washington.
\url{https://www.stat.washington.edu/article/tech-report/software-calculating-gene-survival-and-multigene-descent-probabilities-and}.

Thompson, E. A. (1983)
Gene extinction and allelic origins in complex genealogies (with discussion).
\emph{Proceedings of the Royal Society of London. Series B,
    Biological Sciences}, \bold{219}, 241--251.
\doi{doi:10.1098/rspb.1983.0072}.

Thompson, E. A. (1986)
Ancestry of alleles and extinction of genes in populations with defined
    pedigrees.
\emph{Zoo Biology}, \bold{5}, 161--170.
\doi{doi:10.1002/zoo.1430050210}.
}
\keyword{misc}
