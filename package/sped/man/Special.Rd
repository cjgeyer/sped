\name{Special}
\alias{Special}
\alias{alphas}
\alias{betas}
\alias{gammas}
\alias{inbreeding}
\title{Special Multigene Descent Probabilities}
\concept{genetics}
\description{
  Calculate particular multigene descent probabilities of particular interest.
}
\usage{
alphas(individuals, pedigree)
betas(individuals, pedigree)
gammas(individuals, pedigree)
inbreeding(individuals, pedigree)
}
\arguments{
  \item{individuals}{an atomic vector of names of individuals.  Must be
      of type integer or character.}
  \item{pedigree}{a matrix with three columns. Each row contains the name
      of an individual and its father and mother, in that order.  Must be
      the same type as the preceding argument.}
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
A matrix whose rows are probabilities relating to particular founders
and whose columns are probabilities relating to particular individuals.
The row and column labels say which founders and which individuals.

R function \code{gammas} gives the probability that a gene chosen at random
from the individual comes from either gene of the founder.

R function \code{betas} gives the probability that both genes
from the individual comes from the genes (not necessarily the same gene)
of the founder.

R function \code{alphas} gives inbreeding of the individual relative
to the founder: the probability that both genes
from the individual comes from the same gene (either gene) of the founder.

R function \code{inbreeding} gives inbreeding coefficients of the individuals.
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
\seealso{\code{\link{descent}}}
\examples{
# In alberta pedigree, probability that one gene picked at random
# from individual 1260 is descended from one gene in founder 52.
library(sped)

data(alberta)
gammas(c(1260, 1272), alberta)

data(thompson)
betas(c("U", "V", "Q", "R", "W"), thompson)
alphas(c("U", "V", "Q", "R", "W"), thompson)
inbreeding(c("U", "V", "Q", "R", "W"), thompson)
}
\keyword{misc}
