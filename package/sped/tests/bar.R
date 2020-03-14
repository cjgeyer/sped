
# attempt to reproduce some of the numbers from U of Washington Statistics
# Tech Report 153 (Geyer, 1988) cited in the help page for R function
# descent in this package.

library(sped)
data(alberta)

# In alberta pedigree, probability that one gene picked at random
# from individuals 1085, 1094, 1180, 1260, and 1272 are each descended
# from one gene in founder 52.
# TR says 0.0156250, 0.0156250, 0.0078125, 0.0312500, 0.0078125


vescent <- Vectorize(descent, vectorize.args = "individuals")
foo <- vescent(c(1085, 1094, 1180, 1260, 1272), alberta, c("52"=1))
all(all.equal(foo, c(0.0156250, 0.0156250, 0.0078125, 0.0312500, 0.0078125)))

