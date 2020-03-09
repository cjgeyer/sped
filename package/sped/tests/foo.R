
# attempt to reproduce some of the numbers from U of Washington Statistics
# Tech Report 153 (Geyer, 1988) cited in the help page for R function
# descent in this package.

# In alberta pedigree, probability that one gene picked at random
# from individual 1260 is descended from one gene in founder 52.
# TR says 0.03125

library(sped)
data(alberta)

foo <- descent(1260, alberta, c("52"=1))
all.equal(foo, 0.03125)

# In alberta pedigree, probability that one gene picked at random
# from each of the individuals 1085, 1094, 1180, 1260, and 1272
# are descended from one gene in founder 52.
# TR says 1.91927e-5

foo <- descent(c(1085, 1094, 1180, 1260, 1272), alberta, c("52"=1))
foo
all.equal(foo, 1.91927e-5)

