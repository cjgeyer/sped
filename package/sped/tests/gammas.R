# attempt to reproduce some of the numbers from U of Washington Statistics
# Tech Report 153 (Geyer, 1988) cited in the help page for R function
# descent in this package.

# In alberta pedigree, probability that one gene picked at random
# from each of the individuals 1260 and 1272 are descended from either
# gene in each founder

library(sped)
data(alberta)

foo <- gammas(c(1260, 1272), alberta)

bar <- matrix(c(0.1250000,  0.03125000,
                0.06250000, 0.01562500, 
                0.01562500, 0.0,
                0.00781250, 0.0,
                0.02734375, 0.05078125,
                0.02734375, 0.05078125,
                0.09765625, 0.1210937,
                0.08203125, 0.1054687,
                0.06250000, 0.01562500,
                0.1132812,  0.09765625,
                0.3398437,  0.2929687,
                0.03125000, 0.2187500,
                0.00781250, 0.0),
                ncol = 2, byrow = TRUE)

all(dim(foo) == dim(bar))

dimnames(bar) <- dimnames(foo)

all.equal(foo, bar, tolerance = 1e-6)

