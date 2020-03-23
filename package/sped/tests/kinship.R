
library(sped)
library(pooh)
data(thompson)

individuals <- sort(unique(thompson))
foo <- kinship(individuals, thompson)
individuals.too <- rownames(foo)
setequal(individuals, individuals.too)

founders <- setdiff(individuals, thompson[ , "ind"])

bar <- foo
for (i in 1:nrow(bar))
    for (j in 1:ncol(bar)) {
        baz <- double(length(founders))
        for (k in seq(along = founders)) {
            g <- 1
            names(g) <- founders[k]
            baz[k] <- descent(individuals.too[c(i, j)], thompson, g)
        }
        bar[i, j] <- bar[j, i] <- 2 * sum(baz)
    }

all.equal(foo, bar)

qux <- kinship(c("U", "V", "Q", "R", "W"), thompson)
idx <- match(rownames(qux), rownames(foo))
all.equal(foo[idx, idx], qux)

