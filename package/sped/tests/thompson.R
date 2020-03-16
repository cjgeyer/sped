
library(sped)
data(thompson)

# Attempt to reproduce results from Table 1 in Thompson (1986)
# cited in help pages and design doc

goo <- gammas(c("U", "V", "Q", "R", "W"), thompson)
# goo

# Table 1 only has numbers for founders B, C, F
goo <- goo[c("B", "C", "F"), ]

# transpose to match table 1
goo <- t(goo)
# goo

mygoo <- rbind(c(6, 6, 4) / 16,
               c(20, 20, 8) / 64,
               c(8, 0, 8) / 32,
               c(8, 4, 16) / 32,
               c(72, 56, 80) / 256)

all(dim(goo) == dim(mygoo))

dimnames(mygoo) <- dimnames(goo)

all.equal(goo, mygoo)

# now betas

boo <- betas(c("U", "V", "Q", "R", "W"), thompson)
# boo

# Table 1 only has numbers for founders B, C, F
boo <- boo[c("B", "C", "F"), ]

# transpose to match table 1
boo <- t(boo)
# boo

myboo <- rbind(c(2, 2, 0) / 16,
               c(6, 6, 0) / 64,
               c(2, 0, 0) / 32,
               c(2, 0, 8) / 32,
               c(20, 10, 16) / 256)

all(dim(boo) == dim(myboo))

dimnames(myboo) <- dimnames(boo)

all.equal(boo, myboo)

# now alphas

aoo <- alphas(c("U", "V", "Q", "R", "W"), thompson)
# aoo

# Table 1 only has numbers for founders B, C, F
aoo <- aoo[c("B", "C", "F"), ]

# transpose to match table 1
aoo <- t(aoo)
# aoo

myaoo <- rbind(c(1, 1, 0) / 16,
               c(3, 3, 0) / 64,
               c(1, 0, 0) / 32,
               c(1, 0, 4) / 32,
               c(11, 6, 8) / 256)

all(dim(aoo) == dim(myaoo))

dimnames(myaoo) <- dimnames(aoo)

all.equal(aoo, myaoo)

# now inbreeding

ioo <- inbreeding(c("U", "V", "Q", "R", "W"), thompson)
# ioo

myioo <- c(1/8, 3/32, 1/16, 5/32, 25/256)

length(ioo) == length(myioo)

names(myioo) <- names(ioo)

all.equal(ioo, myioo)

