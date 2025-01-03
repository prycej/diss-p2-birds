library(wiqid)
data(dippers)
str(dippers)
CH <- as.matrix(dippers[, 1:7])
dim(CH)
head(CH)

library(IPMbook)  # for helper functions
f <- getFirst(CH)
table(f)
CH <- CH[f < 7, ]
nrow(CH)  # 255
f <- getFirst(CH)
table(f)

jdata <- list(y = CH, f = f, nind=nrow(CH), n.occasions=ncol(CH),
              z = zKnown(CH))
str(jdata)

wanted = c("phi", "p")

( out_phi._p. <- jags(jdata, NULL, wanted,
                      "models/CJS_phi(.)p(.).jags", DIC=FALSE,
                      n.chains=3, n.adapt=1000, n.iter=10000, parallel = TRUE) )


plot(out_phi._p.)
