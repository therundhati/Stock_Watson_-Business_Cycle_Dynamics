uk <- window(data1.ts, start = c(1969, 1), end = c(2002, 3))
uk_g <- 4*diff(log(uk))

us <- window(data2.ts, start = c(1969, 1), end = c(2002, 3))
us_g <- 4*diff(log(us))

cn <- window(data3.ts, start = c(1969, 1), end = c(2002, 3))
cn_g <- 4*diff(log(cn))

fr <- window(data4.ts, start = c(1969, 1), end = c(2002, 3))
fr_g <- 4*diff(log(fr))

bd <- window(data5.ts, start = c(1969, 1), end = c(2002, 3))
bd_g <- 4*diff(log(bd))

it <- window(data6.ts, start = c(1969, 1), end = c(2002, 3))
it_g <- 4*diff(log(it))

jp <- window(data7.ts, start = c(1969, 1), end = c(2002, 3))
jp_g <- 4*diff(log(jp))

main <- cbind(uk_g, us_g, cn_g, fr_g, bd_g, it_g, jp_g)
# install.packages("vars")
library(vars)
var4 <- VAR(main, p=4)
restrictm <- matrix (c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
                       1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,
                       1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,
                       1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,
                       1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,
                       1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,
                       1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0
),
nrow = 7, ncol = 29, byrow = T)
result <- restrict(var4, method = "manual", resmat = restrictm)

b = matrix(cbind(1, NA, 0, 0, 0, 0, 0,
                 NA, 1, 0, 0, 0, 0, 0,
                 NA, NA, 1, 0, 0, 0, 0,
                 NA, NA, NA, 1, 0, 0, 0,
                 NA, NA, NA, NA, 1, 0, 0,
                 NA, NA, NA, NA, NA, 1, 0,
                 NA, NA, NA, NA, NA, NA, 1), nrow = 7, ncol = 7)

see <- SVAR(result, Bmat = b)
fevd(see)