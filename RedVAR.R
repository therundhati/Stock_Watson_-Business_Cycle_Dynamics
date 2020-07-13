data1 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/ukgdppc.csv", sep = ","))
library(pastecs)
stat.desc(data1[[1]])
data1.ts <- ts(data1, start=c(1950, 1), end=c(2002, 4), frequency=4)
uk <- window(data1.ts, start = c(1959, 4), end = c(1983,3))
uk_g <- 4*diff(log(uk))

data2 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/usgdppc.csv", sep = ","))
data2.ts <- ts(data2, start=c(1950, 1), end=c(2002, 3), frequency=4)
us <- window(data2.ts, start = c(1959, 4), end = c(1983,3))
us_g <- 4*diff(log(us))

data3 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/cngdppc.csv", sep = ","))
data3.ts <- ts(data3, start=c(1950, 1), end=c(2002, 4), frequency=4)
cn <- window(data3.ts, start = c(1959, 4), end = c(1983,3))
cn_g <- 4*diff(log(cn))

data4 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/frgdppc.csv", sep = ","))
data4.ts <- ts(data4, start=c(1950, 1), end=c(2002, 4), frequency=4)
fr <- window(data4.ts, start = c(1959, 4), end = c(1983,3))
fr_g <- 4*diff(log(fr))

data5 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/bdgdppc.csv", sep = ","))
data5.ts <- ts(data5, start=c(1950, 1), end=c(2002, 4), frequency=4)
bd <- window(data5.ts, start = c(1959, 4), end = c(1983,3))
bd_g <- 4*diff(log(bd))

data6 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/itgdppc.csv", sep = ","))
data6.ts <- ts(data6, start=c(1959, 1), end=c(2002, 4), frequency=4)
it <- window(data6.ts, start = c(1959, 4), end = c(1983,3))
it_g <- 4*diff(log(it))

data7 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/jpgdppc.csv", sep = ","))
data7.ts <- ts(data7, start=c(1950, 1), end=c(2002, 4), frequency=4)
jp <- window(data7.ts, start = c(1959, 4), end = c(1983,3))
jp_g <- 4*diff(log(jp))

main <- cbind(uk_g, us_g, cn_g, fr_g, bd_g, it_g, jp_g)
# install.packages("vars")
library(vars)
var4 <- VAR(main, p=4)
restrictm <- matrix (c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,
                       1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,
                       1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,
                       1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,
                       1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,
                       1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,
                       1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1
                      ),
                    nrow = 7, ncol = 29, byrow = T)
result <- restrict(var4, method = "manual", resmat = restrictm)

resultpred <- predict(result)
par(mar = rep(2, 4))
fanchart(resultpred)