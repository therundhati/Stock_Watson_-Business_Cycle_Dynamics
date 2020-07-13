data1 <- as.data.frame(read.csv("~/Desktop/Course/TSF/Paper/Mine/cngdppc.csv", sep = ","))
library(pastecs)
# stat.desc(data1[[1]])
sub1 <- data1[171:211,]     # to retrieve standard deviation of 1990-2002
stat.desc(sub1)


data1 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/ukgdppc.csv", sep = ","))
library(pastecs)
stat.desc(data1[[1]])
data1.ts <- ts(data1, start=c(1950, 1), end=c(2002, 4), frequency=4)
sub1 <- window(data1.ts, start = c(1954, 4), end = c(2002,3))
# plot(sub1)

install.packages("mFilter")
library(mFilter)
library(dplyr)

hp_gdp <- hpfilter(log(sub1), freq = 1600)
# sub1.trans <- mutate(hp = hp_gdp$cycle, lin_cycle = log(sub1) - lin_trend)
plot(hp_gdp$cycle, xlab='Year', ylab='UK')

data2 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/usgdppc.csv", sep = ","))
data2.ts <- ts(data2, start=c(1950, 1), end=c(2002, 3), frequency=4)
hp_gdp1 <- hpfilter(log(data2.ts), freq = 1600)
plot(hp_gdp1$cycle, xlab='Year', ylab='US')

data3 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/cngdppc.csv", sep = ","))
data3.ts <- ts(data3, start=c(1950, 1), end=c(2002, 4), frequency=4)
sub3 <- window(data3.ts, start = c(1959, 4), end = c(2002,3))
hp_gdp3 <- hpfilter(log(sub3), freq = 1600)
plot(hp_gdp3$cycle, xlab='Year', ylab='Canada')

data4 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/frgdppc.csv", sep = ","))
data4.ts <- ts(data4, start=c(1950, 1), end=c(2002, 4), frequency=4)
sub4 <- window(data4.ts, start = c(1959, 4), end = c(2002,3))
hp_gdp4 <- hpfilter(log(sub4), freq = 1600)
plot(hp_gdp4$cycle, xlab='Year', ylab='France')

data5 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/bdgdppc.csv", sep = ","))
data5.ts <- ts(data5, start=c(1950, 1), end=c(2002, 4), frequency=4)
sub5 <- window(data5.ts, start = c(1959, 4), end = c(2002,3))
hp_gdp5 <- hpfilter(log(sub5), freq = 1600)
plot(hp_gdp5$cycle, xlab='Year', ylab='Germany')

data6 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/itgdppc.csv", sep = ","))
data6.ts <- ts(data6, start=c(1950, 1), end=c(2002, 4), frequency=4)
sub6 <- window(data6.ts, start = c(1959, 4), end = c(2002,3))
hp_gdp6 <- hpfilter(log(sub6), freq = 1600)
plot(hp_gdp6$cycle, xlab='Year', ylab='Italy')

data7 <- as.data.frame(read.csv( "~/Desktop/Course/TSF/Paper/Mine/jpgdppc.csv", sep = ","))
data7.ts <- ts(data7, start=c(1950, 1), end=c(2002, 4), frequency=4)
sub7 <- window(data7.ts, start = c(1959, 4), end = c(2002,3))
hp_gdp7 <- hpfilter(log(sub7), freq = 1600)
plot(hp_gdp7$cycle, xlab='Year', ylab='Japan')