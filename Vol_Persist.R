uk <- window(data1.ts, start = c(1983, 4), end = c(2002,3))
uk_g <- 4*diff(log(uk))
uk_v <- c(uk_g)
library(strucchange)
breakpoints(uk_v ~ lag(uk_v, 1) + lag(uk_v, 2) + lag(uk_v, 3) + lag(uk_v, 4), breaks = 1)
ICSS(uk_v)
ar(uk_g, aic = FALSE, order.max = 4)

us <- window(data2.ts, start = c(1983, 4), end = c(2002,3))
us_g <- 4*diff(log(us))
us_v <- c(us_g)
breakpoints(us_v ~ lag(us_v, 1) + lag(us_v, 2) + lag(us_v, 3) + lag(us_v, 4), breaks = 1)
ICSS(us_v)
ar(us_g, aic = FALSE, order.max = 4)

cn <- window(data3.ts, start = c(1959, 4), end = c(1983,3))
cn_g <- 4*diff(log(cn))
cn_v <- c(cn_g)
breakpoints(cn_v ~ lag(cn_v, 1) + lag(cn_v, 2) + lag(cn_v, 3) + lag(cn_v, 4), breaks = 1)
# install.packages("ICSS")
library(ICSS)
ICSS(cn_v)

fr <- window(data4.ts, start = c(1983, 4), end = c(2002,3))
fr_g <- 4*diff(log(fr))
fr_v <- c(fr_g)
breakpoints(fr_v ~ lag(fr_v, 1) + lag(fr_v, 2) + lag(fr_v, 3) + lag(fr_v, 4), breaks = 1)
ICSS(fr_v)
ar(fr_g, aic = FALSE, order.max = 4)

bd <- window(data5.ts, start = c(1983, 4), end = c(2002, 3))
bd_g <- 4*diff(log(bd))
bd_v <- c(bd_g)
breakpoints(bd_v ~ lag(bd_v, 1) + lag(bd_v, 2) + lag(bd_v, 3) + lag(bd_v, 4), breaks = 1)
ICSS(bd_v)
ar(bd_g, aic = FALSE, order.max = 4)

it <- window(data6.ts, start = c(1983, 4), end = c(2002,3))
it_g <- 4*diff(log(it))
it_v <- c(it_g)
breakpoints(it_v ~ lag(it_v, 1) + lag(it_v, 2) + lag(it_v, 3) + lag(it_v, 4), breaks = 1)
ICSS(it_v)
ar(it_g, aic = FALSE, order.max = 4)

jp <- window(data7.ts, start = c(1983, 4), end = c(2002, 3))
jp_g <- 4*diff(log(jp))
jp_v <- c(jp_g)
breakpoints(jp_v ~ lag(jp_v, 1) + lag(jp_v, 2) + lag(jp_v, 3) + lag(jp_v, 4), breaks = 1)
ICSS(jp_v)
ar(jp_g, aic = FALSE, order.max = 4)