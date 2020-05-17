## ----knitr_setup, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, message=FALSE, results = "asis",
  fig.width=7, fig.height=5, fig.asp=0.5,
  collapse = TRUE,
  comment = "#>"
  # dev = 'pdf'
)

## ----setup, include = FALSE---------------------------------------------------
rm(list = ls())                      # clear objects
graphics.off()                       # close graphics windows
#
library(wavemulcor)
data(exchange)
returns <- diff(log(as.matrix(exchange)))
returns <- ts(returns, start=1970, freq=12)
N <- dim(returns)[1]

## ----label=wmc_code,echo=c(-1:-3,-6:-8)---------------------------------------
## Based on data from Figure 7.8 in Gencay, Selcuk and Whitcher (2001)
## plus one random series.

wf <- "d4"
J <- trunc(log2(N))-3

set.seed(140859)

demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
xrand.modwt <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)

xx <- list(demusd.modwt, jpyusd.modwt, xrand.modwt)
names(xx) <- c("DEM.USD","JPY.USD","rand")

Lst <- wave.multiple.correlation(xx)


## ----label=plot_wmc, echo=-1:-2, results='hide'-------------------------------
##Producing correlation plot
Lst <- wave.multiple.regression(xx)
plot_wave.multiple.correlation(Lst)

## ----label=wmr_code-----------------------------------------------------------
Lst <- wave.multiple.regression(xx)

## ----label=plot_wmr, echo=-1, results='hide'----------------------------------
##Producing regression plot
plot_wave.multiple.regression(Lst) # nsig=2)

## ----label=wmcr_code,echo=-5:-6-----------------------------------------------
wf <- "d4"
J <- trunc(log2(N))-3
lmax <- 36

set.seed(140859)

demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)

# ---------------------------

xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
names(xx) <- c("DEM.USD","JPY.USD","rand")

Lst <- wave.multiple.cross.regression(xx, lmax)


## ----label=heat_wmcc, echo=-1, results='hide'---------------------------------
##Producing correlation heat map
heatmap_wave.multiple.cross.correlation(Lst, lmax) #, xaxt="s", pdf.write=NULL)

## ----label=plot_wmcc, echo=-1, results='hide'---------------------------------
##Producing correlation plot
plot_wave.multiple.cross.correlation(Lst, lmax) #, by=2)

## ----label=plot_wmcr, echo=-1, results='hide'---------------------------------
##Producing correlation plot
plot_wave.multiple.cross.regression(Lst, lmax) #, by=2)

## ----label=wlmr_code,echo=-7:-8, results='hide'-------------------------------
wf <- "d4"
M <- 30
window <- "gauss" #uniform"

J <- trunc(log2(N))-3

set.seed(140859)

demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
xrand.modwt <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)

xx <- list(demusd.modwt, jpyusd.modwt, xrand.modwt)
names(xx) <- c("DEM.USD","JPY.USD","rand")

Lst <- wave.local.multiple.regression(xx, M, window=window) #, ymaxr=1)


## ----label=heat_wlmc, echo=-1, results='hide'---------------------------------
##Producing correlation heat map
heatmap_wave.local.multiple.correlation(Lst) #, xaxt="s", pdf.write=NULL)

## ----label=plot_wlmc, echo=-1, results='hide'---------------------------------
##Producing line plots with CI
plot_wave.local.multiple.correlation(Lst) #, xaxt="s")

## ----label=plot_wlmr, echo=-1, results='hide'---------------------------------
##Producing regression plots
plot_wave.local.multiple.regression(Lst) #, xaxt="s")

## ----label=wlmcr_code,echo=-7:-8, results='hide'------------------------------
wf <- "d4"
M <- 30
window <- "gauss" #uniform"
J <- trunc(log2(N))-3
lmax <- 5

set.seed(140859)

demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)

xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
names(xx) <- c("DEM.USD","JPY.USD","rand")

Lst <- wave.local.multiple.cross.regression(xx, M, window=window, lag.max=lmax) #, ymaxr=1)


## ----label=heat_wlmcc_lag, echo=-1, results='hide'----------------------------
##Producing cross-correlation heat map
heatmap_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE) 
             #, xaxt="s", pdf.write=NULL)

## ----label=heat_wlmcc_lev, echo=-1, results='hide'----------------------------
##Producing cross-correlation heat map
heatmap_wave.local.multiple.cross.correlation(Lst, lmax=2, lag.first=TRUE) 
             #, xaxt="s", pdf.write=NULL)

## ----label=plot_wlmcc, echo=-1, eval=FALSE------------------------------------
#  ##Producing cross-correlation plot
#  plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE) #, xaxt="s")

## ----label=plot_wlmcr, echo=-1, eval=FALSE------------------------------------
#  ##Producing cross-regression plot
#  plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2) #, xaxt="s")

