\encoding{latin1}
\name{wave.multiple.cross.correlation}
\alias{wave.multiple.cross.correlation}
\alias{wmcc}
\title{
 Wavelet routine for multiple cross-correlation
}
\description{
 Produces an estimate of the multiscale multiple cross-correlation
 (as defined below).
}
\usage{
 wave.multiple.cross.correlation(xx, lag.max = NULL, ymaxr = NULL)
}
\arguments{
  \item{xx}{
 A list of \eqn{n} (multiscaled) time series, usually the outcomes of dwt or modwt, \emph{i.e.}
 xx <- list(v1.modwt.bw, v2.modwt.bw, v3.modwt.bw)}
  \item{lag.max}{
 maximum lag. If not set, it defaults to half the square root of the length of the
 original series.
}
  \item{ymaxr}{
 index number of the variable whose correlation is calculated against a linear combination
 of the rest, otherwise at each wavelet level wmc chooses
 the one maximizing the multiple correlation.
}
}
\details{
 The routine calculates one single set of wavelet multiple cross-correlations out of \eqn{n} variables
 that can be plotted as one single set of graphs (one per wavelet level), as an
 alternative to trying to make sense out of \eqn{n(n-1)/2 \times J}{n(n-1)/2 x J} sets of wavelet cross-correlations. The code is
 based on the calculation, at each wavelet scale, of the square root of the coefficient of
 determination in a linear combination of variables that includes a lagged variable for which such coefficient of determination
 is a maximum. 
 %The code provided here is based on the wave.correlation routine in Brandon
% Whitcher's waveslim R package Version: 1.6.4, which in turn is based on wavelet methodology
% developed in Percival and Walden (2000); Gencay, Selcuk and Whitcher (2001) and others.
}
\value{
 List of two elements:\cr
 \emph{xy.mulcor:} matrix with as many rows as levels in the wavelet transform object.
 The columns provide the point estimates for the wavelet multiple cross-correlations at different lags.\cr
 \emph{YmaxR:} numeric vector giving, at each wavelet level, the index number of the variable
 whose correlation is calculated against a linear combination of the rest. By default,
 \emph{wmcc} chooses at each wavelet level the variable maximizing the multiple correlation.
}
\references{
 Fernández-Macho, Javier (2011) The wavelet multiple correlation, (mimeo).
}
\author{
 Javier Fernández-Macho, Dpt. of Econometrics and Statistics,
 & Instituto de Economía Pública,
 University of the Basque Country,
 Agirre Lehendakari etorb. 83, E48015 BILBAO, Spain. (email: javier.fernandezmacho@ehu.es).
}
\note{
 Needs \emph{waveslim} package to calculate \emph{dwt} or \emph{modwt} coefficients as inputs to the routine (also for data in the example).
}
\examples{
## Based on data from Figure 7.9 in Gencay, Selcuk and Whitcher (2001)
## plus one random series.

library(wavemulcor)
data(exchange)
returns <- diff(log(exchange))
returns <- ts(returns, start=1970, freq=12)
wf <- "d4"
J <- 6
lmax <- 36
n <- dim(returns)[1]

demusd.modwt <- modwt(returns[,"DEM.USD"], wf, J)
demusd.modwt.bw <- brick.wall(demusd.modwt, wf)
jpyusd.modwt <- modwt(returns[,"JPY.USD"], wf, J)
jpyusd.modwt.bw <- brick.wall(jpyusd.modwt, wf)
rand.modwt <- modwt(rnorm(length(returns[,"DEM.USD"])), wf, J)
rand.modwt.bw <- brick.wall(rand.modwt, wf)

##xx <- list(demusd.modwt.bw, jpyusd.modwt.bw)
xx <- list(demusd.modwt.bw, jpyusd.modwt.bw, rand.modwt.bw)

Lst <- wave.multiple.cross.correlation(xx, lmax)
returns.cross.cor <- as.matrix(Lst$xy.mulcor[1:J,])
YmaxR <- Lst$YmaxR

exchange.names <- c("DEM.USD", "JPY.USD", "RAND")
rownames(returns.cross.cor)<-rownames(returns.cross.cor, 
  do.NULL = FALSE, prefix = "Level ")
lags <- length(-lmax:lmax)

lower.ci <- tanh(atanh(returns.cross.cor) - qnorm(0.975) /
sqrt(matrix(trunc(n/2^(1:J)), nrow=J, ncol=lags)- 3))
upper.ci <- tanh(atanh(returns.cross.cor) + qnorm(0.975) /
sqrt(matrix(trunc(n/2^(1:J)), nrow=J, ncol=lags)- 3))

par(mfrow=c(3,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,0,0))
for(i in J:1) {
matplot((1:(2*lmax+1)),returns.cross.cor[i,], type="l", lty=1, ylim=c(-1,1), 
  xaxt="n", xlab="", ylab="", main=rownames(returns.cross.cor)[[i]][1])
if(i<3) {axis(side=1, at=seq(1, 2*lmax+1, by=12), 
  labels=seq(-lmax, lmax, by=12))}
#axis(side=2, at=c(-.2, 0, .5, 1))
lines(lower.ci[i,], lty=1, col=2) ##Add Connected Line Segments to a Plot
lines(upper.ci[i,], lty=1, col=2)
abline(h=0,v=lmax+1)              ##Add Straight horiz and vert Lines to a Plot
text(1,1, labels=exchange.names[YmaxR[i]], adj=0.25, cex=.8)
}
par(las=0)
mtext('Lag (months)', side=1, outer=TRUE, adj=0.5)
mtext('Wavelet Multiple Cross-Correlation', side=2, outer=TRUE, adj=0.5)
}
\keyword{ statistics }
\keyword{ multivariate }
\keyword{ cross-correlation }
\keyword{ wavelet }
