\encoding{latin1}
\name{wavemulcor-package}
\alias{wavemulcor-package}
\alias{wavemulcor}
\docType{package}
\title{ Wavelet routine for multiple correlation
}
\description{
 Produces an estimate of the multiscale multiple correlation
 (as defined below) along with approximate confidence intervals.
}
\details{
 \tabular{ll}{
 Package: \tab wavemulcor\cr
 Type: \tab Package\cr
 Version: \tab 1.2\cr
 Date: \tab 2011-02-13\cr
 License: \tab GPL (>= 2)\cr
 LazyLoad: \tab yes\cr
 }
 The \emph{wavemulcor} package contains two routines, \emph{wave.multiple.correlation} and\cr \emph{wave.multiple.cross.correlation},
 that calculate single sets of, respectively, wavelet multiple correlations and wavelet multiple cross-correlations out of \eqn{n} variables. They can later be plotted in single graphs, as an alternative to trying to make sense out of \eqn{n(n-1)/2} sets of wavelet correlations or \eqn{n(n-1)/2 \times J}{n(n-1)/2 x J} sets of wavelet cross-correlations. The code is based on the calculation, at each wavelet scale, of the square root of the coefficient of determination in a linear combination of variables for which such coefficient of determination is a maximum.
% The code provided here is based on the wave.correlation routine in Brandon Whitcher's waveslim R package Version: 1.6.4, which in turn is based on wavelet methodology developed in Percival and Walden (2000); Gencay, Selcuk and Whitcher (2001) and others.
}
\note{
 Dependencies: waveslim
}
\author{
 Javier Fernández-Macho (UPV/EHU)\cr
 Maintainer: Javier Fernández-Macho  <javier.fernandezmacho@ehu.es>
}
\references{
 Fernández-Macho, Javier (2011) The wavelet multiple correlation, (mimeo).\cr
 Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An Introduction to Wavelets and Other Filtering
 Methods in Finance and Economics}, Academic Press.\cr
 Whitcher, B. (2010) \emph{`waveslim' R Package.}
}
