\name{wavemulcor-package}
\alias{wavemulcor-package}
\alias{wavemulcor}
\docType{package}
\title{ Wavelet package for global and local multiple correlation
}
\description{
 Produces estimates of multiscale global or local multiple correlation (as defined below) along with approximate confidence intervals.
}
\details{
 \tabular{ll}{
 Package: \tab wavemulcor\cr
 Type: \tab Package\cr
 Version: \tab 2.1.0\cr
 Date: \tab 2017-07-01\cr
 License: \tab GPL (>= 2)\cr
 LazyLoad: \tab yes\cr
 }
 The wavemulcor package contains three routines, wave.multiple.correlation, wave.multiple.cross.correlation and wave.local.multiple.correlation,
that calculate single sets of, respectively, global wavelet multiple correlations, global wavelet multiple cross-correlations and time-localized wavelet multiple correlations out of n variables.
They can later be plotted in single graphs, as an alternative to trying to make sense out of \eqn{n(n-1)/2} sets of global wavelet correlations or \eqn{n(n-1)/2 . J} sets of global wavelet cross-correlations or \eqn{n(n-1)/2 . [J x T]} sets of local wavelet correlations respectively.
The code is based on the calculation, at each wavelet scale, of the square root of the coefficient of determination in that linear combination of either global or locally weighted wavelet coefficients for which such coefficient of determination is a maximum.
}
\note{
 Dependencies: waveslim
}
\author{
 Javier Fernández-Macho (UPV/EHU)\cr
 Maintainer: Javier Fernández-Macho  <javier.fernandezmacho@ehu.eus>
}
\references{
 Fernández-Macho, J., 2012. Wavelet multiple correlation and cross-correlation: A multiscale analysis of Eurozone stock   markets. Physica A: Statistical Mechanics and its Applications 391, 1097–1104. doi:10.1016/j.physa.2011.11.002\cr
 Fernández-Macho, J., 2017. Time-localized wavelet multiple regression and correlation, (mimeo).\cr
 Gençay, R., F. Selçuk and B. Whitcher, 2001. \emph{An Introduction to Wavelets and Other Filtering
 Methods in Finance and Economics}, Academic Press.\cr
 Whitcher, B., 2010. \emph{`waveslim' R Package.}
}
