#----------------------------------------------------------
local.multiple.correlation <-
  function(xx, M, window="gauss", p = .975, ymaxr=NULL) {
    window <- substr(tolower(window),1,4)
    if (window=="unif"){
      #uniform window:
      weights <- function(z,M) {w<-z<=M; w/sum(w)}
    } else if (window=="clev"||window=="tric"){
      #Cleveland tricube window:
      weights <- function(z,M) {w<-z<=M; w<-w*(1-abs(z/M)^3)^3; w/sum(w)}
    } else if (window=="epan"||window=="para"){
      #Epanechnikov parabolic window:
      weights <- function(z,M) {w<-z<=M; w<-w*(1-abs(z/M)^2); w/sum(w)}
    } else if (window=="bart"||window=="tria"){
      #Bartlett triangular window:
      weights <- function(z,M) {w<-z<=M; w<-w*(1-abs(z/M)); w/sum(w)}
    } else if (window=="wend"){
      #Wendland 1995 window:
      weights <- function(z,M) {w<-z<=M; w<-w*(1-abs(z/M))^4*(1+4*abs(z/M)); w/sum(w)}
    } else if (window=="gaus"){
      #gauss window:
      weights <- function(z,M) {w<-1/exp((z/M)^2); w/sum(w)}
    } else {stop("wrong window type")}

    weighted <- function(x,s,M) {z<-abs( (1:length(x))-s ); w<-weights(z,M); w*x} #s=observation number; z=distance |t-s|
    sum.of.squares <- function(x) { sum(x^2, na.rm=TRUE) / sum(!is.na(x)) }
    sum.of.not.squares <- function(x) { sum(x, na.rm=TRUE) / sum(!is.na(x)) }
    d <- length(xx)       #number of series
    dd <- d*(d-1)/2       #number of correlations
    N <- nrow(xx) #number of observations

    val <- lo <- up <- YmaxR <- matrix(nrow=N)
    for (s in 1:N) {

      xy.cor <- vector("list", dd)
      jk <- 0
      for(k in 1:(d-1)) {
        for(j in (k+1):d) {
          jk <- jk+1
          x.w <- weighted( xx[[j]] ,s,M)
          y.w <- weighted( xx[[k]] ,s,M)
          xy.cor[[jk]] <- cor(x.w,y.w)
          if(sum(is.infinite(xy.cor[[jk]]))>0) browser()
        }}
      # browser()
      r <- unlist(xy.cor)
      if (is.na(sum(r))||is.nan(sum(r))){ xy.mulcor <- NA
      } else{
        P <- diag(d)/2
        P[lower.tri(P)] <- r
        P <- P+t(P)
        if (qr(P)$rank < d) {xy.mulcor <- 1; Pimax <- NA}
        else {
          Pidiag <- diag(solve(P))
          if(is.null(ymaxr)) {
            Pimax <- which.max(Pidiag) ## detect i | x[i] on rest x gives max R2
          } else {Pimax <- ymaxr}
          xy.mulcor <- sqrt(1-1/Pidiag[Pimax]) ## max(sqrt(1-1/diag(solve(P))))
          if (xy.mulcor>1) browser()
        }}
      #}
      val[s] <- xy.mulcor
      lo[s]  <- tanh(atanh(xy.mulcor)-qnorm(p)/sqrt(N-3))
      up[s]  <- tanh(atanh(xy.mulcor)+qnorm(p)/sqrt(N-3))
      YmaxR[s] <- Pimax
    }
    names(val) <- names(lo) <- names(up) <- names(YmaxR) <- names(xx[[1]])
    out <- list(val=as.data.frame(val),lo=as.data.frame(lo),up=as.data.frame(up),YmaxR=as.data.frame(YmaxR))
    return(out)
  }
#--------------------------------------------------------------
