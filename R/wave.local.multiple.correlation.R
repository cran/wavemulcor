#----------------------------------------------------------
wave.local.multiple.correlation <-
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
    l <- length(xx[[1]])  #number of scales J+1 (wavelet coefficients at levels 1 to J plus the scaling coeffs at level J+1)
    N <- length(xx[[1]][[1]]) #number of observations
    wav <- attr(xx[[1]],"wavelet")
    if(substr(wav,1,1)=="h") {L<-2
    } else if(substr(wav,1,1)=="d") {L<-as.numeric(substring(wav,2))
    } else if(substr(wav,1,1)=="l") {L<-as.numeric(substring(wav,3))} else L<-8

    ###IMPORTANT: must shift wave coeffs to left by L_i/2 at each level to be in phase with time series!!! (Getal, p.145)
    for(j in 1:d) {
      xx[[j]] <- phase.shift(xx[[j]], wav)
	}

    val <- lo <- up <- YmaxR <- list()

    cat("\nlev:")
    for(i in 1:l) { cat(sprintf("%s",i))
      val[[i]] <- lo[[i]] <- up[[i]] <- YmaxR[[i]] <- matrix(nrow=N)
      for (s in 1:N) {

        x.var <- vector("list", d)
        for(j in 1:d) {
          x.w <- weighted( xx[[j]][[i]] ,s,M)  #xx.w[[j]][[i]][[s]]
          x.var[[j]] <- sum.of.squares(x.w)
        }

        xy.cor <- vector("list", dd)
        jk <- 0
        for(k in 1:(d-1)) {
          for(j in (k+1):d) {
            jk <- jk+1
            x.w <- weighted( xx[[j]][[i]] ,s,M)  #xx.w[[j]][[i]][[s]]
            y.w <- weighted( xx[[k]][[i]] ,s,M)  #xx.w[[k]][[i]][[s]]
            xy.w <- as.vector( x.w * y.w )
            xy.cov <- sum.of.not.squares(xy.w)
            xy.cor[[jk]] <- xy.cov / sqrt(x.var[[j]] * x.var[[k]])
            if(sum(is.infinite(xy.cor[[jk]]))>0) browser()
          }}
        r <- unlist(xy.cor)
        if (is.na(sum(r))||is.nan(sum(r))){ xy.mulcor[[i]] <- NA
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
        #not Ntilde but Nhat ### NOTE!!!! as in waveslim for DWT and approx for MODWT (Getal p.241,260)
        Lprime <- ceiling( (L-2)*(1-1/2^i) )
        Nhat <- 2^(ceiling(log2(N))-i)-Lprime
        val[[i]][s] <- xy.mulcor
        lo[[i]][s]  <- tanh(atanh(xy.mulcor)-qnorm(p)/sqrt(Nhat-3))
        up[[i]][s]  <- tanh(atanh(xy.mulcor)+qnorm(p)/sqrt(Nhat-3))
        YmaxR[[i]][s] <- Pimax
      }}
    names(val) <- names(lo) <- names(up) <- names(YmaxR) <- names(xx[[1]])
    out <- list(val=as.data.frame(val),lo=as.data.frame(lo),up=as.data.frame(up),YmaxR=as.data.frame(YmaxR))
    return(out)
  }
#--------------------------------------------------------------
