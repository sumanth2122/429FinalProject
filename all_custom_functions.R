#################################################################

check.mvnorm.plot <- function(data) {
  data <- as.matrix(data)
  d2 <- 1:nrow(data)
  for (i in 1:nrow(data)) {
    d2[i] <- t(as.matrix(data[i,] -
                           apply(data,2,mean)))%*%solve(cov(data))%*%(as.matrix(data[i,]-apply(data,2,mean)))
  }
  p <- ((1:nrow(data))-.5)/nrow(data)
  plot(qchisq(p,ncol(data)),sort(d2),ylab="D^2 Values",xlab="Theoretical Chi-Square Quantiles",
       main="Chi-Square Probability Plot")
  abline(0,1)
}

#################################################################

one.sample.hotelling <- function(Data,mu.0){
  #### One Sample Hotelling's T squared test
  #### Data: matrix of the 2 vectors you want tested
  #### mu.0:vector of length p, with hypothesized values of mu
  Data=Data
  n=nrow(Data)
  p=ncol(Data)
  S <- cov(Data)
  y.bar <- apply(Data,2,mean)
  ybar.minus.mu <- y.bar - mu.0
  T.square <- n*t(ybar.minus.mu)%*%solve(S)%*%ybar.minus.mu
  F.stat <- round(((n-p)/(p*(n-1)))*T.square,3)
  P.value <- round(1-pf(F.stat,p,n-p),digits=5)
  Output <- paste("F=",F.stat,"P-Value=", P.value,sep=" ")
  return(Output)
}

#################################################################

two.sample.hotelling <- function(formula,data){
  as.formula(formula)
  as.matrix(data)
  #### formula: Use y~x format, where y is a matrix of variables
  #### whose means are to be tested, and x represents a factor
  #### data: Name of the data frame
  manova.obj <- manova(formula,data=data)
  return(summary(manova.obj))
}

#################################################################

paired.hotelling = function(Data){
  #### Data: Data matrix-all y variables first p/2 columns, all x variables [(p/2)-p]
  n=nrow(Data)
  col=ncol(Data)
  half=col/2
  D=matrix(0,n,half)
  for(i in 1:half)
  {
    d <- Data[,i]-Data[,half+i]
    D[,i]=d
  }
  paired.Tsq <- manova(D~1)
  return(summary(paired.Tsq,intercept=T,test="Hotelling"))
}

#################################################################

k.sample.profile = function(Data,Test="P",K=1)
{
  ### Data: first column categorical variable, next columns, sequential measurements
  ### Test: which type of test, p=parallelism,c=coincidence,f=flatness
  ### K: number of groups, default is 1
  nrow=nrow(Data)
  ncol=ncol(Data)
  if(K==1) {
    Data <- as.matrix(Data)
    C=matrix(0,ncol-1,ncol)
    for(itor in 1:ncol-1){
      for(itorat in 1:ncol){
        if(itorat-itor==0){
          C[itor,itorat]=1
        }
        if(itorat-itor==1){
          C[itor,itorat]=-1
        }}}
    trans.data <- Data%*%t(C)
    one.samp.profile <- manova(trans.data~1)
    output=(summary(one.samp.profile,intercept=T,test="Hotelling"))
  }
  if((Test=="P" | Test=="p") & K!= 1){
    crow=ncol-2
    ccol=ncol-1
    C=matrix(0,crow,ccol)
    for(itor in 1:crow){
      for(itorat in 1:ccol){
        if(itorat-itor==0){
          C[itor,itorat]=1
        }
        if(itorat-itor==1){
          C[itor,itorat]=-1
        }}}
    Y <- as.matrix(Data[,2:ncol])
    trans.data <- Y%*%t(C)
    two.samp.profile <- manova(trans.data~as.factor(Data[,1]))
    output=summary(two.samp.profile,test="Wilks")
  }
  if((Test=="C" | Test=="c") & K!= 1){
    j <- rep(1,times=ncol-1)
    Y <- as.matrix(Data[,2:ncol])
    trans.data <- Y%*%j
    coin.test <- aov(trans.data~as.factor(Data[,1]))
    output=summary(coin.test,test=F)
  }
  if((Test=="F" | Test=="f") & K!= 1){
    crow=ncol-2
    ccol=ncol-1
    C=matrix(0,crow,ccol)
    for(itor in 1:crow){
      for(itorat in 1:ccol){
        if(itorat-itor==0){
          C[itor,itorat]=1
        }
        if(itorat-itor==1){
          C[itor,itorat]=-1
        }}}
    Y <- as.matrix(Data[,2:ncol])
    trans.data <- Y%*%t(C)
    test.flat <- manova(trans.data~1+as.factor(Data[,1]))
    output=summary(test.flat,intercept=T,test="Hotelling")
  }
  return(output)
}

#################################################################

test.specific <- function(S, sig.0, N, k=1){
  # S is the sample covariance matrix or the pooled covariance matrix, S.pl
  # sig.0 is the hypothesized covariance matrix
  # N is the sample size of either the one group, or the total
  # sample size of all of the groups combined, in the case of S.pl
  # k is the number of groups, and is by default 1.
  p <- ncol(S)
  nu <- N - k
  evals <- eigen(S %*% solve(sig.0))$values
  u <- nu * (sum(evals - log(evals)) - p)
  u.prime <- (1-(1/(6*nu -1))*(2*p + 1 - 2/(p+1))) * u
  chi.df <- .5*p*(p+1)
  p.val.u <- 1 - pchisq(u, chi.df)
  p.val.u.prime <- 1 - pchisq(u.prime, chi.df)
  return(list(u = u, u.prime = u.prime, p.value.u = p.val.u,
              p.value.u.prime = p.val.u.prime, df = chi.df))
}

#################################################################

Box.M <- function(samp.covs,n){
  # n is a vector of sample sizes
  # samp.covs is a list containing sample covariance matrices
  k <- length(n)
  pp <- dim(samp.covs[[1]])
  df <- n-1
  S.dets <- 1:k
  Sp <- matrix(nrow=pp[1],ncol=pp[1],0)
  for (i in 1:k) {
    S.dets[i] <- det(samp.covs[[i]])
    Sp <- Sp + (df[i]*samp.covs[[i]])
  }
  Sp <- Sp/sum(df)
  det.Sp <- det(Sp)
  M <- prod((S.dets^(df/2)))/(det.Sp^(.5*sum(df)))
  c <- (sum(1/df)-1/sum(df))*((2*pp[1]^2+3*pp[1]-1)/(6*(pp[1]+1)*(k-1)))
  u <- -2*(1-c)*log(M) #(Using 7.23, p.257)
  u.df <- .5*(k-1)*pp[1]*(pp[1]+1)
  p.value <- 1-pchisq(u,u.df)
  return(list("u"=u,"p.value.u"=p.value,"chi.df"=u.df))
}

#################################################################

discrim <- function(Y, group){
  Y <- data.matrix(Y)
  group <- as.factor(group)
  m1 <- manova(Y ~ group)
  nu.h <- summary(m1)$stats[1]
  nu.e <- summary(m1)$stats[2]
  p <- ncol(Y)
  SS <- summary(m1)$SS
  E.inv.H <- solve(SS$Residuals) %*% SS$group
  eig <- eigen(E.inv.H)
  s <- min(nu.h, p)
  lambda <- Re(eig$values[1:s])
  a <- Re(eig$vectors[,1:s])
  a.star <- (sqrt(diag(SS$Residuals/nu.e)) * a)
  return(list("a"=a, "a.stand"=a.star))
}

#################################################################

discr.sig <- function(Y, group){
  Y <- data.matrix(Y)
  group <- as.factor(group)
  m1 <- manova(Y ~ group)
  sums <- summary(m1)
  evals <- sums$Eigenvalues
  nu.e <- m1$df
  nu.h <- m1$rank-1
  k <- nu.h + 1
  p <- ncol(m1$coef)
  N <- nu.e + nu.h + 1
  s <- min(p, nu.h)
  lam <- numeric(s)
  dfs <- numeric(s)
  for(m in 1:s){
    lam[m] <- prod(1/(1+evals[m:s]))
    dfs[m] <- (p-m+1)*(k-m)
  }
  V <- -(N - 1 - .5*(p+k))*log(lam)
  p.val <- 1 - pchisq(V, dfs)
  out <- cbind(Lambda=lam, V, p.values=p.val)
  dimnames(out)[[1]] <- paste("LD",1:s,sep="")
  return(out)
}

#################################################################

partial.F <- function(Y, group){
  Y <- data.matrix(Y)
  group <- as.factor(group)
  p <- ncol(Y)
  m1 <- manova(Y ~ group)
  nu.e <- m1$df
  nu.h <- m1$rank-1
  Lambda.p <- summary(m1,test="Wilks")$stats[3]
  Lambda.p1 <- numeric(p)
  for(i in 1:p){
    dat <- data.matrix(Y[,-i])
    m2 <- manova(dat ~ group)
    Lambda.p1[i] <- summary(m2,test="Wilks")$stats[3]
  }
  Lambda <- Lambda.p / Lambda.p1
  F.stat <- ((1 - Lambda) / Lambda) * ((nu.e - p + 1)/nu.h)
  p.val <- 1 - pf(F.stat, nu.h, nu.e - p + 1)
  out <- cbind(Lambda, F.stat, p.value = p.val)
  dimnames(out)[[1]] <- dimnames(Y)[[2]]
  ord <- rev(order(out[,2]))
  return(out[ord,])
}

#################################################################

discr.plot <- function(Y, group, leg = NULL){
  a <- discrim(Y, group)$a
  z <- data.matrix(Y) %*% a
  plot(z[,1], z[,2], type = "n", xlab = "LD1", ylab = "LD2")
  for(i in 1:length(unique(group))){
    points(z[group == unique(group)[i],1],
           z[group == unique(group)[i],2], pch = i)
  }
  if(is.null(leg)) leg <- as.character(unique(group))
  legend("topright", legend = leg, pch = 1:length(unique(group)))
}

#################################################################

lin.class <- function(Y,group){
  # Install MASS package if not already installed
  if (!require("MASS")) install.packages("MASS")
  library(MASS)
  Y <- data.matrix(Y)
  group <- as.factor(group)
  p <- ncol(Y)
  m1 <- manova(Y ~ group)
  nu.e <- m1$df
  nu.h <- m1$rank-1
  Sp <- summary(m1)$SS$Residual/(nu.e)
  cio <- 1:m1$rank
  c.mat <- matrix(nrow=m1$rank,ncol=p,0)
  for (i in 1:m1$rank) {
    cio[i] <- -.5*t(lda(Y,group)$means[i,])%*%solve(Sp)%*%
      lda(Y,group)$means[i,]
    c.mat[i,] <- t(lda(Y,group)$means[i,])%*%solve(Sp)
  }
  return(list("coefs"=c.mat,"c.0"=cio))
}

#################################################################

rates <- function(data,group,method="l") {
  if (!require("MASS")) install.packages("MASS")
  library(MASS)
  data <- as.matrix(data)
  group <- as.matrix(group)
  da.obj <- lda(data,group)
  if (method=="q") {
    da.obj <- qda(data,group)
    method <- "QDA"
  }
  tab <- table(original=group,predicted=predict(da.obj)$class)
  if (method=="l") method <- "LDA"
  cor.rate <- sum(predict(da.obj)$class==group)/nrow(data)
  er.rate <- 1-cor.rate
  return(list("Correct Class Rate"=cor.rate,"Error Rate"=er.rate,
              "Method"=method,"Confusion Matrix"=tab))
}

#################################################################

cancor.test <- function(x,y){
  r <- cancor(x,y)$cor
  p <- ncol(y)
  q <- ncol(x)
  n <- nrow(x)
  # Wilk's Lambda calculations
  L <- prod(1-r^2)
  t.num <- (p^2)*(q^2) - 4
  t.den <- p^2 + q^2 - 5
  t <- sqrt(t.num/t.den)
  w <- n - .5*(p + q +3)
  F.df1 <- p*q
  F.df2 <- w*t - .5*p*q + 1
  L1t <- L^(1/t)
  F.stat <- ((1-L1t)/L1t)*(F.df2/F.df1)
  p.val.F <- 1 - pf(F.stat, F.df1, F.df2)
  L.sum <- cbind(signif(L,4), round(F.stat,3), df1=F.df1,
                 round(F.df2,2), round(p.val.F,4))
  # Pillai's V(s) calculations
  V <- sum(r^2)
  s <- min(p,q)
  m <- .5*(abs(q-p) -1)
  N <- .5*(n - q - p -2)
  F.stat <- (2*N + s + 1)*V / ((2*m + s + 1)*(s-V))
  F.df1 <- s * (2*m+s+1)
  F.df2 <- s * (2*N+s+1)
  p.val.F <- 1 - pf(F.stat, F.df1, F.df2)
  V.sum <- cbind(signif(V,4), round(F.stat,3), F.df1, F.df2,
                 round(p.val.F, 4))
  # Lawley-Hotelling U(s) calculations
  U <- sum(r^2 / (1-r^2))
  F.stat <- 2*(s*N+1)*U/(s^2*(2*m + s + 1))
  F.df1 <- s*(2*m + s+ 1)
  F.df2 <- 2*(s*N + 1)
  p.val.F <- 1 - pf(F.stat, F.df1, F.df2)
  U.sum <- cbind(signif(U,4), round(F.stat,3), F.df1, F.df2,
                 round(p.val.F, 4))
  theta <- r[1]^2
  F.stat <- NA
  F.df1 <- NA
  F.df2 <- NA
  p.val.F <- NA
  theta.sum <- cbind(signif(theta,4), round(F.stat,3),
                     F.df1, F.df2,
                     round(p.val.F, 4))
  out <- rbind(L.sum, V.sum, U.sum, theta.sum)
  dimnames(out)[[1]] <- c("Wilks Lambda", "Pillai V(s)",
                          "Lawley-Hotelling U(s)", "Roy's Largest Root")
  dimnames(out)[[2]] <- c("Statistic", "Approx F", "df1", "df2",
                          "p-value")
  return(out)
}

#################################################################

cancor.step.test <- function(x,y){
  r <- cancor(x,y)$cor
  p. <- ncol(y)
  q. <- ncol(x)
  n. <- nrow(x)
  s <- min(p.,q.)
  out <- numeric(5)
  for(k in 1:s){
    L <- prod(1-r[k:s]^2)
    p <- p. - k + 1
    q <- q. - k + 1
    n <- n. - k + 1
    t.num <- (p^2)*(q^2) - 4
    t.den <- p^2 + q^2 - 5
    t <- sqrt(t.num/t.den)
    w <- n - .5*(p + q +3)
    F.df1 <- p*q
    F.df2 <- w*t - .5*p*q + 1
    L1t <- L^(1/t)
    F.stat <- ((1-L1t)/L1t)*(F.df2/F.df1)
    p.val.F <- 1 - pf(F.stat, F.df1, F.df2)
    L.sum <- cbind(signif(L,4), round(F.stat,3), F.df1,
                   round(F.df2,2), round(p.val.F,4))
    out <- rbind(out, L.sum)
  }
  out <- out[-1,]
  out <- cbind(1:s, out)
  dimnames(out)[[2]] <- c("k", "Wilk's Lambda", "Approx F",
                          "df1", "df2", "p-value")
  return(out)
}

#################################################################

princomp.test <- function(dataset){
  lambda <- prcomp(dataset)$sdev^2
  n <- nrow(dataset)
  p <- length(lambda)
  mult <- n - ((2*p + 11)/6)
  u <- matrix(0,p-1)
  df <- matrix(0,p-1)
  p.val <- matrix(0,p-1)
  for(k in 2:p){
    l.bar <- mean(lambda[(p-k+1):p])
    u[k-1] <- mult * (k*log(l.bar) - sum(log(lambda[(p-k+1):p])))
    df[k-1] <- .5*(k-1)*(k+2)
    p.val[k-1] <- 1 - pchisq(u[k-1], df[k-1])
  }
  out <- cbind(lambda[1:(p-1)], p:2, rev(u), rev(df), rev(p.val))
  dimnames(out)[[2]] <- c("Eigenvalue", "k", "u","df","p-value")
  return(out)
}

#################################################################

