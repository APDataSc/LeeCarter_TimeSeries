
rm( list = ls() )
graphics.off( )
options(scipen=999)

# Lista de paquetes 
.packages <-  c( "data.table",  
                 "tidyverse", 
                 "lubridate", 
                 "ggplot2",  
                 "MortalityLaws",
                 "fpp2",
                 "urca")

# Instalacion de paquetes no instalados
.inst <- .packages %in% installed.packages()
if( length( .packages[ !.inst ] ) > 0 ){
  install.packages( .packages[ !.inst ], dependencies = TRUE )
}

# Carga de paquetes
lapply( .packages, require, character.only = TRUE )

## lectura datos 
lt_ecu_f <- 
  fread('dat/ECU_TABMORT_LCLIM_X1_1950_2020.csv') %>%
  .[ sex=="f" & year!=2020 , .( year, date_ref, sex, age, mx, ex) ]


# Next we compute the log of the rates. We will store them in a matrix of 55 
# years by 19 age groups.
rates <- lt_ecu_f$mx
M <- matrix(log(rates), 70, 101, byrow = TRUE)

"Fitting the Model"
# The first thing we need is the mean log-rate for each age group. This is easily 
# computed using colMeans(). The results are remarkably close to Table 1 in the original paper.
a <- colMeans(M)
a
plot(a)

# The next step is to subtract the average age pattern a from all years.
for(j in 1:101) M[,j] <- M[,j] - a[j]

# We are now ready to compute the Singular Value Decomposition (SVD), which writes 
# M = U D V, where U and V’ are orthogonal matrices and D is a diagonal matrix of 
# singular values. In fact we need just the first left and right singular vectors. 
# The first column of U times D1,1 times the first row of V’ has the best rank-1 
# approximation to the input matrix.

d <- svd(M, 1, 1)
plot(d$d/sum(d$d)*100, type="l")

# Lee and Carter normalize the first row of V so it sums to one and call it b. 
# This vector models how the different age groups react to mortality change.

b <- d$v/sum(d$v)
head(b, 5) 
plot(b)

# These values are remarkably close to the b’s published in Table 2. Lee and Carter 
# also take the first column of U, multiply by D1,1 and multiply by the sum of the 
# first row of V’ (to cancel the division) and call that k. This vector captures 
# overall mortality change over time.

k <- d$u * sum(d$v) * d$d[1]
head(k, 5)
tail(k, 5)

plot(k)



"Plotting Parameters and Fits"
# The next task is to compute the Lee-Carter fits to the mortality rates in 1933 
# and 1987, reproducing part of Figure 2 in the original paper.

usmx2 <- filter(lt_ecu_f, year==1950 | year == 1985 | year == 2019) |> 
  mutate(year = factor(year),
         fit = c(exp(a + b * k[1]), a + b * k[35], exp(a + b * k[70])))

ggplot(usmx2, aes(age, mx, color=year)) + geom_point() + 
  geom_line(aes(age,fit,color=year)) + scale_y_log10() +   
  ggtitle("Lee-Carter Fits for 1933 and 1987")

# Here’s the trajectory of k
trend <- data.frame(year = 1950:2019, k = k)
ggplot(trend, aes(year, k)) + geom_line() +
  ggtitle("Lee-Carter k for 1933-1987")


"ARIMA Model"

k <- ts(k)
autoplot(k)

k %>% ur.kpss() %>% summary()

k %>% diff() %>% ur.kpss() %>% summary()

ndiffs(k)

fit <- auto.arima(k, seasonal=FALSE)

fit %>% forecast(h=30) %>% autoplot(include=80)

ggAcf(k)
ggPacf(k)

k %>% diff() %>% ggtsdisplay(main="")

checkresiduals(fit)

autoplot(fit)

forecast_k <- fit %>% forecast(h=30)
k_f <- forecast_k$mean[1:30]  

## function to compute LC log-rates from LC parameters (Lee & Carter, 1992)
LCeta <- function(ax,bx,kt){
  n <- length(kt)
  One <- matrix(1,n,1)
  ETA <- ax%*%t(One) + bx%*%t(kt)
  return(ETA)
}


LMX.lc1 <- LCeta(a, b, k)
LMX.lcf <- LCeta(a, b, k_f)


matplot(unique(lt_ecu_f$age), LMX.lc1,t="l", lty=1, 
        col=rainbow(ncol(M)), ylab = expression(ln(m[x])), 
        xlab = "x", ylim = c(-10,0))

matplot(unique(lt_ecu_f$age), LMX.lcf,t="l", lty=1, 
         ylab = "ln(m_x)", 
        xlab = "Age", add = T)

## observed (log-)rates
MX.obs <- matrix(rates,101,70)
LMX.obs <- log(MX.obs)


## function for constructing a classic (& rather general) lifetable
## from mortality rates (Chapter 3 in Preston et al. 2001)
lifetable.mx <- function(x, mx, sex="M", ax=NULL){
  m <- length(x)
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

## function to derive life expectancy (e0) 
## from mx (based on previous function)
e0.mx <- function(x, mx, sex="M", ax=NULL){
  lt <- lifetable.mx(x,mx, sex,ax)
  return.ex <- lt$ex[1]
  return(return.ex)
}


## computing e0 for each year
e0obs <- apply(MX.obs, 2, e0.mx, x=0:100)
e0lc <- apply(exp(LMX.lc1), 2, e0.mx, x=0:100)

## plot observed and LC e0
plot(1950:2019, e0obs,t="p",ylim=range(e0obs,e0lc))
lines(1950:2019,e0lc,col=2)


## computing PIs of kt
sd <- sqrt(fit$sigma2) 
tF <- 1950:2019

c.val <- qnorm(0.975)
kt.up <- k + c.val*sd
kt.low <- k - c.val*sd
plot(tF,k,t="l",ylim = range(kt.up,kt.low))
lines(tF,kt.up,lty=2,col=2)
lines(tF,kt.low,lty=2,col=2)

## compute lower and upper log-rates 
LMX.lc.up <- LCeta(a, b, kt.up)
LMX.lc.low <- LCeta(a, b, kt.low)

## compute lower and upper e0
e0lc.up <- apply(exp(LMX.lc.low),2,e0.mx,x=0:100)
e0lc.low <- apply(exp(LMX.lc.up),2,e0.mx,x=0:100)

plot(tF,e0lc,t="l",ylim = range(e0lc.up,e0lc.low))
lines(tF,e0lc.up,lty=2,col=2)
lines(tF,e0lc.low,lty=2,col=2)
lines(tF,e0obs,lty=2,col=4)


## panel 1: LE
t <- 1950:2019
cex.main <- 1.65
col.obs <- "grey30"
col.lc <- 4
col.lcT <- adjustcolor(col.lc,alpha.f = 0.3)
cex.axis.title <- 1.4

plot(t,e0obs,t="n",ylim=ylimE0,axes=F,xlab="",ylab="")
axis(1,padj = -0.5);axis(2,las=2,hadj = .75);grid();box()
title("Life expectancy at birth",cex.main=cex.main,line = 0.35)
points(t,e0obs,pch=16,col=col.obs)
xx <- c(tF,rev(tF))
yy <- c(e0lc.up,rev(e0lc.low))
polygon(xx,yy,col=col.lcT,border = col.lcT)
lines(tF,e0lc,col=col.lc,lwd=2)
mtext("year", 1, line=1.75,cex=cex.axis.title)
mtext(expression(e[0]), 2, line=1.5,cex=cex.axis.title,las=2)

legend("topleft",legend=c("Observed","Lee-Carter (1992)"),lwd=2,cex=1.25,
       pch=c(16,NA),lty=c(NA,1),
       col=c(col.obs,col.lc),bg="white",inset=0.05)
legend("topleft",legend=c("","", ""),pch=c(NA,15),col=c(NA,col.lcT),
       bty="n",lwd=NA,cex=1.25,lty=NA,inset=0.05)


#Comentarios adicionales
