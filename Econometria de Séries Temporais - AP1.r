# @hidden_cell
options(warn=-1)
rm(list=ls())

library("tseries")
library(forecast)
library(ggplot2)
library(readr)
theme_set(theme_minimal())

################################################################################
## LOAD SERIES                                                                ##
################################################################################

series <- read_csv("series.csv",col_types = cols())

serie <- series$SER19
nobs <- length(serie)
nyears <- floor(nobs/12)
nmonths <- nobs -  nyears*12

endYear <- 2018
endMonth <- 12

## transforma série em timeseries do R - para 90 obs 7,5 anos
#tserie <- ts(serie, start=c(endYear - nyears, endMonth - nmonths+1), end=c(endYear, endMonth), frequency=12)
tserie <- ts(serie,  frequency=1)

#tserie <- ts(serie, frequency=1)
s <- summary(tserie)
u <-mean(serie)
sd <- sd(serie)

options(repr.plot.width=8, repr.plot.height=4)

ts.plot(tserie, xlab="Tempo", ylab="", main="Série 19")
abline(h = u, col="red", lwd=1, lty=1)
abline(h = u+sd, col="blue", lwd=1, lty=2)
abline(h = u-sd, col="blue", lwd=1, lty=2)

abline(v=12,col="green", lwd=1, lty=3)
abline(v=24,col="green", lwd=1, lty=3)
abline(v=30,col="green", lwd=1, lty=3)
abline(v=42,col="green", lwd=1, lty=3)

print(s)

################################################################################
## TESTE ESTACIONARIEDADE                                                     ##
################################################################################

adf.test(tserie)
kpss.test(tserie)



################################################################################
## FUNCOES AUTO CORRELACAO                                                    ##
################################################################################

tserie.acf <- acf(tserie,lag.max=24)
tserie.pacf <- pacf(tserie)


################################################################################
## ARMA - AR(1,2) && MA(1)                                                    ##
################################################################################

for (p in c(0,1,2)){
  for (q in c(0,1)){ 
    ARMAN <- arima(tserie, order = c(p,0,q))

    cat(sprintf("Verificando ARIMA(%i, 0, %i) ", p, q))     
    cat(sprintf("AIC: %f\t",AIC(ARMAN)))
    cat(sprintf("BIC: %f\n",BIC(ARMAN)))
    
    cat(sprintf('p-value:\n'))
    print(
      (1-pnorm(abs(ARMAN$coef)/sqrt(diag(ARMAN$var.coef))))*2)

    cat(sprintf('R^2 = %f\n\n', cor(fitted(ARMAN),serie)^2))
      
      
    cat(sprintf('--------------------------------\n'))
    
  }
}

################################################################################
## MELHOR MODELO ARIMA(0,0,1) ou MA(1)                                        ##
################################################################################

ARMA <- arima(tserie, order = c(0,0,1))
ARMA_fit <- tserie - residuals(ARMA)

ts.plot(tserie, xlab="Tempo", ylab="", main="Série MA(1) ")
points(ARMA_fit, type = "l", col = 2, lty = 2)


acf(residuals(ARMA))
pacf(residuals(ARMA))


################################################################################
## MELHOR MODELO ARIMA(0,0,5) ou MA(5)                                        ##
################################################################################

MA_6 <- arima(tserie, order = c(0,0,6))
MA_6_fit <- tserie - residuals(MA_6)

cat(sprintf("Verificando ARIMA(0, 0, 5) "))     
cat(sprintf("AIC: %f\t",AIC(MA_6)))
cat(sprintf("BIC: %f\n",BIC(MA_6)))

cat(sprintf('p-value:\n'))
print(
  (1-pnorm(abs(MA_6$coef)/sqrt(diag(MA_6$var.coef))))*2)

cat(sprintf('R^2 = %f\n\n', cor(fitted(MA_6),serie)^2))
      
      
cat(sprintf('--------------------------------\n'))




ts.plot(tserie, xlab="Tempo", ylab="", main="Série MA(6) ")
points(MA_6_fit, type = "l", col = 2, lty = 2)

     
acf(residuals(MA_6))
pacf(residuals(MA_6))

################################################################################
## MA(6) x MA(1)                                                              ##
################################################################################
ts.plot(tserie, xlab="Tempo", ylab="", main="Série MA(1) x MA(6) ")
points(ARMA_fit,  type="l",col = 2 ,lty=2)
points(MA_6_fit, type = "o", pch=3, col =4, lty = 2)

tserie.training <- subset(tserie, end=60) 
tserie.test <- subset(tserie, start=61) 



for (p in c(0,1,2)){
  for (q in c(0,1)){ 
    cat(sprintf("Verificando ARIMA(%i, 0, %i)\n\n ", p, q))    
    tserie.trainingModel <- arima(tserie.training, order = c(p,0,q))
    tserie.testModel <- Arima(tserie.test, model=tserie.trainingModel)
    print(accuracy(tserie.testModel))
    
    cat(sprintf('--------------------------------\n'))
    
  }
}

tserie.trainingModel <- Arima(tserie.training, order=c(0,0,1))
tserie.testModel <- Arima(tserie.test, model=tserie.trainingModel)
autoplot(tserie.test, series="Conjunto de Teste") +
  autolayer(fitted(tserie.testModel, h=1),
            series="Valores previstos")

################################################################################
## RESIDUOS DE MA(6)                                                          ##
################################################################################


ts.plot(residuals(MA_6), xlab="Tempo", ylab="", main="Série resíduos ")
ts.plot(residuals(MA_6)^2, xlab="Tempo", ylab="", main="Série resíduos^2")



acf(residuals(MA_6))  ##NAO HA EFEITO ARCH
acf(residuals(MA_6)^2)  ##NAO HA EFEITO ARCH
pacf(residuals(MA_6))  ##NAO HA EFEITO ARCH


shapiro.test(residuals(MA_6))
mean(residuals(MA_6))

checkresiduals(MA_6)


tserie.garch <- garch(tserie, c(1,1))
summary(tserie.garch) 



ht.arch08=tserie.garch$fit[,1]^2
plot(ht.arch08,xlab="Tempo", ylab="", main='Variância condicional')

