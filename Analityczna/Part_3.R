#Pakiety
if(!require("mvtnorm")) {install.packages("mvtnorm"); library(mvtnorm)}

rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")

data <- read.csv("C:/Users/Marek/Documents/R_Scripts/Ekonometria_bayesowska/ogorki.csv")

data <- data[, c("Q_cuc", "P_cuc",  "P_tom", "P_oni")]

data_log <- log(data)
#data_log

model <- lm(Q_cuc ~ ., data = data_log)
summary(model)

N.data <- nrow(data)
X <- as.matrix(cbind(matrix(1, nrow = N.data, ncol = 1), data[, c("P_cuc",  "P_tom", "P_oni")]))
Beta.ols.data <- model$coefficients
v.data <- model$df.residual
XTX.data <- t(X) %*% X
s2.data <- sum((model$residuals) ^ 2) / v.data

#2. Parametry a priori
Beta.prior <- c(10, -0.57, -0.04, -0.11)
sm2.prior <- 16
U.prior <- diag(4) 
U.prior[1, 1] <- 1000
U.prior[2, 2] <- 0.39
U.prior[3, 3] <- 0.81
U.prior[4, 4] <- 0.81
v.prior <- 30

#3. Parametry a posteriori
Beta.posterior <- solve(solve(U.prior) + XTX.data) %*% (solve(U.prior) %*% Beta.prior + XTX.data %*% Beta.ols.data)
U.posterior <- solve(solve(U.prior) + XTX.data)
v.posterior <- v.prior + N.data
vs2.posterior <- v.prior / sm2.prior + v.data * s2.data + t(Beta.ols.data - Beta.prior) %*% solve(U.prior + solve(XTX.data)) %*% (Beta.ols.data - Beta.prior)
sm2.posterior <- 1 / (vs2.posterior / v.posterior)

#4. Parametry 
x.tau <- t(as.matrix(c(1,
                       3.321432,
                       3.990734,
                       2.921797)))  #obs nr 20

predictive.mean <- x.tau %*% Beta.posterior
predictive.scale <- 1 / sm2.posterior * (1 + x.tau %*% U.posterior %*% t(x.tau))
predictive.df <- v.posterior

#5. Gestosc predykcyjna
y.space <- seq(from = 1, to = 15, by = 0.1)
n_eval_points <- length(y.space)
library(metRology)
posterior.pred.dens.y <- dt.scaled(y.space, 
                                   mean = predictive.mean, sd = predictive.scale ^ 0.5, df = predictive.df)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)

plot(y.space, posterior.pred.dens.y, las = 1, lwd = 2, bty = "n", col = green_area,
     ylim = c(0, max(c(max(posterior.pred.dens.y), max(posterior.pred.dens.y))) + 0.00001), type = "l", 
     ylab = "gestosc predykcyjna", main = "gestosc predykcyjna")

#6. Predykcja punktowa
abline(v = predictive.mean,col = green_line, lwd = 3)
text(x = predictive.mean, y = max(posterior.pred.dens.y), 
     paste("prognoza punktowa = ", round(predictive.mean, digits = 4)), col = green_line)

#7. Predykcja przedzia?owa
lower90 <- qt.scaled(0.05, df = predictive.df, mean = predictive.mean, sd = (predictive.scale) ^ 0.5, lower.tail = TRUE, log.p = FALSE)
upper90 <- qt.scaled(0.95, df = predictive.df, mean = predictive.mean, sd = (predictive.scale) ^ 0.5, lower.tail = TRUE, log.p = FALSE)
print(paste("90-procentowy przedzial ufno?ci dla prognozy: ", round(lower90, 2), "-", round(upper90, 2)))

lower50 <- qt.scaled(0.25, df = predictive.df, mean = predictive.mean, sd = (predictive.scale) ^ 0.5, lower.tail = TRUE, log.p = FALSE)
upper50 <- qt.scaled(0.75, df = predictive.df, mean = predictive.mean, sd = (predictive.scale) ^ 0.5, lower.tail = TRUE, log.p = FALSE)
print(paste("50-procentowy przedzial ufno?ci dla prognozy: ", round(lower50, 2), "-", round(upper50, 2)))

