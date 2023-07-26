rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")

data <- read.csv("C:/Users/Marek/Documents/R_Scripts/Ekonometria_bayesowska/ogorki.csv")

data <- data[, c("Q_cuc", "P_cuc",  "P_tom", "P_oni")]
### funkcja eval nie lubi spacji w nazwach zmiennych, dla obliczenia czynnika bayesa trzeba _
#names(data) <- c("Wolumen ogorki", "Cena ogorki", "Cena pomidory", "Cena cebula")
data_log <- log(data)
model <- lm(Q_cuc ~ ., data = data_log)
summary(model)



if(!require("manipulate")) {install.packages("manipulate"); library(manipulate)}
if(!require("mvtnorm")) {install.packages("mvtnorm"); library(mvtnorm)}

y <- as.matrix(data$Q_cuc)
N.data <- length(y)
X <- cbind(as.matrix(rep(1, N.data)), 
           as.matrix(data[, c("P_cuc", "P_tom", "P_oni")]))
Beta.ols.data <- model$coefficients
v.data <- model$df.residual
XTX.data <- t(X) %*% X
s2.data <- sum((model$residuals) ^ 2) / v.data

#3. Parametry a priori

Beta.prior <- c(10, -0.57, -0.04, -0.11)
sm2.prior <- 16
U.prior <- diag(4) 
U.prior[1, 1] <- 1000
U.prior[2, 2] <- 0.39
U.prior[3, 3] <- 0.81
U.prior[4, 4] <- 0.81
v.prior <- 30


#4. Parametry a posteriori

Beta.posterior <- solve(solve(U.prior) + XTX.data) %*% (solve(U.prior) %*% Beta.prior + XTX.data %*% Beta.ols.data)
U.posterior <- solve(solve(U.prior) + XTX.data)
v.posterior <- v.prior + N.data
vs2.posterior <- v.prior / sm2.prior + v.data * s2.data + t(Beta.ols.data-Beta.prior) %*% solve(U.prior + solve(XTX.data)) %*% (Beta.ols.data - Beta.prior)
sm2.posterior <- 1 / (vs2.posterior / v.posterior)
vs2.prior <- v.prior / sm2.prior

n_var <- 4

if (n_var == 2){
  beta.space <- seq(from = -1.1, to = -0.7, by = 0.001)
}

if (n_var == 3){
  beta.space <- seq(from = -0.87, to = -0.75, by = 0.001)
}

if (n_var == 4){
  beta.space <- seq(from = -0.7, to = -0.55, by = 0.001)
}

n_eval_points <- length(beta.space)
n_parameters <- length(Beta.posterior)
prior.marg.dens.beta <- matrix(NA, nrow = n_parameters, ncol = n_eval_points)
posterior.marg.dens.beta <- matrix(NA, nrow = n_parameters, ncol = n_eval_points)
for(ii in 1:length(Beta.posterior)) {
  prior.marg.dens.beta[ii, ] <- apply(as.matrix(beta.space), 1, dmvt,
                                      delta = Beta.prior[ii], sigma = as.matrix(U.prior[ii, ii] / sm2.prior), df = v.prior, log = FALSE)
  posterior.marg.dens.beta[ii, ] <- apply(as.matrix(beta.space), 1, dmvt,
                                          delta = Beta.posterior[ii], sigma = as.matrix(U.posterior[ii, ii] / sm2.posterior), df = v.posterior, log = FALSE)
}

# 6. Szkicujemy HPDI wok?? wybranego wsp??czynnika kierunkowego r?wnania regresji
###  Wybierz:
###  2 <- elastyczno?? cenowa popytu w?asna
###  3 <- elastyczno?? cenowa popytu krzy?owa wzgl?dem samochod?w
###  4 <- elastyczno?? dochodowa popytu

grey_area <- rgb(160, 160, 160, 80, names = NULL, maxColorValue = 255)
grey_line <- rgb(80, 80, 80, 160, names = NULL, maxColorValue = 255)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)
red_area <- rgb(255, 100, 123, 80, names = NULL, maxColorValue = 255)
red_line <- rgb(200, 0, 30, 160, names = NULL, maxColorValue = 255)

par(mfrow = c(1, 1))
manipulate( 
  { credible_set_indicator <- as.vector(as.integer(posterior.marg.dens.beta[n_var, ] > line_level))
    credible_set_begin <- match(1, credible_set_indicator)
    credible_set_end <- length(credible_set_indicator) - match(1, rev(credible_set_indicator))
    x1 <- beta.space[credible_set_begin]
    x2 <- beta.space[credible_set_end]
    posterior.cs <- posterior.marg.dens.beta[n_var, ] * credible_set_indicator
    HPDI_probab <- sum(posterior.cs) * 0.01
    
    plot(beta.space, posterior.marg.dens.beta[n_var, ], las = 1, lwd = 2, bty = "n", col = green_line,
         ylim = c(0, max(posterior.marg.dens.beta[n_var, ] + 1)), type = "l", ylab = "gestosc", main = colnames(data)[n_var])
    polygon(c(beta.space, rev(beta.space)), 
            c(posterior.marg.dens.beta[n_var, ], rep(0, length(beta.space))), 
            col = green_area, border = NA)
    text(Beta.posterior[n_var], max(posterior.marg.dens.beta[n_var, ]) + 0.6, paste("E(beta) a posteriori = ", round(Beta.posterior[n_var], digits = 4)), col = green_line)
    abline(v = Beta.posterior[n_var], col = green_line, lwd = 1)
    abline(h = line_level, col = red_line, lwd = 3)
    polygon(c(beta.space, rev(beta.space)), 
            c(posterior.cs, rep(0, length(beta.space))), 
            col = red_area, border = NA)
    
    text(Beta.posterior[n_var], max(posterior.marg.dens.beta[n_var, ]) - 3, paste(round(HPDI_probab * 10, digits = 0), "% przedzial HPDI: (", round(x1, digits = 2), " , ", round(x2, digits = 2), ")"), col = "black")
  },
  line_level = slider(0, max(posterior.marg.dens.beta[n_var, ]) + 0.002, step = 0.001, initial = max(posterior.marg.dens.beta[n_var, ]) + 0.001)
)



P_y_M1 <- ((det(U.posterior) ^ 0.5) * gamma(v.posterior / 2) * ((vs2.posterior) ^ (- v.posterior / 2))) / ((pi ^ (N.data / 2)) * (det(U.prior) ^ 0.5) * gamma(v.prior / 2) * ((vs2.prior) ^ (- v.prior / 2)))

P_y_M2 <- rep(NA, 3)
for (ii in 2:4) {
  X_2 <- X[, -c(ii)]
  eval(parse(text = paste("OLS_results_2 <- lm(Q_cuc ~ ", paste(colnames(X_2), collapse = "+"), ", data = data_log)", sep = "")))
  Beta.ols.data_2 <- OLS_results_2$coefficients
  v.data_2 <- OLS_results_2$df.residual
  XTX.data_2 <- t(X_2) %*% X_2
  s2.data_2 <- sum((OLS_results_2$residuals) ^ 2) / v.data_2
  
  Beta.prior_2 <- Beta.prior[-c(ii)]
  U.prior_2 <- U.prior[- c(ii), -c(ii)]
  sm2.prior_2 <- sm2.prior
  v.prior_2 <- v.prior
  vs2.prior_2 <- vs2.prior
  
  Beta.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2) %*% (solve(U.prior_2) %*% Beta.prior_2 + XTX.data_2 %*% Beta.ols.data_2)
  U.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2)
  v.posterior_2 <- v.prior_2 + N.data
  vs2.posterior_2 <- v.prior_2 / sm2.prior_2 + v.data_2 * s2.data_2 + t(Beta.ols.data_2 - Beta.prior_2) %*% solve(U.prior_2 + solve(XTX.data_2)) %*% (Beta.ols.data_2 - Beta.prior_2)
  sm2.posterior_2 <- 1 / (vs2.posterior_2 / v.posterior_2)
  
  P_y_M2[ii - 1] <- ((det(U.posterior_2) ^ 0.5) * gamma(v.posterior_2 / 2) * ((vs2.posterior_2) ^ (- v.posterior_2 / 2))) / ((pi ^ (N.data / 2)) * (det(U.prior_2) ^ 0.5) * gamma(v.prior_2 / 2) * ((vs2.prior_2) ^ (- v.prior_2 / 2)))
}

# 2. Wyznaczamy czynniki Bayesa "dla poszczególnych zmiennych"
BF_1_2 <- rep(P_y_M1, 3) / P_y_M2
BF_1_2_table <- data.frame(names(Beta.ols.data[2:4]), BF_1_2)
colnames(BF_1_2_table) <- c("zmienna", "czynnik Bayesa (analitycznie)")
(BF_1_2_table)

