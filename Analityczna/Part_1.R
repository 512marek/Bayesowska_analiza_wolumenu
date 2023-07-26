rm(list = ls())

data <- read.csv("C:/Users/Marek/Documents/R_Scripts/Ekonometria_bayesowska/ogorki.csv")

data <- data[, c("Q_cuc", "P_cuc","P_tom", "P_oni")]
names(data) <- c("Wolumen ogorki", "Cena ogorki", "Cena pomidory", "Cena cebula")
data_log <- log(data)

model <- lm(`Wolumen ogorki` ~ ., data = data_log)
summary(model)

if(!require("manipulate")) {install.packages("manipulate"); library(manipulate)}
if(!require("mvtnorm")) {install.packages("mvtnorm"); library(mvtnorm)}

y <- as.matrix(data$`Wolumen ogorki`)
N.data <- length(y)
X <- cbind(as.matrix(rep(1, N.data)), 
           as.matrix(data[, c("Cena ogorki", "Cena pomidory", "Cena cebula")]))
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


grey_area <- rgb(160, 160, 160, 80, names = NULL, maxColorValue = 255)
grey_line <- rgb(80, 80, 80, 160, names = NULL, maxColorValue = 255)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)

par(mfrow = c(1, 1))
for(jj in 1:length(Beta.posterior)) {
  
  if (jj == 1){
  beta.space <- seq(from = 8, to = 18, by = 0.001)
}
  
  if (jj == 2){
    beta.space <- seq(from = -1.3, to = 0, by = 0.001)
}
  
  if (jj == 3){
    beta.space <- seq(from = -1, to = 0.2, by = 0.001)
}
  
  if (jj == 4){
    beta.space <- seq(from = -0.8, to = 0.2, by = 0.001)
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
  
  title <- ifelse(jj==1, "Stala", colnames(data)[jj])
  plot(beta.space, prior.marg.dens.beta[jj, ], las = 1, lwd = 2, bty = "n", col = grey_area,
       ylim = c(0, max(c(max(prior.marg.dens.beta[jj, ]),max(posterior.marg.dens.beta[jj, ]))) + 1), type = "l", ylab = "gęstość", main = title)
  polygon(c(beta.space, rev(beta.space)), c(prior.marg.dens.beta[jj, ], 
                                            rep(0, length(beta.space))), col = grey_area, border = NA)
  abline(v = Beta.prior[jj], col = grey_line, lwd = 1)
  text(Beta.prior[jj], max(posterior.marg.dens.beta[jj, ]) * 0.4, paste("E(B) a priori = ", Beta.prior[jj]), col = grey_line)
  abline(v = Beta.ols.data[jj], col = rgb(0, 0, 0, 1), lwd = 1)
  text(Beta.posterior[jj], max(posterior.marg.dens.beta[jj, ]) * 0.8, paste("OLS = ", round(Beta.ols.data[jj], 4)), col = rgb(0, 0, 0, 1))
  lines(beta.space, posterior.marg.dens.beta[jj, ], lwd = 2, col = green_line)
  polygon(c(beta.space, rev(beta.space)), c(posterior.marg.dens.beta[jj, ], 
                                            rep(0, length(beta.space))), col = green_area, border = NA)
  abline(v = Beta.posterior[jj], col = green_line, lwd = 1)
  text(Beta.posterior[jj], max(posterior.marg.dens.beta[jj, ]) + 0.6, paste("E(B) a posteriori = ", round(Beta.posterior[jj], digits = 4)), col = green_line)
}

