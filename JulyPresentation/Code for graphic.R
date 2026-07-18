library(latex2exp)

set.seed(22122022)

xval <- runif(50, 0, 10)

sq <- function(x){
	val <- (x - 5)^2
	return(val)
}

yval <- sq(xval) + rnorm(50, sd = 4) + 5

mod <- lm(yval ~ xval)

fit <- function(x, coeffs){
	y <- coeffs[2] * x + coeffs[1]
	return(y)
}

{
plot(xval, yval, xlab = TeX("$exp(f(t, \\Psi_i))$"),
		 ylab = TeX("$h(t, \\Psi_i)$"), pch = 16)
		 
lines(seq(0,10,0.01), sq(seq(0,10,0.01)) + 5, col = "blue", lwd = 2)

lines(seq(0,10,0.01), fit(seq(0,10,0.01), mod$coefficients),
			lwd = 2, col = "red")

text(3, 26,
		 labels = TeX("$h(t, \\Psi_i) = e^{f(t, \\Psi_i)^2}$"),
		 col = "blue", cex = 1.5)

text(5, 17,
		 labels = TeX("$h(t, \\Psi_i) = e^{\\beta \\times f(t, \\Psi_i)\\ +\\ \\gamma}$"),
		 col = "red", cex = 1.5)
}

summary(mod)
