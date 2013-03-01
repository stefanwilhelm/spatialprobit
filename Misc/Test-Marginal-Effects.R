library(spatialprobit)

# Matlab
#% Marginal Effects for n=400/m=1    beta = [ 0;  1;  -1 ];  % (3 x 1)
#%Direct           lower 01      Coefficient         upper 99         mean/std           t-prob
#%x1               0.156380         0.200033         0.233789        13.647970         0.000000
#%x2              -0.247969        -0.211052        -0.168369       -13.983016         0.000000
#
#%Indirect         lower 01      Coefficient         upper 99         mean/std           t-prob
#%x1               0.182758         0.339019         0.557173         4.758287         0.000003
#%x2              -0.591885        -0.358437        -0.196150        -4.590112         0.000006
#
#%Total            lower 01      Coefficient         upper 99         mean/std           t-prob
#%x1               0.361572         0.539052         0.736049         7.204819         0.000000
#%x2              -0.827782        -0.569489        -0.383547        -6.850810         0.000000
#

set.seed(2)
fit <- LeSagePaceExperiment(n=400, beta = c(0,  1, -1), ndraw=1000, burn.in = 200, m=1,
 computeMarginalEffects = TRUE)
summary(fit)
impacts(fit)

set.seed(2)
n <- 400
beta <- c(0,  1, -1)
rho <- 0.75
X <- cbind(intercept = 1, x = rnorm(n), y = rnorm(n))
I_n <- sparseMatrix(i = 1:n, j = 1:n, x = 1)
W <- kNearestNeighbors(x = rnorm(n), y = rnorm(n), k = 6)
eps <- rnorm(n = n, mean = 0, sd = 1)
z <- solve(qr(I_n - rho * W), X %*% beta + eps)
y <- as.double(z >= 0)
S <- (I_n - rho * W)

# Theoretische Impacts
# S_r(W) <- phi(S^(-1) * x_r * beta_r) %*% S^(-1) * I_n * beta_r
mu <- as.double(solve(qr(S), X %*% beta))
S_r <- diag(dnorm(mu)) %*% solve(S)
# Average Direct Impact =  n^(-1) tr(S_r(W)) beta_r
direct <- mean(diag(S_r)) * beta[-1]   # 0.2086822 -0.2086822
# Average Total Impact =  n^(-1) 1' S_r(W) 1 beta_r
total <- mean(rowSums(S_r)) * beta[-1] # 0.7168747 -0.7168747

# Average Indirect Impacts =  Average Total Impact - Average Direct Impact
indirect <- total - direct
cbind(direct, indirect, total)

#> cbind(direct, indirect, total)
#         direct   indirect      total
#[1,]  0.2086822  0.5081925  0.7168747
#[2,] -0.2086822 -0.5081925 -0.7168747





direct <- cbind(
  lower_01 = apply(fit$direct, 2, quantile, prob=0.01),
  Coefficient = apply(fit$direct, 2, mean),
  upper_01 = apply(fit$direct, 2, quantile, prob=0.99),
  mean_sd = apply(fit$direct, 2, mean)/apply(fit$direct, 2, sd)
)
direct


indirect <- cbind(
  lower_01 = apply(fit$indirect, 2, quantile, prob=0.01),
  Coefficient = apply(fit$indirect, 2, mean),
  upper_01 = apply(fit$indirect, 2, quantile, prob=0.99),
  mean_sd = apply(fit$indirect, 2, mean)/apply(fit$indirect, 2, sd)
)
indirect

total <- cbind(
  lower_01 = apply(fit$total, 2, quantile, prob=0.01),
  Coefficient = apply(fit$total, 2, mean),
  upper_01 = apply(fit$total, 2, quantile, prob=0.99),
  mean_sd = apply(fit$total, 2, mean)/apply(fit$total, 2, sd)
)
total


mfx <- marginal.effects(fit)

direct = cbind(
  lower_01 = apply(mfx$direct, 2, quantile, prob=0.01),
  Coefficient = apply(mfx$direct, 2, mean),
  upper_01 = apply(mfx$direct, 2, quantile, prob=0.99),
  mean_sd = apply(mfx$direct, 2, mean)/apply(mfx$direct, 2, sd)
)
direct

# Direct
#    lower_01 Coefficient  upper_01   mean_sd
#x  0.1685635   0.2134651  0.247906 12.590938
#y -0.2191582  -0.1783309 -0.133471 -9.764394

indirect = cbind(
  lower_01 = apply(mfx$indirect, 2, quantile, prob=0.01),
  Coefficient = apply(mfx$indirect, 2, mean),
  upper_01 = apply(mfx$indirect, 2, quantile, prob=0.99),
  mean_sd = apply(mfx$indirect, 2, mean)/apply(mfx$indirect, 2, sd)
)
indirect
# --> indirekte Effekte sehen fast wie totale Effekte in Matlab aus

#    lower_01 Coefficient   upper_01   mean_sd
#x  0.3110425   0.5205660  0.8325988  4.769106
#y -0.7898519  -0.5096722 -0.3059449 -4.644033

total = cbind(
  lower_01 = apply(mfx$total, 2, quantile, prob=0.01),
  Coefficient = apply(mfx$total, 2, mean),
  upper_01 = apply(mfx$total, 2, quantile, prob=0.99),
  mean_sd = apply(mfx$total, 2, mean)/apply(mfx$total, 2, sd)
);

#    lower_01 Coefficient   upper_01   mean_sd
#x  0.4942853   0.7212815  1.0256794  6.293185
#y -1.0006702  -0.7062377 -0.4900613 -6.020970