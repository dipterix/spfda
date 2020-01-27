

require(spfda)
require(dipsaus)
require(splines)
require(fields)
require(refund)
require(ROCR)

set.seed(NULL)
seed = ceiling(runif(1) * 10000); seed
n = 100
n_timepoints = 100
ncores = 8
b = 0.003
K = 30
ord = 4

dipsaus::make_forked_clusters(workers = 8)

set_palette <- function(pal){
  pal %?<-% c("orange", "dodgerblue3", "darkgreen", "orangered", "brown",  "purple3", "black")
  grDevices::palette(pal)
}
set_palette()

# step 1 generate true data
env = simulate_data(n = n, n_timepoints = n_timepoints, seed = seed)
X = env$X
Y = env$Y
time = env$time
beta = env$beta
# step 2 get weight
W = get_weights(X, Y, b = b, parts = list(c(0, 0.4), c(0.4,0.8), c(0.8,1)))
# eigen(W)$val
# step 3 tune model
alphas = seq(0.05,0.95,0.05)
lambdas = 10^seq(-1, 1.5, length.out = 100)

a = get0('aa', ifnotfound = array(0, c(length(lambdas), length(alphas), 2)))

for(ii in rev(seq_along(alphas))){
  alpha = alphas[[ii]]
  re <- future.apply::future_sapply(lambdas, function(lambda){

    res = spfda(Y = Y, X = X, time = time, nsp = K, ord = ord, lambda = lambda,
                alpha = alpha, W = W, init = NULL, max_iter = 100, CI = FALSE)
    coef = res$get_coef()
    BIC = BIC(res)

    c(BIC, sqrt(mean((beta - coef) ^ 2)))

  }, future.seed = TRUE)

  print(alpha)
  a[,ii,] = t(re)
  par(mfrow = c(1,2))
  par(mar=c(2,2,2,2))
  for(ii in 1:2){
    try({
      fields::image.plot(as.matrix(a[,,ii]), x = lambdas, y = alphas)
    })
  }
}; aa = a


y = apply(a[,,1],2, which.min)
x = seq_along(alphas)

tbl = data.frame(
  alpha = alphas,
  lambda = lambdas[y],
  BIC = a[,,1][cbind(y,x)],
  RMSE = sqrt((a[,,2][cbind(y,x)])^2 / length(beta))
)
tbl
idx = which.min(tbl$BIC); idx
lambda = tbl$lambda[idx]
alpha = tbl$alpha[idx]
print(c(lambda, alpha))

ols_res = spfda(Y = Y, X = X, time = time, nsp = K, ord = 4, lambda = 0, alpha = 1, W = NULL)
res = spfda(Y = Y, X = X, time = time, nsp = K, ord = 4, lambda = lambda,
            alpha = alpha, W = W, init = ols_res$gamma, max_iter = 100, CI = TRUE)

plot_coef(res$get_coef(), res$raw$f_sd, beta = beta)

model_grbridge = list(coefficients=res$get_coef(), se = res$raw$f_sd)
model_mcp = run_grmcp(env$X, env$Y, K)
model_lasso = run_grlasso(env$X, env$Y, K)
model_for2s = run_for2s(env$X, env$Y, K)
model_ols = run_ols(env$X, env$Y, K)


time_partition = list(c(0,0.4), c(0.4,0.8), c(0.8,1), c(0,1))
models = list(
  'Proposed' = model_grbridge,
  'Group MCP' = model_mcp,
  'Group Lasso' = model_lasso,
  '2-Step FoS' = model_for2s,
  'OLS' = model_ols
)

# 1. Accuracy beta1, beta2, beta3, overall
RMSE = lapply(time_partition, function(part){
  sel = time >= part[1] & time < part[2]

  RMSE = sapply(models, function(m){
    re = sapply(1:3, function(j){
      sqrt(mean((m$coefficients[j, sel] - beta[j,sel])^2))
    })
    c(re, sqrt(mean((m$coefficients[, sel] - beta[,sel])^2)))
  })
  RMSE = as.data.frame(RMSE)
  colnames(RMSE) = names(models)

  RMSE$Time = c(sprintf('%.1f - %.1f', part[1], part[2]))
  RMSE$Beta = c(1:3, 'avg')
  RMSE
})
RMSE
# Compare 2 TPR beta_3 FPR beta_1 and 3

SPARSITY = lapply(time_partition, function(part){
  eps = 0

  sel = time >= part[1] & time < part[2]

  SPARSITY = sapply(models, function(m){

    coef = m$coefficients
    TP = abs(coef) > eps & (beta != 0)
    FN = abs(coef) <= eps & beta != 0
    FP = abs(coef) > eps & beta == 0
    TN = abs(coef) <= eps & beta == 0

    # FP / (FP + TN) for beta1
    c(
      FPR_1 = sum(FP[1,sel]) / ( sum(FP[1,sel]) + sum(TN[1,sel]) ),
      FPR_3 = sum(FP[3,sel]) / ( sum(FP[3,sel]) + sum(TN[3,sel]) ),
      TPR_2 = sum(TP[2,sel]) / ( sum(TP[2,sel]) + sum(FN[2,sel]) ),
      TPR_3 = sum(TP[3,sel]) / ( sum(TP[3,sel]) + sum(FN[3,sel]) )
    )
  })
  SPARSITY = as.data.frame(SPARSITY)
  colnames(SPARSITY) = names(models)

  SPARSITY$Time = c(sprintf('%.1f - %.1f', part[1], part[2]))
  SPARSITY
})
SPARSITY

# 3. coverage
COVERAGE = lapply(time_partition, function(part){
  sel = time >= part[1] & time < part[2]

  COVERAGE = sapply(models, function(m){
    re = sapply(1:3, function(j){
      mean(abs(m$coefficients[j, sel] - beta[j,sel]) < 2*m$se[j, sel])
    })
    c(re, mean(abs(m$coefficients[, sel] - beta[,sel]) < 2*m$se[, sel]))
  })
  COVERAGE = as.data.frame(COVERAGE)
  colnames(COVERAGE) = names(models)

  COVERAGE$Time = c(sprintf('%.1f - %.1f', part[1], part[2]))
  COVERAGE$Beta = c(1:3, 'avg')
  COVERAGE
})
COVERAGE


# save
saveRDS(list(
  seed = seed,
  X=X, Y=Y,
  lambda = lambda, alpha=alpha,
  models = models,
  tbl = tbl,
  RMSE = RMSE, SPARSITY = SPARSITY,
  COVERAGE = COVERAGE, K = K,
  beta = beta
), file = sprintf('%d.rds', seed))

