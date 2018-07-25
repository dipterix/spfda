require(glmnet)
require(magrittr)
require(splines)
require(rave)
require(stringr)


# Generate X, Y
re = spfda.sample$sine(init_phase = 1, n = 200, ar = c(0), sd = 0.4)
X = re$X; dim(X)
Y = re$Y; length(Y)
time = re$time

i = 2
res = sparse_bs(X[Y==i, ], Y[Y==i], time, nknots = 150, niter = 10,  alpha = 0.01, lambda = 0.01)


## ECoG


# load subject
utils = rave_preprocess_tools()
utils$load_subject('YAB', 'Congruency1')
wave_info = utils$get_wavelet_log()
#utils$get_wavelet_log()$fre
epoch = load_meta('epoch', subject_id = utils$get_subject_id(), meta_name = 'YABa')
time_ind = seq(-100, 300)
blocks = utils$get_blocks()

load_EcoG = function(chl){

  re = lapply(blocks, function(block){
    # par(mfrow=c(1,2))
    re = utils$load_wave_data(block = block, chl = chl, complex = F)
    power = re$coef^2
    sapply(round(epoch$Time[epoch$Block == block] * wave_info$target_srate), function(ii){
      i = ii + time_ind
      x = power[6, i]
      x = x / median(x) - 1
    }) ->
      s
    sapply(str_split(epoch$Label[epoch$Block == block], '_'), function(x){
      if(length(x) == 2){
        x[[2]]
      }else{
        ''
      }
    }) ->
      label
    list(
      data = s,
      label = label
    )
  })

  X = t(do.call(cbind, lapply(re, function(x){x$data})))
  Y = unlist(lapply(re, function(x){x$label}))
  time = time_ind / wave_info$target_srate
  return(list(
    X=X,Y=Y,time=time
  ))
}
chls = c(13:16, 18:44, 47:63, 65:84)
alpha = 0.2
lambda = 0.2
niter = 2
re = load_EcoG(17)
time = re$time
lapply_async(chls, function(chl, ...){
  re = load_EcoG(chl)
  X = re$X; dim(X)
  Y = re$Y; length(Y)
  time = re$time
  sel_a = which(Y == 'a')
  sel_v = which(Y == 'v')
  sel_av = which(Y == 'av')
  n = min(length(sel_a), length(sel_v), length(sel_av))
  replicate(10, {
    X_delta = X[sample(sel_a, n), ] + X[sample(sel_v, n), ] - X[sample(sel_av, n), ]
    res = sparse_bs(X_delta, rep(1, n), time, nknots = 150, niter = 10,  alpha, lambda, plot = F, quiet = T)
    res$B %*% res$w
  }) -> recalls
  recalls = recalls[,1,]
  apply(recalls, 1, median)
}, .call_back = print) -> res

res = sapply(rave:::dropNulls(res), I)
rave:::image_plot(x = time, y = chls, z = res, symmetric = T)

unique(Y)
alpha = 0.2
lambda = 0.2
i = 'v'; res = sparse_bs(X[Y==i, ], Y[Y==i], time, nknots = 150, niter = 10,  alpha, lambda)
i = 'a'; res = sparse_bs(X[Y==i, ], Y[Y==i], time, nknots = 150, niter = 10,  alpha, lambda)
i = 'av'; res = sparse_bs(X[Y==i, ], Y[Y==i], time, nknots = 150, niter = 10,  alpha, lambda)

re = load_EcoG(38)
X = re$X; dim(X)
Y = re$Y; length(Y)
time = re$time
sel_a = which(Y == 'a')
sel_v = which(Y == 'v')
sel_av = which(Y == 'av')
n = min(length(sel_a), length(sel_v), length(sel_av))

X_delta = X[sample(sel_a, n), ] + X[sample(sel_v, n), ] - X[sample(sel_av, n), ]
res = sparse_bs(X, rep(1, n), time, nknots = 150, niter = 10,  alpha, lambda)
res$B %*% res$w
