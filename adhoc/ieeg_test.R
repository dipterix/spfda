require(spfda)
require(stringr)
require(rave)

m = unlist(load_modules())[[3]]
init_app(m)
rm(list = local({l = ls(all.names = T); l[l!= 'm']}));
debug_module(m)

d = t(collapse$over_frequency(bl_power))

cond = module_tools$get_meta('trials')$Condition
cond

av = str_detect(cond, '_av$')
a = str_detect(cond, '_a$')
v = str_detect(cond, '_v$')

sel = av | a | v
Y = d[sel, ]

X = cbind(a,v,av)[sel,]
X[X[,3], 1:2] = T


alpha = 0.4
NN = nrow(Y)
L = 100
(log(L))^(alpha / (alpha-1))

lambda = 1 / (T * (L / sqrt(NN * T))^alpha); lambda
re = spfda_admm(Y / 100,X,lambda = lambda, alpha = alpha,max_iter = 10, nknots = L, use_admm = T, W=NULL)
matplot(t(as.matrix(re$eta %*% re$B)) * 100, type = 'l', lty=1)

electrode
