require(spfda)
# generate signals
N = 4096

f1 = spfda.sample$sine_thred(0, 2*N, sd = 0)$X; # matplot(t(f1), type='l')
f2 = spfda.sample$sine(0, 2*N, sd = 0)$X; # matplot(t(f2), type='l')
# f3 = sine_var(32)$X[, round(seq(1,4000, length.out = 401))];  # matplot(t(f3), type='l')
f3 = spfda.sample$sine_var(N, sd = 0)$X;
sd = 1; ar = 1:4

f1 = f1[, 24 + (1:377)] * 10
f2 = f2[, 0 + (1:377)] * 3
f3 = f3[, 10 + seq(1, by = 10, length.out=377)] * 2


s1 = f1[1:N, ] + spfda.sample$gen_nois(377, N, ar = ar, sd = sd)
s2 = f2[1:N, ] + spfda.sample$gen_nois(377, N, ar = ar, sd = sd)
s3 = f1[N + seq_len(N), ] + f2[N + seq_len(N), ] + f3 + spfda.sample$gen_nois(377, N, ar = ar, sd = sd)

e = spfda.sample$gen_nois(377, N, ar = ar, sd = sd)
cov = cov(e)
eigen = eigen(cov)
W = (eigen$vectors) %*% diag((sqrt(eigen$values))^-1) %*% t(eigen$vectors)

X = rbind(s1,s2,s3)

Y = rep(c(1,2,3), each = N)

Y=t(sapply(Y,function(ii){list(c(1,0,0),c(0,1,0),c(1,1,1))[[ii]]}))
z=X;X=Y;Y=z


alpha = 0.5
NN = N * 3
L = 100

(log(L))^(alpha / (alpha-1))

lambda = 10 / (T * (L / sqrt(NN * T))^alpha); lambda

# lambda = (NN / L)^(alpha/2); lambda
re = spfda_admm(Y,X,lambda = lambda, alpha = alpha,max_iter = 10, nknots = L, use_admm = T, W=NULL)
matplot(t(as.matrix(re$eta %*% re$B)), type = 'n', lty=1)

matpoints(t(as.matrix(re$eta %*% re$B)), type = 'l', lty=1)
matpoints(cbind(colMeans(f1),colMeans(f2),colMeans(f3)), type='l', lty=2)


plot(re, std = 1)

re$comps[[3]]$recall

plot(re$W, type = 'l')

re$comps[[1]]$recall

plot(re$B %*% re$W[1:96], type='l')
points(colMeans(f1), type= 'l', col= 'blue')

plot(re$B %*% re$W[(1:96) + 96], type='l')
points(colMeans(f2), type= 'l', col= 'blue')
plot(re$B %*% re$W[(1:96) + 96*2], type='l')
points(colMeans(f3), type= 'l', col= 'blue')





