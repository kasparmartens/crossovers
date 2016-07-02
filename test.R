library(reshape)
library(ggplot2)

source("generate_data.R")
source("hmm.R")
Rcpp::sourceCpp("forward_backward.cpp")

# generate HMM with x in {1, 2, 3} and y in {1, 2, 3, 4}
pi = c(1/3, 1/3, 1/3)
A = rbind(c(0.95, 0, 0.05), c(0.05, 0.95, 0), c(0, 0.05, 0.95))
# B_true = rbind(c(0, 0.1, 0.45-0.1, 0.45+0.1), c(0, 0.1, 0.45, 0.45), c(0.6, 0.3, 0.1, 0))
B = rbind(c(0.5, 0, 0, 0.5), c(0, 0, 0.5, 0.5), c(0.6, 0.4, 0, 0))


n = 2000
hmm_obs = generate_discrete_HMM(n, pi, A, B)
y = hmm_obs$y
df = data.frame(t = 1:n, x = hmm_obs$x, y = hmm_obs$y)
df.m = melt(df, id.vars = "t")
ggplot(df.m, aes(t, value, group=1)) + 
  geom_path() + geom_point(aes(col=value)) + 
  facet_wrap(~ variable, nrow=2) + 
  theme_bw()


res = gibbs_sampling_hmm(y, n_hidden_states = nrow(A), alpha0 = 0.1, max_iter = 1000)
lapply(res$trace_A[950:1000], round, 2)
lapply(res$trace_B[950:1000], round, 2)


latent.m = melt(res$trace_x[900:1000, ])
ggplot(latent.m, aes(X2, value)) + geom_path() + facet_wrap(~ X1)

# average number of matches with the true x (calculate for each mcmc iteration)
apply(res$trace_x, 1, function(x)mean(x == hmm_obs$x))

