library(reshape2)
library(ggplot2)
library(dplyr)

source("generate_data.R")
source("hmm.R")
source("helpers_diagnostics.R")
Rcpp::sourceCpp("forward_backward.cpp")

# generate HMM with x in {1, 2, 3} and y in {1, 2, 3, 4}
pi = c(1/3, 1/3, 1/3)
A = rbind(c(0.9, 0, 0.1), c(0.1, 0.9, 0), c(0, 0.1, 0.9))
# B_true = rbind(c(0, 0.1, 0.45-0.1, 0.45+0.1), c(0, 0.1, 0.45, 0.45), c(0.6, 0.3, 0.1, 0))
B = rbind(c(0.5, 0, 0, 0.5), c(0, 0, 0.5, 0.5), c(0.6, 0.4, 0, 0))


n = 500
hmm_obs = generate_discrete_HMM(n, pi, A, B)
y = hmm_obs$y
df_true = data.frame(t = 1:n, x = hmm_obs$x)
ggplot(df_true, aes(t, x, group=1)) + 
  geom_path() + geom_point(aes(col = x)) + 
  theme_bw()


res = gibbs_sampling_hmm(y, n_hidden_states = nrow(A), alpha0 = 0.1, max_iter = 1500)
res_relabelled = match_states(res$trace_x, res$trace_A, res$trace_B, B)

# visualise draws from the posterior of hidden states
print(visualise_hidden_states(res_relabelled$trace_x, hmm_obs$x))

# posterior distribution of transition probabilities
print(visualise_transition_mat(res_relabelled$trace_A, A))
print(visualise_transition_mat(res_relabelled$trace_B, B))

# accuracy and switch rate
print(visualise_accuracy_and_switch_rate(res_relabelled$trace_x, hmm_obs$x))

