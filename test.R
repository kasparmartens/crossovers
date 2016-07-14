library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)

source("generate_data.R")
source("hmm.R")
source("helpers_diagnostics.R")
source("helpers_relabelling.R")
Rcpp::sourceCpp("forward_backward.cpp")

# generate HMM with x in {1, 2, 3} and y in {1, 2, 3, 4}
pi = c(1/3, 1/3, 1/3)
A = rbind(c(0.9, 0, 0.1), c(0.1, 0.9, 0), c(0, 0.1, 0.9))
# B = rbind(c(0, 0.1, 0.45-0.1, 0.45+0.1), c(0, 0.1, 0.45, 0.45), c(0.6, 0.3, 0.1, 0))
B = rbind(c(0.5, 0, 0, 0.5), c(0, 0, 0.5, 0.5), c(0.6, 0.4, 0, 0))


n = 500
hmm_obs = generate_discrete_HMM(n, pi, A, B)
# hmm_obs = generate_nonmarkov_seq(n, B, n_breakpoints = 50)
df_true = data.frame(t = 1:n, x = hmm_obs$x)
ggplot(df_true, aes(t, x, group=1)) + 
  geom_path() + geom_point(aes(col = x)) + 
  theme_bw() + xlab("") + ggtitle("True hidden sequence")


res = gibbs_sampling_hmm(y = hmm_obs$y, n_hidden_states = nrow(A), alpha0 = 0.1, max_iter = 1500, burnin=500)
res_relabelled = match_states(res$trace_x, res$trace_A, res$trace_B, true_x = hmm_obs$x, k = nrow(A))

# visualise draws from the posterior of hidden states
print(visualise_hidden_states(res_relabelled$trace_x, hmm_obs$x, n_display = 100))

# posterior distribution of transition probabilities
print(visualise_transition_mat(res_relabelled$trace_A, A))
print(visualise_transition_mat(res_relabelled$trace_B, B))

# accuracy and switch rate
print(visualise_accuracy_and_switch_rate(res_relabelled$trace_x, hmm_obs$x))

# parallel tempering
temperatures = c(1, 4, 8, 16)
res = parallel_tempering_hmm(y = hmm_obs$y, n_hidden_states = nrow(A), temperatures = temperatures, 
                             alpha0 = 0.1, max_iter = 1500, burnin = 500)
res_relabelled = match_states_parallel_tempering(res$trace_x, res$trace_A, res$trace_B, true_x = hmm_obs$x, k = nrow(A))

print(visualise_hidden_states_parallel_tempering(res_relabelled$trace_x, hmm_obs$x, temperatures, n_display = 100))
print(visualise_transition_mat_parallel_tempering(res_relabelled$trace_A, A, temperatures))
print(visualise_transition_mat_parallel_tempering(res_relabelled$trace_B, B, temperatures))
