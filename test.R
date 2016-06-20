source("generate_data.R")
source("hmm.R")

# generate HMM with x in {1, 2, 3} and y in {1, 2, 3, 4}
alpha = c(1, 1, 1)
# starting_probs
pi = c(0.3, 0.5, 0.2)
# hidden_transition
A = rdirichlet(3, alpha)
# obs_transition
B = rdirichlet(length(pi), rep(1, 4))

n = 100
alpha0 = 1
hmm_obs = generate_discrete_HMM(n, pi, A, B)
y = hmm_obs$y

# starting values for the parameters pi, A, B
pi_params = rep(alpha0, 3)
A_params = matrix(alpha0, 3, 3)
B_params = matrix(alpha0, 3, 4)
pi = rdirichlet(1, pi_params)
A = rdirichlet_mat(A_params)
B = rdirichlet_mat(B_params)

res = gibbs_sampling_hmm(y, pi, A, B, max_iter = 1000)
# average number of matches with the true x (calculate for each mcmc iteration)
apply(res$trace_x, 1, function(x)mean(x == hmm_obs$x))
# marginal distribution
res$marginal_distr
