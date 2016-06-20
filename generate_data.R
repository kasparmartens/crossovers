generate_discrete_HMM = function(n, starting_probs, hidden_transition, obs_transition){
  x = rep(NA, n)
  k_hidden = length(starting_probs)
  k_obs = ncol(obs_transition)
  x[1] = sample(1:k_hidden, 1)
  for(i in 2:n){
    x[i] = sample(1:k_hidden, 1, prob = hidden_transition[x[i-1], ])
  }
  y = rep(NA, n)
  for(i in 1:n){
    y[i] = sample(1:k_obs, 1, prob = obs_transition[x[i], ])
  }
  return(list(x = factor(x, levels = 1:k_hidden), y = factor(y, levels = 1:k_obs)))
}
