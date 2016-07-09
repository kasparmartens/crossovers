library(gtools)

initial_state_update = function(x){
  return(table(x[1]))
}

transition_mat_update = function(x){
  x1 = x[-length(x)]
  x2 = x[-1]
  return(table(x1, x2))
}


forward_backward_rcpp = function(y, pi, A, B, k, n){
  # forward
  obj = forward_backward_fast(pi, A, B, y, k, n)
  return(list(P = obj$P, x_draw = factor(obj$x_draw, levels = 1:k)))
}

forward_backward = function(y, pi, A, B, k, n){
  # forward
  P = vector("list", n)
  temp = matrix(pi, k, k) * matrix(B[, y[1]], k, k, byrow=T)
  P[[1]] = temp / sum(temp)
  for(t in 2:n){
    temp = matrix(colSums(P[[t-1]]), k, k) * A * matrix(B[, y[t]], k, k, byrow=T)
    P[[t]] = temp / sum(temp)
  }
  # backward 1
  x_draw = rep(NA, n)
  x_draw[n] = sample(1:k, 1, prob = rowSums(P[[n]]))
  for(t in (n-1):1){
    x_draw[t] = sample(1:k, 1, prob = P[[t+1]][, x_draw[t+1]])
  }
  # backward 2
  # Q = vector("list", n)
  # Q[[n]] = P[[n]]
  # for(t in (n-1):2){
    # Q[[t]] = P[[t]] * matrix(colSums(Q[[t+1]]) / colSums(P[[t]]), k, k, byrow=T)
  # }
  # most likely hidden state
  return(list(P = P, x_draw = factor(x_draw, levels = 1:k)))
}

rdirichlet_mat = function(dirichlet_params){
  t(apply(dirichlet_params, 1, function(x)rdirichlet(1, x)))
}

compute_marginal_distribution = function(P_list, k, n){
  marginal_distr = matrix(0, k, n)
  for(t in length(P_list):2){
    marginal_distr[, t] = colSums(P_list[[t]])
  }
  marginal_distr[, 1] = rowSums(P_list[[2]])
  return(marginal_distr)
}

gibbs_sampling_hmm = function(y, n_hidden_states, alpha0 = 0.1, max_iter = 1000, burnin = 500){
  if(!is.factor(y)) stop("y must be a factor variable!")
  if(burnin >= max_iter) stop("burnin too large!")
  n = length(y)
  n_symbols = length(unique(y))
  
  # starting values for the parameters pi, A, B
  k = n_hidden_states
  pi = rep(1/k, k)
  A = 0.5 * diag(k) + 0.5*matrix(1/k, k, k)
  B = matrix(1/n_symbols, k, n_symbols)
  
  trace_x = matrix(0, max_iter-burnin, n)
  trace_A = list()
  trace_B = list()
  for(i in 1:max_iter){
    # update hidden states
    res = forward_backward_rcpp(y, pi, A, B, k=k, n=n)
    # update transition matrices
    pi = as.numeric(rdirichlet(1, 1 + initial_state_update(res$x_draw)))
    A = rdirichlet_mat(alpha0 + transition_mat_update(res$x_draw))
    B = rdirichlet_mat(alpha0 + table(res$x_draw, y))
    # marginal_distr = marginal_distr + 1/(max_iter)*compute_marginal_distribution(res$Q, k, n)
    if(i > burnin){
      trace_x[i-burnin, ] = res$x_draw
      trace_A[[i-burnin]] = A
      trace_B[[i-burnin]] = B
    }
    if(i %% 100 == 0) cat("iter", i, "\n")
  }
  return(list(trace_x = trace_x, trace_A = trace_A, trace_B = trace_B))
}

match_states = function(trace_x, trace_A, trace_B, true_B){
  for(i in 1:length(trace_B)){
    states = identify_states_KL(true_B, trace_B[[i]])
    # relabel matrix B
    trace_B[[i]][states, ] = trace_B[[i]]
    # relabel matrix A
    trace_A[[i]][states, states] = trace_A[[i]]
    # relabel hidden sequence x
    trace_x[i, ] = relabel_seq(trace_x[i, ], states)
  }
  return(list(trace_x = trace_x, trace_A = trace_A, trace_B = trace_B))
}

relabel_seq = function(x, states){
  out = x
  for(i in 1:length(unique(x))){
    out[x == i] = states[i]
  }
  return(out)
}

KL_distance = function(p, q) sum(ifelse(p == 0, 0, p * log(p / (q+1e-16))))

identify_states_KL = function(P, Q){
  n_states = nrow(Q)
  states = rep(NA, n_states)
  for(i in 1:n_states){
    distances = apply(P, 1, function(p)KL_distance(p, Q[i, ]))
    states[i] = which.min(distances)
  }
  # if some of the states coincide, change them randomly
  if(length(setdiff(1:n_states, states)) > 0){
    print(states)
    states[duplicated(states)] = sample(setdiff(1:n_states, states))
  }
  return(states)
}
