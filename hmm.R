library(gtools)

initial_state_update = function(x){
  return(table(x[1]))
}

transition_mat_update = function(x){
  x1 = x[-length(x)]
  x2 = x[-1]
  return(table(x1, x2))
}

forward_backward = function(y, pi, A, B, k, n){
  # forward
  P = vector("list", n)
  temp = matrix(pi, k, k) * A * matrix(B[, y[2]], k, k, byrow=T)
  P[[2]] = temp / sum(temp)
  for(i in 3:n){
    temp = matrix(rowSums(P[[i-1]]), k, k) * A * matrix(B[, y[i]], k, k, byrow=T)
    P[[i]] = temp / sum(temp)
  }
  # backward 1
  x_draw = rep(NA, n)
  x_draw[n] = sample(1:k, 1, prob = colSums(P[[n]]))
  for(t in (n-1):1){
    x_draw[t] = sample(1:k, 1, prob = P[[t+1]][, x_draw[t+1]])
  }
  # backward 2
  Q = vector("list", n)
  Q[[n]] = P[[n]]
  for(t in (n-1):2){
    Q[[t]] = P[[t]] * matrix(colSums(Q[[t+1]]) / colSums(P[[t]]), k, k, byrow=T)
  }
  # most likely hidden state
  return(list(P = P, Q = Q, x_draw = factor(x_draw, levels = 1:k)))
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

gibbs_sampling_hmm = function(y, pi, A, B, alpha0 = 1, max_iter = 1000){
  n = length(y)
  k = nrow(A)
  marginal_distr = matrix(0, k, n)
  trace_x = matrix(0, max_iter, n)
  for(i in 1:max_iter){
    # update hidden states
    res = forward_backward(y, pi, A, B, k=k, n=n)
    # update transition matrices
    pi = as.numeric(rdirichlet(1, alpha0 + initial_state_update(res$x_draw)))
    A = rdirichlet_mat(alpha0 + transition_mat_update(res$x_draw))
    B = rdirichlet_mat(alpha0 + table(res$x_draw, y))
    marginal_distr = marginal_distr + 1/(max_iter)*compute_marginal_distribution(res$Q, k, n)
    trace_x[i, ] = res$x_draw
  }
  return(list(trace_x = trace_x, marginal_distr = marginal_distr))
}

