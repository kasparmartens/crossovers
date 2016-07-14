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
  x_draw[n] = sample(1:k, 1, prob = colSums(P[[n]]))
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

# likelihood p(y|x, A, B) = p(y|x, B)
loglikelihood = function(y, x, B){
  sum(log(B[cbind(x, y)]))
}

MH_step_between_chains = function(y, x, B, temperatures){
  # pick an index j, the other one will be j+1
  j = sample(1:(length(temperatures)-1), 1)
  x1 = x0[[j]]
  x2 = x0[[j+1]]
  beta1 = 1/temperatures[j]
  beta2 = 1/temperatures[j+1]
  if(runif(1) < exp((beta1-beta2)*(loglikelihood(y, x2, B)-loglikelihood(y, x1, B)))){
    x0[[j]] = x2
    x0[[jj]] = x1
  }
  return(x0)
}

parallel_tempering_hmm = function(y, n_hidden_states, temperatures, alpha0 = 0.1, max_iter = 1000, burnin = 500){
  if(!is.factor(y)) stop("y must be a factor variable!")
  if(burnin >= max_iter) stop("burnin too large!")
  n = length(y)
  n_symbols = length(unique(y))
  
  n_chains = length(temperatures)
  betas = 1/temperatures
  
  # starting values for the parameters pi, A, B
  k = n_hidden_states
  pi = lapply(1:n_chains, function(j)rep(1/k, k))
  A = lapply(1:n_chains, function(j)0.5 * diag(k) + 0.5*matrix(1/k, k, k))
  B = lapply(1:n_chains, function(j)matrix(1/n_symbols, k, n_symbols))
  
  trace_x = lapply(1:n_chains, function(x)matrix(0, max_iter-burnin, n))
  trace_A = lapply(1:n_chains, function(x)list())
  trace_B = lapply(1:n_chains, function(x)list())
  for(i in 1:max_iter){
    # update hidden states
    x = lapply(1:n_chains, function(j){
      forward_backward_rcpp(y, pi[[j]], A[[j]], B[[j]], k=k, n=n)$x_draw
    })
    # update transition matrices
    pi = lapply(1:n_chains, function(j){
      as.numeric(rdirichlet(1, 1 + initial_state_update(x[[j]])))
    })
    A = lapply(1:n_chains, function(j){
      rdirichlet_mat(alpha0 + transition_mat_update(x[[j]]))
    })
    B = lapply(1:n_chains, function(j){
      rdirichlet_mat(alpha0 + 1/temperatures[j]*table(x[[j]], y))
    })
    if(i > burnin){
      for(j in 1:n_chains){
        trace_x[[j]][i-burnin, ] = x[[j]]
        trace_A[[j]][[i-burnin]] = A[[j]]
        trace_B[[j]][[i-burnin]] = B[[j]]
      }
    }
    if(i %% 100 == 0) cat("iter", i, "\n")
  }
  return(list(trace_x = trace_x, trace_A = trace_A, trace_B = trace_B))
}
