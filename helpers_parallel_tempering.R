
# evaluate mixture of MV gaussian densities at point x
gaussian_mixture_log_density = function(x, lambda, mu, Sigma){
  n_components = length(lambda)
  log_contribution = rep(NA, n_components)
  for(k in 1:n_components){
    log_contribution[k] = log(lambda[[k]]) + dmvnorm(x, mu[[k]], Sigma[[k]], log=TRUE)
  }
  log_density = log(sum(exp(log_contribution)))
  return(log_density)
}

calculate_density_grid = function(xval, yval, beta){
  values = outer(xval, yval)
  for(i in 1:length(xval)){
    for(j in 1:length(yval)){
      values[i, j] = beta * gaussian_mixture_log_density(c(xval[i], yval[j]), lambda, mu, Sigma)
    }
  }
  rownames(values) = xval
  colnames(values) = yval
  return(exp(values))
}


MH_step = function(x0, loglik, logprior, beta){
  proposal = x0 + rnorm(length(x0), 0, sqrt(1))
  if(runif(1) < exp(beta * loglik(proposal) + logprior(proposal) - beta * loglik(x0) - logprior(x0))){
    return(proposal)
  } else{
    return(x0)
  }
}

MH_step_between_chains = function(x0, loglik, temperatures, create_plots = TRUE){
  # pick an index j
  j = sample(1:length(temperatures), 1)
  # pick an adjacent one
  jj = ifelse(j == 1, j+1, ifelse(j == length(temperatures), j-1, sample(c(j-1, j+1), 1)))
  # propose to switch j and jj
  x1 = x0[[j]]
  x2 = x0[[jj]]
  beta1 = 1/temperatures[j]
  beta2 = 1/temperatures[jj]
  if(runif(1) < exp((beta1-beta2)*(loglik(x2)-loglik(x1)))){
    x0[[j]] = x2
    x0[[jj]] = x1
    if(create_plots) {
      df1 = data.frame(do.call("rbind", x0[c(jj, j)]), temperature = temperatures[c(j, jj)])
      print(p2 + geom_point(aes(X1, X2), data=df1, size=4, col="red"))
      df2 = data.frame(do.call("rbind", x0[c(j, jj)]), temperature = temperatures[c(j, jj)])
      df3 = data.frame(x = df1$X1, y = df1$X2, xend = df2$X1, yend=df2$X2, temperature = temperatures[c(j, jj)])
      print(p2 + geom_point(aes(X1, X2), data=df1, size=4, col="red") + 
              geom_point(aes(X1, X2), data=df2, size=4, col="red") + 
              geom_segment(aes(x, y, xend=xend, yend=yend), data=df3, col="red", linetype="dashed"))
      print(p2 + geom_point(aes(X1, X2), data=df2, size=4, col="red"))
    }
  }
  return(x0)
}
