library(combinat)

match_states = function(trace_x, trace_A, trace_B, true_x, k){
  perms = permn(1:k)
  for(i in 1:length(trace_B)){
    relabelled_x = lapply(perms, function(current_permutation){
      suppressMessages(plyr::mapvalues(trace_x[i, ], 1:k, current_permutation))
    })
    pointwise_acc = sapply(relabelled_x, function(x) mean(true_x == x))
    states = perms[[which.max(pointwise_acc)]]
    
    # relabel matrix B
    trace_B[[i]][states, ] = trace_B[[i]]
    # relabel matrix A
    trace_A[[i]][states, states] = trace_A[[i]]
  }
  return(list(trace_x = trace_x, trace_A = trace_A, trace_B = trace_B))
}

# match_states = function(trace_x, trace_A, trace_B, true_B){
#   for(i in 1:length(trace_B)){
#     states = identify_states_KL(true_B, trace_B[[i]])
#     # relabel matrix B
#     trace_B[[i]][states, ] = trace_B[[i]]
#     # relabel matrix A
#     trace_A[[i]][states, states] = trace_A[[i]]
#     # relabel hidden sequence x
#     trace_x[i, ] = relabel_seq(trace_x[i, ], states)
#   }
#   return(list(trace_x = trace_x, trace_A = trace_A, trace_B = trace_B))
# }

# relabel_seq = function(x, states){
#   out = x
#   for(i in 1:length(unique(x))){
#     out[x == i] = states[i]
#   }
#   return(out)
# }
# 
# KL_distance = function(p, q) sum(ifelse(p == 0, 0, p * log(p / (q+1e-16))))
# 
# sample_fix <- function(x, ...) x[sample(length(x), ...)]
# 
# identify_states_KL = function(P, Q){
#   n_states = nrow(Q)
#   states = rep(NA, n_states)
#   for(i in 1:n_states){
#     distances = apply(P, 1, function(p)KL_distance(p, Q[i, ]))
#     states[i] = which.min(distances)
#   }
#   # if some of the states coincide, change them randomly
#   if(length(setdiff(1:n_states, states)) > 0){
#     states[duplicated(states)] = sample_fix(setdiff(1:n_states, states))
#   }
#   return(states)
# }
