visualise_hidden_states = function(trace_x, true_x){
  n = length(true_x)
  df.m = data.frame(t(trace_x)) %>%
    mutate(t = 1:n) %>%
    melt(id.vars = "t")
  df_true = data.frame(t = 1:n, x = true_x)
  
  p = ggplot() + 
    geom_path(aes(t, factor(value), group=variable), df.m, col="grey50") + 
    # geom_path(aes(t, x, group=1), df_true, linetype="dashed") + 
    geom_point(aes(t, x, col=x), df_true) + 
    theme_bw()
  return(p)
}

summarise_transition_mat = function(trace_mat_list){
  n_row = nrow(trace_mat_list[[1]])
  n_col = ncol(trace_mat_list[[1]])
  trace_mat_list = lapply(trace_mat_list, function(x){
    rownames(x) <- 1:n_row
    colnames(x) <- 1:n_col
    return(x)
  })
  df.m = do.call("rbind", trace_mat_list) %>%
    melt() %>%
    group_by(Var1, Var2) %>%
    mutate(iter = 1:length(Var1))
  return(df.m)
}

visualise_transition_mat = function(trace_mat, true_mat){
  df = summarise_transition_mat(trace_mat)
  df2 = summarise_transition_mat(list(true_mat))
  
  p = ggplot(df, aes(iter, value, group=1)) + 
    geom_path() + 
    geom_hline(aes(yintercept=value), data=df2, col="red", linetype="dashed") +
    facet_grid(Var1 ~ Var2) + 
    theme_bw()
  return(p)
}

prop_transitions = function(x) mean(x[-length(x)] != x[-1])
correct_states = function(x, x_true) mean(x == x_true)

visualise_accuracy_and_switch_rate = function(trace_x, x_true){
  switch_rate = apply(trace_x, 1, prop_transitions)
  accuracy = apply(trace_x, 1, correct_states, x_true)
  df = data.frame(accuracy, switch_rate)
  p = ggplot(df, aes(switch_rate, accuracy)) + 
    geom_jitter() + theme_bw()
  return(p)
}
