wc_eval_total_sum_of_squares <- function(data)
{
  return(sum(apply(X = data, 2, function(x) (x - mean(x, na.rm = TRUE))^2)))
}

wc_eval_within_sum_of_squares <- function(data, centroids, assignment)
{
  intra_cluster_distances <- c()


  for(c in 1:nrow(centroids))
  {
    intra_cluster_distances <- c(intra_cluster_distances, mean(apply(X = t(apply(X = data[assignment == c, ], MARGIN = 1, FUN = '-', as.matrix(centroids[c, !grepl('WCCluster', colnames(centroids))]))^2), MARGIN = 1, FUN = sum)))
  }
  return(intra_cluster_distances)
}

wc_eval_between_sum_of_squares <- function(data, centroids, assignment)
{
  return(wc_eval_total_sum_of_squares(data) - sum(wc_eval_within_sum_of_squares(data, centroids, assignment)))
}

wc_eval_ball_hall <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'ball_hall')[[1]])
}

wc_eval_banfeld_raftery <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'banfeld_raftery')[[1]])
}

wc_eval_c_index <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'c_index')[[1]])
}

wc_eval_calinski_harabasz <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'calinski_harabasz')[[1]])
}

wc_eval_davies_bouldin <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'davies_bouldin')[[1]])
}

wc_eval_det_ratio <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'det_ratio')[[1]])
}

wc_eval_dunn <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'dunn')[[1]])
}

wc_eval_gamma <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'gamma')[[1]])
}

wc_eval_g_plus <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'g_plus')[[1]])
}

wc_eval_silhouette <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'silhouette')[[1]])
}

wc_eval_xie_beni <- function(data, assignment)
{
  return(intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'xie_beni')[[1]])
}
