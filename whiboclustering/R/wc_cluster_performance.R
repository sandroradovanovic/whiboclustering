#' Calculate total sum of squares
#'
#' @param data A dataset for which total sum of squared should be calculated.
#' @return A number which shows total sum of squares.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_eval_total_sum_of_squares <- function(data)
{
  return(sum(apply(X = data, 2, function(x) (x - mean(x, na.rm = TRUE))^2)))
}

#' Calculate within (Cluster) sum of squares
#'
#' @param data A dataset for which within sum of squared should be calculated.
#' @param centroids A data frame of cluster representatives.
#' @param assignment Vector of assignments.
#' @return A vector of number which shows within (cluster) sum of squares.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_eval_within_sum_of_squares <- function(data, centroids, assignment)
{
  intra_cluster_distances <- c()


  for(c in 1:nrow(centroids))
  {
    if (sum(assignment == c) == 0)
    {
      intra_cluster_distances <- c(intra_cluster_distances, 0)
    }
    intra_cluster_distances <- c(intra_cluster_distances, mean(apply(X = t(apply(X = data[assignment == c, ], MARGIN = 1, FUN = '-', as.matrix(centroids[c, !grepl('WCCluster', colnames(centroids))]))^2), MARGIN = 1, FUN = sum)))
  }
  return(intra_cluster_distances)
}

#' Calculate between (Clusters) sum of squares
#'
#' @param data A dataset for which between sum of squared should be calculated.
#' @param centroids A data frame of cluster representatives.
#' @param assignment Vector of assignments.
#' @return A vector of number which shows between (clusters) sum of squares.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_eval_between_sum_of_squares <- function(data, centroids, assignment)
{
  return(wc_eval_total_sum_of_squares(data) - sum(wc_eval_within_sum_of_squares(data, centroids, assignment)))
}

#' Calculate Ball-Hall internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_ball_hall <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'ball_hall')[[1]])
}

#' Calculate Banfeld-Raftery internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_banfeld_raftery <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'banfeld_raftery')[[1]])
}

#' Calculate C index internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_c_index <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'c_index')[[1]])
}

#' Calculate Calinski-Harabasz internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_calinski_harabasz <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'calinski_harabasz')[[1]])
}

#' Calculate Davies-Bouldin internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_davies_bouldin <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'davies_bouldin')[[1]])
}

#' Calculate Det ratio internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_det_ratio <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'det_ratio')[[1]])
}

#' Calculate Dunn index internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_dunn <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'dunn')[[1]])
}

#' Calculate Gamma index internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_gamma <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'gamma')[[1]])
}

#' Calculate G+ index internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_g_plus <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'g_plus')[[1]])
}

#' Calculate Silhouette score internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_silhouette <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'silhouette')[[1]])
}

#' Calculate Xie-Beni internal Cluster evaluation measure
#'
#' @param data A dataset for which internal cluster quality should be calculated.
#' @param assignment Vector of assignments.
#' @return A value of internal cluster quality evaluation measure.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom clusterCrit intCriteria
wc_eval_xie_beni <- function(data, assignment)
{
  return(clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignment), crit = 'xie_beni')[[1]])
}
