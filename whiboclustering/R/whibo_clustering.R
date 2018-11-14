#' As S4 class to represent WhiBo Cluster model
#'
#' @slot whibo_cluster Whibo Clustering object - list of objects for White-Box Clustering
#' @rdname whibo_cluster
#' @name whibo_cluster-class
#' @exportClass whibo_cluster
#' @author Sandro Radovanovic
methods::setClass(Class = 'whibo_cluster', representation = 'list')

source(file = 'R/wc_normalization.R')
source(file = 'R/wc_initialization.R')
source(file = 'R/wc_assignment.R')
source(file = 'R/wc_recalculate.R')
source(file = 'R/wc_cluster_performance.R')

#' Find Cluster model using White-Box Cluster Algorithm Design.
#'
#' @param data Data on which clustering should be performed.
#' @param k Number of Cluster Representatives.
#' @param normalization_type Which normalization should be used (look at \code{wc_norm_types} for possible values). Default value is \code{No}.
#' @param cluster_initialization_type Which initialization of Cluster Representatives should be used (look at \code{wc_init_types} for possible values). Default value is \code{Random}.
#' @param assignment_type Which assignment function should be used (look at \code{wc_assign_types} for possible values). Default value is \code{Euclidean}.
#' @param recalculation_type Which function for updating Cluster Representatives should be used (look at \code{wc_recalculate_types} for possible values). Default value is \code{Mean}.
#' @param max_iteration Number of iterations. Default value is 20.
#' @param no_of_restarts Number of restarts of whole clustering procedure. Default value is 1.
#' @return Object of type \code{whibo_cluster} which include Cluster Representatives (\code{centroids}), number of elements per cluster (\code{elements_per_cluster}), assignments (\code{assignments}), measures of cluster quality (\code{within_sum_of_squares}, \code{between_ss_div_total_ss} and \code{internal_measures_of_quality}), cluster models per iterations (\code{model_history}), iterations (\code{iterations}) and parameters used (\code{params})
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @seealso \code{plot.whibo_cluster}, \code{predict.whibo_cluster}
#' @importFrom clusterCrit intCriteria
#' @examples
#' data <- iris[, 1:4] #Take only numerical columns
#'
#' #Perform k-means clustering
#' model <- whibo_clustering(data = data, k = 3)
#' model
#'
#' #Perform some unorthodox clustering
#' model <- whibo_clustering(data = data, k = 3,
#'      normalization_type = 'Z', cluster_initialization_type = 'Ward',
#'      assignment_type = 'Correlation', recalculation_type = 'Trimean')
#'
#' @export
whibo_clustering <- function(data, k = 3,
                             normalization_type = 'No', cluster_initialization_type = 'Random',
                             assignment_type = 'Euclidean', recalculation_type = 'Mean',
                             max_iteration = 20, no_of_restarts = 1)
{
  #DATA NORMALIZATION
  model <- wc_normalize(data = data, normalization_type = normalization_type)
  data <- model$data

  #PARAMETERS
  params <- list('normalization_type' = normalization_type, 'cluster_initialization_type' = cluster_initialization_type,
                 'assignment_type' = assignment_type, 'recalculation_type' = recalculation_type, 'normalization_model' = model$model)

  history_of_cluster_models <- list()
  best_model <- NA
  best_performance <- 0

  #VARIABLES
  history_of_cluster_models <- list()
  best_model <- NA
  best_performance <- 0

  #TESTNG MULTIPLE TIMES OR ONCE IF SAID SO
  for(m in 1:no_of_restarts)
  {
    #HELP VARIABLES
    iteration <- 1
    assignment_hist <- list()
    centroids_hist <- list()

    #GENERATE INITIAL REPRESENTATIVES
    centroids <- wc_initialize(data = data, k = k,initialization_type = cluster_initialization_type)
    centroids_hist[[iteration]] <- centroids

    #ASSIGNING EXAMPLES TO CENTROIDS
    assignments <- wc_assignment(data = data, centroids = centroids, assignment_type = assignment_type)
    assignment_hist[[iteration]] <- assignments

    last_centroids <- centroids

    stoping <- FALSE
    while(!stoping)
    {
      iteration <- iteration + 1

      #RECALCULATE CLUSTER REPRESENTATIVES
      centroids <- wc_recalculate(data = data, assignment = assignments, recalculate_type = recalculation_type, assignment_type = assignment_type)
      centroids_hist[[iteration]] <- centroids

      #ASSIGNING EXAMPLES TO CENTROIDS
      assignments <- wc_assignment(data = data, centroids = centroids, assignment_type = assignment_type)
      assignment_hist[[iteration]] <- assignments

      #IF ASSIGNMENT WAS EMPTY
      for(cent in 1:k)
      {
        if(nrow(centroids[centroids$WCCluster == cent, ]) == 0)
        {
          centroids <- rbind.data.frame(centroids, last_centroids[last_centroids$WCCluster == cent, ])
        }
        centroids <- centroids[order(centroids$WCCluster, decreasing = FALSE), ]
      }

      #CENTROIDS STOPPING
      if (sum(centroids != last_centroids) == 0)
      {
        stoping <- TRUE
      }

      #ITERATION STOPPING
      if(iteration == max_iteration)
      {
        stoping <- TRUE
      }

      last_centroids <- centroids
    }
    #END WHILE
    elements_per_cluster <- table(assignments)

    within_ss <- wc_eval_within_sum_of_squares(data = data, centroids = last_centroids, assignment = assignments)
    between_ss <- wc_eval_between_sum_of_squares(data = data, centroids = last_centroids, assignment = assignments)
    total_ss <- wc_eval_total_sum_of_squares(data = data)

    evaluation <- clusterCrit::intCriteria(traj = as.matrix(data), part = as.integer(assignments), crit = 'all')

    output <- list('centroids' = centroids, 'elements_per_cluster' = elements_per_cluster, 'assignments' = assignments,
                   'within_sum_of_squares' = within_ss, 'between_ss_total_ss' = sum(between_ss)/total_ss,
                   'internal_measures_of_quality' = evaluation,
                   'history' = list('centroids_history' = centroids_hist, 'assignment_history' = assignment_hist, 'num_of_iterations' = iteration))

    history_of_cluster_models[[m]] <- output

    if(sum(between_ss)/total_ss > best_performance)
    {
      best_performance <- sum(between_ss)/total_ss
      best_model <- output
    }
  }
  #END FOR
  model_output <- list('centroids' = best_model$centroids, 'elements_per_cluster' = best_model$elements_per_cluster, 'assignments' = best_model$assignments,
                       'within_sum_of_squares' = best_model$within_sum_of_squares, 'between_ss_div_total_ss' = best_model$between_ss_total_ss,
                       'internal_measures_of_quality' = best_model$internal_measures_of_quality,
                       'model_history' = best_model$history, 'iterations' = history_of_cluster_models, 'params' = params)

  class(model_output) <- 'whibo_cluster'
  # model_output <- new(Class = 'whibo_cluster', model_output)

  return(model_output)
}

#' Show White-Box Cluster Algorithm model
#'
#' @param x WhiBo Cluster model.
#' @param ... None of those will be used.
#' @return Summary text about Cluster model.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @seealso \code{summary.whibo_cluster}
#' @examples
#' data <- iris[, 1:4] #Numerical data only
#'
#' model <- whibo_clustering(data = data, k = 3)
#' print(model)
#'
#' @rdname print
#' @export
print.whibo_cluster <- function(x, ...)
{
  model <- x

  cat('----------WhiBo cluster model----------')
  cat('\n')
  cat(sprintf('Normalization type:\t\t%s\n', model$params$normalization_type))
  cat(sprintf('Initialization type:\t\t%s\n', model$params$cluster_initialization_type))
  cat(sprintf('Assignment type:\t\t%s\n', model$params$assignment_type))
  cat(sprintf('Update repr. type:\t\t%s\n', model$params$recalculation_type))
  cat('---------------------------------------')

  cat('\n')
  cat('\n')

  cat(sprintf('Centroids:\n'))
  print(model$centroids[, !grepl('WCCluster', colnames(model$centroids))])

  cat('\n')
  cat('\n')

  cat(sprintf('Assignments:\n'))
  print(model$assignments)

  cat('\n')

  cat(sprintf('Number of elements per cluster:\n'))
  print(as.data.frame(t(as.matrix(model$elements_per_cluster))), row.names = F)

  cat('\n')
  cat('\n')
  cat(sprintf('Finished in %s iterations', model$model_history$num_of_iterations))

  cat('\n')
  cat('---------------------------------------')
  cat('\n')

  cat('within sum of squares per cluster\n')
  cat(model$within_sum_of_squares)
  cat('\n')
  cat('\n')
  cat(sprintf('Between SS / Total SS: \t%2.2f%%\n', round(model$between_ss_div_total_ss * 100, digits = 2)))
  cat(sprintf('Davies-Bouldin index: \t%1.3f\n', model$internal_measures_of_quality$davies_bouldin))
  cat(sprintf('Silhoutte index: \t%1.3f\n', model$internal_measures_of_quality$silhouette))
  cat(sprintf('Dunn index: \t\t%1.3f\n', model$internal_measures_of_quality$dunn))
  cat(sprintf('C index: \t\t%1.3f\n', model$internal_measures_of_quality$c_index))
  cat('Many more internal cluster evaluation metric can be found in internal_measures_of_quality list...')

  cat('\n\n')
  cat('You can access history of cluster model (each iteration) and each restart phase')
}

#' Show White-Box Cluster Algorithm model
#'
#' @param object WhiBo Cluster model.
#' @param ... None of those will be used.
#' @return Summary text about Cluster model.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @seealso \code{print.whibo_cluster}
#' @examples
#' data <- iris[, 1:4] #Numerical data only
#'
#' model <- whibo_clustering(data = data, k = 3)
#' summary(model)
#'
#' @rdname summary
#' @export
summary.whibo_cluster <- function(object, ...)
{
  model <- object

  cat('----------WhiBo cluster model----------')
  cat('\n')
  cat(sprintf('Normalization type:\t\t%s\n', model$params$normalization_type))
  cat(sprintf('Initialization type:\t\t%s\n', model$params$cluster_initialization_type))
  cat(sprintf('Assignment type:\t\t%s\n', model$params$assignment_type))
  cat(sprintf('Update repr. type:\t\t%s\n', model$params$recalculation_type))
  cat('---------------------------------------')

  cat('\n')
  cat('\n')

  cat(sprintf('Centroids:\n'))
  print(model$centroids[, !grepl('WCCluster', colnames(model$centroids))])

  cat('\n')
  cat('\n')

  cat(sprintf('Assignments:\n'))
  print(model$assignments)

  cat('\n')

  cat(sprintf('Number of elements per cluster:\n'))
  print(as.data.frame(t(as.matrix(model$elements_per_cluster))), row.names = F)

  cat('\n')
  cat('\n')
  cat(sprintf('Finished in %s iterations', model$model_history$num_of_iterations))

  cat('\n')
  cat('---------------------------------------')
  cat('\n')

  cat('within sum of squares per cluster\n')
  cat(model$within_sum_of_squares)
  cat('\n')
  cat('\n')
  cat(sprintf('Between SS / Total SS: \t%2.2f%%\n', round(model$between_ss_div_total_ss * 100, digits = 2)))
  cat(sprintf('Davies-Bouldin index: \t%1.3f\n', model$internal_measures_of_quality$davies_bouldin))
  cat(sprintf('Silhoutte index: \t%1.3f\n', model$internal_measures_of_quality$silhouette))
  cat(sprintf('Dunn index: \t\t%1.3f\n', model$internal_measures_of_quality$dunn))
  cat(sprintf('C index: \t\t%1.3f\n', model$internal_measures_of_quality$c_index))
  cat('Many more internal cluster evaluation metric can be found in internal_measures_of_quality list...')

  cat('\n\n')
  cat('You can access history of cluster model (each iteration) and each restart phase')
}

#' Plot WhiBo Cluster Representatives
#'
#' @param x WhiBo Cluster model.
#' @param ... None of those will be used.
#' @return Line plot with Cluster representatives
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @examples
#' data <- iris[, 1:4] #Numerical data only
#'
#' model <- whibo_clustering(data = data, k = 3)
#' plot(model)
#'
#' @rdname plot
#' @export
plot.whibo_cluster <- function(x, ...)
{
  model <- x

  #Create basic plot without axes
  graphics::matplot(t(model$centroids[, -1]), type = 'l', lty = rep(1, nrow(model$centroids)), xlab = 'Cluster Representatives', ylab = 'Value', axes = FALSE)

  #Add axes
  graphics::axis(2)
  graphics::axis(1 ,at = 1:ncol(model$centroids[, -1]),labels = colnames(model$centroids[, -1]))

  #Add box around plot
  graphics::box()
}

# registerS3method("plot","whibo_cluster","plot.whibo_cluster", envir = getNamespace("whiboclustering"))

#' Plot WhiBo Cluster Representatives
#'
#' @param model WhiBo Cluster model.
#' @param data Data used for clustering (optional).
#' @return Ploting pairs plot where Cluster representatives are presented with data (if provided).
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @examples
#' data <- iris[, 1:4] #Numerical data only
#'
#' model <- whibo_clustering(data = data, k = 3)
#' plot_pairs(model) #Ploting Cluster Representatives only
#'
#' plot_pairs(model, data) #Ploting Cluster Representatives and Data
#'
#' @rdname plot_pairs
#' @export
plot_pairs <- function(model, data)
{

  #Plot pairs
  if (missing(data))
  {
    graphics::plot(model$centroids[ , !grepl('WCCluster', colnames(model$centroids))], cex = 2, pch = 3, col = seq(1:nrow(model$centroids)))

    print('Data points are ommited from plot')
  }
  else
  {
    if(nrow(data) != length(model$assignments))
    {
      stop('There is discrepency between data and model')
    }

    new_data <- eval(call(name = as.character(wc_norm_types$Method[tolower(wc_norm_types$Type) == tolower(model$params$normalization_type)]), data, model$params$normalization_model))
    new_data <- rbind.data.frame(new_data$data, model$centroids[, !grepl('WCCluster', colnames(model$centroids))])

    graphics::plot(new_data, cex = c(rep(0.6, nrow(new_data) - nrow(model$centroids)), rep(1.4, nrow(model$centroids))), pch = c(rep(1, nrow(new_data) - nrow(model$centroids)), rep(3, nrow(model$centroids))), col = c(model$assignments, seq(1:nrow(model$centroids))))
  }
}

#' Predict to which Cluster new data belongs
#'
#' @param object WhiBo Cluster model.
#' @param data Data for which Cluster should be obtained.
#' @param ... None of those will be used.
#' @return Vector of assignments.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @examples
#' data <- iris[1:100, 1:4] #Numerical data only and first 100 rows
#'
#' model <- whibo_clustering(data = data, k = 3)
#' predict(object = model, data = iris[101:150, 1:4])
#'
#' @rdname predict
#' @export
predict.whibo_cluster <- function(object, data, ...)
{
  model <- object
  if(missing(data))
  {
    return(model$assignments)
  }
  else
  {
    new_data <- eval(call(name = as.character(wc_norm_types$Method[tolower(wc_norm_types$Type) == tolower(model$params$normalization_type)]), data, model$params$normalization_model))
    return(wc_assignment(data = new_data$data, centroids = model$centroids, assignment_type = model$params$assignment_type))
  }
}

#Generate Manual file - Commented, but not forgoten
#system("R CMD Rd2pdf . --title=WhiBoClustering yourpackagename --output=././manual.pdf --force --no-clean --internals")
