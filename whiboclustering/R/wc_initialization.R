#' General Component for Initialization of Cluster Representatives.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives.
#' @param initialization_type String which signal which initialization type to be used. Check \code{wc_init_types} for possible values.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_initialize <- function(data, k = 3, initialization_type)
{
  # CHECKING FOR ERRORS
  if (!(tolower(initialization_type) %in% tolower(wc_init_types$Type)))
  {
    stop('Please enter initialization function that is available in wc_init_types data frame')
  }

  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  #INITIALIZE CLUSTER REPRESENTATIVES
  centroids <- eval(call(name = as.character(wc_init_types$Method[tolower(wc_init_types$Type) == tolower(initialization_type)]), data, k))

  return(centroids)
}

#' Random Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_init_random <- function(data, k = 3)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  centroids <- as.data.frame(apply(X = data, MARGIN = 2, FUN = function(x) {stats::runif(n = k, min = min(x), max = max(x))}))
  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' Forgy algorithm Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_init_forgy <- function(data, k = 3)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  centroids <- as.data.frame(data[sample(x = 1:dim(data)[1], size = k, replace = FALSE), ])
  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' K-Means++ Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_init_kmeansplusplus <- function(data, k = 3)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  centroids <- numeric(k)
  prob <- rep(1, dim(data)[1])

  for(c in 1:k)
  {
    centroids[c] <- sample.int(n = dim(data)[1], size = 1, prob = prob)
    prob <- apply(X = data, MARGIN = 1, FUN = function(x) {sqrt(sum((x - data[centroids[c], ])^2))})
  }

  centroids <- as.data.frame(data[centroids, ])

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' KKZ Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_init_kkz <- function(data, k)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  centroids <- data[0, ]
  centroids <- rbind.data.frame(centroids, data[which.max(rowSums(apply(X = data, MARGIN = 2, FUN = function(x) sqrt(x^2/sum(x^2))))), ])

  for(c in 2:k)
  {
    el <- which.max(apply(X = data[!row.names(data) %in% row.names(centroids), ], MARGIN = 1, FUN = function(x) {sqrt(sum(x - centroids)^2)}))
    centroids <- rbind.data.frame(centroids, data[el, ])
  }
  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' PCA Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_init_pca <- function(data, k)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  pca_model <- stats::prcomp(data)

  centroids <- as.data.frame(t(apply(X = data, MARGIN = 2, FUN = mean) + pca_model$rotation[, 1:k] * apply(X = data, MARGIN = 2, FUN = stats::sd)))

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' DIANA Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom cluster diana
wc_init_diana <- function(data, k)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  #require(cluster)
  model <- cluster::diana(x = data, diss = FALSE, stop.at.k = k)

  cuts <- stats::cutree(tree = model, k = k)
  centroids <- stats::aggregate(data, by = list(cuts), FUN = mean)[, 2:(dim(data)[2] + 1)]

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' AGNES Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom cluster agnes
wc_init_agnes <- function(data, k)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  #require(cluster)
  model <- cluster::agnes(x = data, diss = FALSE, method = 'average')

  cuts <- stats::cutree(tree = model, k = k)
  centroids <- stats::aggregate(data, by = list(cuts), FUN = mean)[, 2:(dim(data)[2] + 1)]

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' Ward algorithm Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
#' @importFrom cluster agnes
wc_init_ward <- function(data, k)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  #require(cluster)
  model <- cluster::agnes(x = data, diss = FALSE, method = 'ward')

  cuts <- stats::cutree(tree = model, k = k)
  centroids <- stats::aggregate(data, by = list(cuts), FUN = mean)[, 2:(dim(data)[2] + 1)]

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' Quantile Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_init_quantile <- function(data, k)
{
  if(!(class(data) %in% c('data.frame', 'matrix')))
  {
    stop('Data should be data.frame or matrix')
  }

  if(!is.numeric(k))
  {
    stop('k should be numeric')
  }

  centroids <- as.data.frame(apply(X = data, MARGIN = 2, FUN = function(x) {stats::quantile(x = x, probs = (2 * 1:k - 1)/(2 * k))}))

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' CCIA Cluster Representatives initialization.
#'
#' @param data A dataset for which Cluster Representatives needs to be initialized.
#' @param k A number of Cluster Representatives to be initialized.
#' @return As a result initial Cluster Representatives are obtained. Result is in for of data.frame or data.matrix.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_init_ccia <- function(data, k)
{
  centroids <- as.data.frame(apply(X = data, MARGIN = 2, FUN = function(x) {mean(x, na.rm = TRUE) + stats::sd(x, na.rm = TRUE) * stats::qnorm((2 * 1:k - 1)/(2 * k))}))

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

#' Data frame for possible values of initialization types.
#' @author Sandro Radovanovic \email{sandro.radovanovic@@gmail.com}
wc_init_types <- data.frame()

wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'Random', 'Method' = 'wc_init_random'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'Forgy', 'Method' = 'wc_init_forgy'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'KMeans++', 'Method' = 'wc_init_kmeansplusplus'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'KKZ', 'Method' = 'wc_init_kkz'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'PCA', 'Method' = 'wc_init_pca'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'AGNES', 'Method' = 'wc_init_agnes'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'DIANA', 'Method' = 'wc_init_diana'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'Ward', 'Method' = 'wc_init_ward'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'Quantile', 'Method' = 'wc_init_quantile'))
wc_init_types <- rbind.data.frame(wc_init_types, data.frame('Type' = 'CCIA', 'Method' = 'wc_init_ccia'))
