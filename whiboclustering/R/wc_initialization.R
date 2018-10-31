wc_initialize <- function(data, k = 3, initialization_type)
{
  # CHECKING FOR ERRORS
  if (!(tolower(initialization_type) %in% tolower(wc_init_types$Type)))
  {
    stop('Please enter initialization function that is available in wc_init_types data frame')
  }

  #INITIALIZE CLUSTER REPRESENTATIVES
  centroids <- eval(call(name = as.character(wc_init_types$Method[tolower(wc_init_types$Type) == tolower(initialization_type)]), data, k))

  return(centroids)
}

wc_init_random <- function(data, k = 3)
{
  centroids <- as.data.frame(apply(X = data, MARGIN = 2, FUN = function(x) {runif(n = k, min = min(x), max = max(x))}))
  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

wc_init_forgy <- function(data, k = 3)
{
  centroids <- data[sample(x = 1:dim(data)[1], size = k, replace = FALSE), ]
  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

wc_init_kmeansplusplus <- function(data, k = 3)
{
  centroids <- numeric(k)
  prob <- rep(1, dim(data)[1])

  for(c in 1:k)
  {
    centroids[c] <- sample.int(n = dim(data)[1], size = 1, prob = prob)
    prob <- apply(X = data, MARGIN = 1, FUN = function(x) {sqrt(sum((x - data[centroids[c], ])^2))})
  }

  centroids <- data[centroids, ]

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

wc_init_kkz <- function(data, k)
{
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

wc_init_pca <- function(data, k)
{
  pca_model <- prcomp(data)

  centroids <- t(apply(X = data, MARGIN = 2, FUN = mean) + pca_model$rotation[, 1:k] * apply(X = data, MARGIN = 2, FUN = sd))

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

wc_init_diana <- function(data, k)
{
  require(cluster)
  model <- diana(x = data, diss = FALSE, stop.at.k = k)

  cuts <- cutree(tree = model, k = k)
  centroids <- aggregate(data, by = list(cuts), FUN = mean)[, 2:(dim(data)[2] + 1)]

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

wc_init_agnes <- function(data, k)
{
  require(cluster)
  model <- agnes(x = data, diss = FALSE, method = 'average')

  cuts <- cutree(tree = model, k = k)
  centroids <- aggregate(data, by = list(cuts), FUN = mean)[, 2:(dim(data)[2] + 1)]

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

wc_init_ward <- function(data, k)
{
  require(cluster)
  model <- agnes(x = data, diss = FALSE, method = 'ward')

  cuts <- cutree(tree = model, k = k)
  centroids <- aggregate(data, by = list(cuts), FUN = mean)[, 2:(dim(data)[2] + 1)]

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

wc_init_quantile <- function(data, k)
{
  centroids <- apply(X = data, MARGIN = 2, FUN = function(x) {quantile(x = x, probs = (2 * 1:k - 1)/(2 * k))})

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

wc_init_ccia <- function(data, k)
{
  centroids <- apply(X = data, MARGIN = 2, FUN = function(x) {mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE) * qnorm((2 * 1:k - 1)/(2 * k))})

  centroids$WCCluster <- 1:k
  row.names(centroids) <- 1:k
  centroids <- centroids[c('WCCluster', setdiff(names(centroids), 'WCCluster'))]
  return(centroids)
}

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
