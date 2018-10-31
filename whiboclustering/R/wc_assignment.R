wc_assignment <- function(data, centroids, assignment_type)
{
  # CHECKING FOR ERRORS
  if (!(tolower(assignment_type) %in% tolower(wc_assign_types$Type)))
  {
    stop('Please enter assignment function that is available in wc_assign_types data frame')
  }

  #ASSIGN EXAMPLES TO CENTROIDS
  assignment <- eval(call(name = as.character(wc_assign_types$Method[tolower(wc_assign_types$Type) == tolower(assignment_type)]), data, centroids))

  return(assignment)
}


wc_assign_euclidean <- function(data, centroids)
{
  assignment <- apply(X = as.matrix(dist(rbind(data, centroids[, !grepl("WCCluster", colnames(centroids))]), method = 'euclidean'))[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(centroids))], MARGIN = 1, FUN = function(x) {which.min(x)})
  assignment <- as.numeric(assignment)
  return(assignment)
}

wc_assign_squared_euclidean <- function(data, centroids)
{
  assignment <- apply(X = as.matrix(dist(rbind(data^2, centroids[, !grepl("WCCluster", colnames(centroids))]^2), method = 'euclidean'))[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(centroids))], MARGIN = 1, FUN = function(x) {which.min(x)})
  assignment <- as.numeric(assignment)
  return(assignment)
}

wc_assign_manhattan <- function(data, centroids)
{
  assignment <- apply(X = as.matrix(dist(rbind(data, centroids[, !grepl("WCCluster", colnames(centroids))]), method = 'manhattan'))[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(centroids))], MARGIN = 1, FUN = function(x) {which.min(x)})
  assignment <- as.numeric(assignment)
  return(assignment)
}

wc_assign_canberra <- function(data, centroids)
{
  assignment <- apply(X = as.matrix(dist(rbind(data, centroids[, !grepl("WCCluster", colnames(centroids))]), method = 'canberra'))[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(centroids))], MARGIN = 1, FUN = function(x) {which.min(x)})
  assignment <- as.numeric(assignment)
  return(assignment)
}

wc_assign_chebyshev <- function(data, centroids)
{
  assignment <- apply(X = as.matrix(dist(rbind(data, centroids[, !grepl("WCCluster", colnames(centroids))]), method = 'maximum'))[1:nrow(data), (nrow(data) + 1):(nrow(data) + nrow(centroids))], MARGIN = 1, FUN = function(x) {which.min(x)})
  assignment <- as.numeric(assignment)
  return(assignment)
}

wc_assign_cosine <- function(data, centroids)
{
  distances <- (as.matrix(data) %*% as.matrix(t(centroids[, !grepl("WCCluster", colnames(centroids))])))/(as.matrix(apply(X = data, MARGIN = 1, FUN = function(x) {sqrt(sum(x^2))})) %*% as.matrix(t(apply(X = centroids, MARGIN = 1, FUN = function(x) {sqrt(sum(x^2))}))))
  assignment <- apply(X = distances, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_correlation <- function(data, centroids)
{
  distances <- 1 - (as.matrix(data - 1/ncol(data) * rowSums(data)) %*% as.matrix(t(centroids[, !grepl("WCCluster", colnames(centroids))] - 1/ncol(centroids[, !grepl("WCCluster", colnames(centroids[, !grepl("WCCluster", colnames(centroids))]))]) * rowSums(centroids[, !grepl("WCCluster", colnames(centroids))])))) / (sqrt(as.matrix(data - 1/ncol(data) * rowSums(data))^2) %*% sqrt(as.matrix(t(centroids[, !grepl("WCCluster", colnames(centroids))] - 1/ncol(centroids[, !grepl("WCCluster", colnames(centroids))]) * rowSums(centroids[, !grepl("WCCluster", colnames(centroids))]))^2)))
  assignment <- apply(X = distances, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_sorensen <- function(data, centroids)
{
  upper_matrix <- matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(abs(x))})})), nrow = nrow(data), byrow = TRUE)
  lower_matrix <- matrix(data = rep(x = rowSums(data), times = k), nrow = nrow(data), byrow = FALSE) + matrix(data = rep(x = t(rowSums(centroids[, !grepl("WCCluster", colnames(centroids))])), times = nrow(data)), nrow = nrow(data))

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_soergel <- function(data, centroids)
{
  upper_matrix <- matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(abs(x))})})), nrow = nrow(data), byrow = TRUE)
  lower_matrix <- t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(pmax(r, c))})}))

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_kulczynski <- function(data, centroids)
{
  upper_matrix <- matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(abs(x))})})), nrow = nrow(data), byrow = TRUE)
  lower_matrix <- t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(pmin(r, c))})}))

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_lorentzian <- function(data, centroids)
{
  dist_matrix <- matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(log(1 + abs(x)))})})), nrow = nrow(data), byrow = TRUE)

  assignment <- apply(X = dist_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_gower <- function(data, centroids)
{
  dist_matrix <- 1/ncol(data) * matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(abs(x))})})), nrow = nrow(data), byrow = TRUE)

  assignment <- apply(X = dist_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_intersection <- function(data, centroids)
{
  dist_matrix <- 0.5 * matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(abs(x))})})), nrow = nrow(data), byrow = TRUE)

  assignment <- apply(X = dist_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_czekanowski <- function(data, centroids)
{
  upper_matrix <- t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(pmin(r, c))})}))
  lower_matrix <- matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '+', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(x)})})), nrow = nrow(data), byrow = TRUE)

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_motika <- function(data, centroids)
{
  upper_matrix <- t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(pmax(r, c))})}))
  lower_matrix <- matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '+', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(x)})})), nrow = nrow(data), byrow = TRUE)

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_ruzicka <- function(data, centroids)
{
  upper_matrix <- t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(pmin(r, c))})}))
  lower_matrix <- t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(pmax(r, c))})}))

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_tanimoto <- function(data, centroids)
{
  upper_matrix <- t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(pmax(r, c) - pmin(r, c))})}))
  lower_matrix <- t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(pmax(r, c))})}))

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_inner_product <- function(data, centroids)
{
  assignment <- apply(X = as.matrix(data) %*% as.matrix(t(centroids[, !grepl("WCCluster", colnames(centroids))])), MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_jaccard_numerical <- function(data, centroids)
{
  upper_matrix <- matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(x^2)})})), nrow = nrow(data), byrow = TRUE)
  lower_matrix <- matrix(data = rep(x = rowSums(data^2), times = k), nrow = nrow(data), byrow = FALSE) + matrix(data = rep(x = t(rowSums(centroids[, !grepl("WCCluster", colnames(centroids))]^2)), times = nrow(data)), nrow = nrow(data)) - as.matrix(data) %*% as.matrix(t(centroids[, !grepl("WCCluster", colnames(centroids))]))

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_dice_numerical <- function(data, centroids)
{
  upper_matrix <- matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(x^2)})})), nrow = nrow(data), byrow = TRUE)
  lower_matrix <- matrix(data = rep(x = rowSums(data^2), times = k), nrow = nrow(data), byrow = FALSE) + matrix(data = rep(x = t(rowSums(centroids[, !grepl("WCCluster", colnames(centroids))]^2)), times = nrow(data)), nrow = nrow(data))

  assignment <- apply(X = upper_matrix / lower_matrix, MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_fidelity_numerical <- function(data, centroids)
{
  assignment <- apply(X = t(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(sqrt(r * c))})})), MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_bhattacharyya_numerical <- function(data, centroids)
{
  assignment <- apply(X = t(-1 * log(apply(X = data, MARGIN = 1, FUN = function(r) {apply(X = centroids[, !grepl("WCCluster", colnames(centroids))], MARGIN = 1, FUN = function(c) {sum(sqrt(r * c))})}))), MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_hellinger_numerical <- function(data, centroids)
{
  assignment <- apply(X = 2 * sqrt(matrix(unlist(lapply(X = apply(X = data, MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(x^2)})})), nrow = nrow(data), byrow = TRUE)), MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_whittaker <- function(data, centroids)
{
  assignment <- apply(X = matrix(unlist(lapply(X = apply(X = data/apply(X = data, MARGIN = 2, FUN = mean), MARGIN = 1, FUN = '-', centroids[, !grepl("WCCluster", colnames(centroids))]/apply(X = data, MARGIN = 2, FUN = mean)), FUN = function(x) {apply(X = x, MARGIN = 1, FUN = function(x) {sum(x^2)})})), nrow = nrow(data), byrow = TRUE), MARGIN = 1, FUN = function(x) {which.min(x)})
  return(assignment)
}

wc_assign_types <- data.frame()

wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Euclidean', 'Method' = 'wc_assign_euclidean'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Squared Euclidean', 'Method' = 'wc_assign_squared_euclidean'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Manhattan', 'Method' = 'wc_assign_manhattan'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Canberra', 'Method' = 'wc_assign_canberra'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Chebyshev', 'Method' = 'wc_assign_chebyshev'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Cosine', 'Method' = 'wc_assign_cosine'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Correlation', 'Method' = 'wc_assign_correlation'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Sorensen', 'Method' = 'wc_assign_sorensen'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Soergel', 'Method' = 'wc_assign_soergel'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Kulczynski', 'Method' = 'wc_assign_kulczynski'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Lorentzian', 'Method' = 'wc_assign_lorentzian'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Gower', 'Method' = 'wc_assign_gower'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Inersection', 'Method' = 'wc_assign_intersection'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Czekanowski', 'Method' = 'wc_assign_czekanowski'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Motika', 'Method' = 'wc_assign_motika'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Ruzicka', 'Method' = 'wc_assign_ruzicka'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Tanimoto', 'Method' = 'wc_assign_tanimoto'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Inner Product', 'Method' = 'wc_assign_inner_product'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Jaccard', 'Method' = 'wc_assign_jaccard_numerical'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Dice', 'Method' = 'wc_assign_dice_numerical'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Fidelity', 'Method' = 'wc_assign_fidelity_numerical'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Bhattacharyya', 'Method' = 'wc_assign_bhattacharyya_numerical'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Hellinger', 'Method' = 'wc_assign_hellinger_numerical'))
wc_assign_types <- rbind.data.frame(wc_assign_types, data.frame('Type' = 'Whittaker', 'Method' = 'wc_assign_whittaker'))
