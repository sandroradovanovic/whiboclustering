wc_recalculate <- function(data, assignment, recalculate_type)
{
  # CHECKING FOR ERRORS
  if (!(tolower(recalculate_type) %in% tolower(wc_recalculate_types$Type)))
  {
    stop('Please enter assignment function that is available in wc_assign_types data frame')
  }

  #RECALCULATE CLUSTER REPRESENTATIVE
  centroids <- eval(call(name = as.character(wc_recalculate_types$Method[tolower(wc_recalculate_types$Type) == tolower(recalculate_type)]), data, assignment))

  return(centroids)
}

wc_recalc_mean <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = mean)
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_median <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = median)
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_geomtric_mean <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = function(x) {prod(x)^(1/length(x))})
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_harmonic_mean <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = function(x) {length(x)/(sum(1/x))})
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_quadratic_mean <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = function(x) {sqrt(1/length(x) * sum(x^2))})
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_trimmed_mean <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = mean, trim = 0.05)
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_trimean <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = function(x) {0.5 * (quantile(x)[[3]] + (quantile(x)[[2]] + quantile(x)[[4]])/2)})
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_midhinge <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = function(x) {(quantile(x)[[2]] + quantile(x)[[4]])/2})
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_midrange <- function(data, assignment)
{
  new_centroids <- aggregate(x = data, by = list(assignment), FUN = function(x) {(quantile(x)[[1]] + quantile(x)[[5]])/2})
  colnames(new_centroids)[1] <- "WCCluster"
  return(new_centroids)
}

wc_recalc_online_mean <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}

wc_recalc_online_median <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}

wc_recalc_online_trimmed_mean <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}

wc_recalc_online_geometric_mean <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}


wc_recalc_online_harmonic_mean <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}

wc_recalc_online_quadratic_mean <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}

wc_recalc_online_trimean <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}

wc_recalc_online_midhinge <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}

wc_recalc_online_midrange <- function(data, assignment, assignment_type)
{
  assignments <- assignment
  for(i in 1:nrow(data))
  {
    if(i == 1)
    {
      new_centroids <- wc_recalc_mean(data = data, assignment = assignments)
    }
    new_assignment <- wc_assignment(data = data, centroids = new_centroids, assignment_type = assignment_type)
    new_centroids <- wc_recalc_mean(data = data, assignment = new_assignment)

    #IF ASSIGNMENT WAS EMPTY
    for(cent in 1:k)
    {
      if(nrow(new_centroids[new_centroids$WCCluster == cent, ]) == 0)
      {
        new_centroids <- rbind.data.frame(new_centroids, last_centroids[last_centroids$WCCluster == cent, ])
      }
      new_centroids <- new_centroids[order(new_centroids$WCCluster, decreasing = FALSE), ]
    }

    last_centroids <- new_centroids

    if(sum(new_assignment == assignments) == length(assignments))
    {
      break
    }
    else
    {
      assignments <- new_assignment
    }
  }
  return(new_centroids)
}

wc_recalculate_types <- data.frame()
assign(x = "wc_recalculate_types", value = wc_recalculate_types, envir = globalenv())

wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Mean', 'Method' = 'wc_recalc_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Median', 'Method' = 'wc_recalc_median'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Trimmed mean', 'Method' = 'wc_recalc_trimmed_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Geometric mean', 'Method' = 'wc_recalc_geomtric_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Harmonic mean', 'Method' = 'wc_recalc_harmonic_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Quadratic mean', 'Method' = 'wc_recalc_quadratic_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Trimean', 'Method' = 'wc_recalc_trimean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Midhinge', 'Method' = 'wc_recalc_midhinge'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Midrange', 'Method' = 'wc_recalc_midrange'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online mean', 'Method' = 'wc_recalc_online_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online median', 'Method' = 'wc_recalc_online_median'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online trimmed mean', 'Method' = 'wc_recalc_online_trimmed_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online geometric mean', 'Method' = 'wc_recalc_online_geometric_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online harmonic mean', 'Method' = 'wc_recalc_online_harmonic_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online Quadratic mean', 'Method' = 'wc_recalc_online_quadratic_mean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online trimean', 'Method' = 'wc_recalc_online_trimean'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online midhinge', 'Method' = 'wc_recalc_online_midhinge'))
wc_recalculate_types <- rbind.data.frame(wc_recalculate_types, data.frame('Type' = 'Online midrange', 'Method' = 'wc_recalc_online_midrange'))
