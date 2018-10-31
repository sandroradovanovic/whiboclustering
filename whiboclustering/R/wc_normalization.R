wc_normalize <- function(data, normalization_type)
{
  # CHECKING FOR ERRORS
  if (!(tolower(normalization_type) %in% tolower(wc_norm_types$Type)))
  {
    stop('Please enter normalization function that is available in wc_norm_types data frame')
  }

  #PREPERANG DATA IF NECCESSARY
  if(sum(is.na(data)) > 0)
  {
    data <- as.data.frame(lapply(data, function(x) ifelse(is.na(x), mean(x, na.rm = T), x)))
    print("NA values present in data. These elements will be replaced with average of column")
  }

  if(sum(apply(data, 2, is.infinite)) > 0)
  {
    data <- as.data.frame(lapply(data, function(x) ifelse(is.infinite(x), NA, x)))
    data <- as.data.frame(lapply(data, function(x) ifelse(is.na(x), mean(x, na.rm = T), x)))

    print("Inf values present in data. These elements will be replaced with average of column")
  }

  if(sum(apply(data, 2, is.nan)) > 0)
  {
    data <- as.data.frame(lapply(data, function(x) ifelse(is.nan(x), mean(x, na.rm = T), x)))

    print("NaN values present in data. These elements will be replaced with average of column")
  }

  #NORMALIZING DATA
  norm_data <- eval(call(name = as.character(wc_norm_types$Method[tolower(wc_norm_types$Type) == tolower(normalization_type)]), data))

  return(norm_data)
}

wc_norm_no <- function(data, model = NULL)
{
  norm_data <- data

  output <- list('data' = norm_data, 'model' = NULL)

  return(output)
}

wc_norm_z <- function(data, model = NULL)
{
  if(is.null(model))
  {
    norm_data <- scale(x = data, center = TRUE, scale = TRUE)
    output <- list('data' = norm_data, 'model' = list('center' = attr(x = norm_data, which = 'scaled:center'), 'scale' = attr(x = norm_data, which = 'scaled:scale')))
  }
  else
  {
    norm_data <- scale(x = data, center = model$center, scale = model$scale)
    output <- list('data' = norm_data, 'model' = model)
  }


  return(output)
}

wc_norm_l2 <- function(data, model = NULL)
{
  if(is.null(model))
  {
    model <- c()
    norm_data <- apply(X = data, MARGIN = 2, FUN = function(x)
    {
      model <<- append(model, sqrt(sum(x^2, na.rm = TRUE)))
      x/sqrt(sum(x^2, na.rm = TRUE))
    })
  }
  else
  {
    norm_data <- sweep(x = data, MARGIN = 2, STATS = model, FUN = '/')
  }

  output <- list('data' = norm_data, 'model' = model)
  return(output)
}

wc_norm_l1 <- function(data, model = NULL)
{
  if(is.null(model))
  {
    model <- c()
    norm_data <- apply(X = data, MARGIN = 2, FUN = function(x)
    {
      model <<- append(model, sum(x, na.rm = TRUE))
      x/sum(x, na.rm = TRUE)
    })
  }
  else
  {
    norm_data <- sweep(x = data, MARGIN = 2, STATS = model, FUN = '/')
  }

  output <- list('data' = norm_data, 'model' = model)
  return(output)
}

wc_norm_linf <- function(data, model = NULL)
{
  if(is.null(model))
  {
    model <- c()
    norm_data <- apply(X = data, MARGIN = 2, FUN = function(x)
    {
      model <<- append(model, max(x, na.rm = TRUE))
      x/max(x, na.rm = TRUE)
    })
  }
  else
  {
    norm_data <- sweep(x = data, MARGIN = 2, STATS = model, FUN = '/')
  }

  output <- list('data' = norm_data, 'model' = model)
  return(output)
}

wc_norm_max_min <- function(data, model = NULL)
{
  if(is.null(model))
  {
    mins <- c()
    diffs <- c()
    norm_data <- apply(X = data, MARGIN = 2, FUN = function(x)
    {
      mins <<- append(mins, min(x, na.rm = TRUE))
      diffs <<- append(diffs, diff(range(x, na.rm = TRUE)))
      (x - min(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE))
    })
    output <- list('data' = norm_data, 'model' = list('mins' = mins, 'diffs' = diffs))
  }
  else
  {
    norm_data <- sweep(x = sweep(x = data, MARGIN = 2, STATS = model$mins, FUN = '-'), MARGIN = 2, STATS = model$diffs, FUN = '/')
    output <- list('data' = norm_data, 'model' = model)
  }
  return(output)
}

wc_norm_mean <- function(data, model = NULL)
{
  if(is.null(model))
  {
    means <- c()
    diffs <- c()
    norm_data <- apply(X = data, MARGIN = 2, FUN = function(x)
    {
      means <<- append(means, mean(x, na.rm = TRUE))
      diffs <<- append(diffs, diff(range(x, na.rm = TRUE)))
      (x - mean(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE))
    })
    output <- list('data' = norm_data, 'model' = list('mins' = mins, 'diffs' = diffs))
  }
  else
  {
    norm_data <- sweep(x = sweep(x = data, MARGIN = 2, STATS = model$means, FUN = '-'), MARGIN = 2, STATS = model$diffs, FUN = '/')
    output <- list('data' = norm_data, 'model' = model)
  }
  return(output)
}

wc_norm_log <- function(data, model = NULL)
{
  if(is.null(model))
  {
    model <- c()
    norm_data <- apply(X = data, MARGIN = 2, FUN = function(x)
    {
      model <<- append(model, log(x = prod(x, na.rm = TRUE)))
      log(x = x)/log(x = prod(x, na.rm = TRUE))
    })
  }
  else
  {
    norm_data <- sweep(x = log(data), MARGIN = 2, STATS = model, FUN = '/')
  }

  output <- list('data' = norm_data, 'model' = model)
  return(output)
}

wc_norm_non_monotonic <- function(data, model = NULL)
{
  if(is.null(model))
  {
    maxs <- c()
    vars <- c()
    norm_data <- apply(X = data, MARGIN = 2, FUN = function(x)
    {
      maxs <<- append(maxs, max(x, na.rm = TRUE))
      vars <<- append(vars, var(x, na.rm = TRUE))
      exp(-(x - max(x, na.rm = TRUE))^2/var(x, na.rm = TRUE))
    })
    output <- list('data' = norm_data, 'model' = list('maxs' = maxs, 'vars' = vars))
  }
  else
  {
    norm_data <- exp(-sweep(x = sweep(x = data, MARGIN = 2, STATS = model$maxs, FUN = '-')^2, MARGIN = 2, STATS = model$vars, FUN = '/'))
    output <- list('data' = norm_data, 'model' = model)
  }

  return(output)
}

wc_norm_comprehensive <- function(data, model = NULL)
{
  if(is.null(model))
  {
    mins <- c()
    maxs <- c()
    norm_data <- apply(X = data, MARGIN = 2, FUN = function(x)
    {
      mins <<- append(mins, min(x, na.rm = TRUE))
      maxs <<- append(maxs, max(x, na.rm = TRUE))
      1 - exp(abs(x - min(x, na.rm = TRUE))/(min(x, na.rm = TRUE) - max(x, na.rm = TRUE)))
    })
    output <- list('data' = norm_data, 'model' = list('mins' = mins, 'maxs' = maxs))
  }
  else
  {
    norm_data <- 1 - exp(sweep(x = abs(sweep(x = data, MARGIN = 2, STATS = model$mins, FUN = '-')), MARGIN = 2, STATS = model$mins - model$maxs, FUN = '/'))
    output <- list('data' = norm_data, 'model' = model)
  }

  return(output)
}

wc_norm_decimal_scaling <- function(data, model = NULL)
{
  if(is.null(model))
  {
    norm_data <- scale(x = data, center = FALSE, scale = 10^(ceiling(log10(apply(X = abs(data), MARGIN = 2, FUN = max)))))
    output <- list('data' = norm_data, 'model' = attr(x = norm_data, which = 'scaled:scale'))
  }
  else
  {
    norm_data <- scale(x = data, center = FALSE, scale = model)
    output <- list('data' = norm_data, 'model' = model)
  }

  return(output)
}

wc_norm_sigmoid <- function(data, model = NULL)
{
  if(is.null(model))
  {
    norm_data <- (1 - exp(-scale(x = data, center = TRUE, scale = TRUE)))/(1 + exp(-scale(x = data, center = TRUE, scale = TRUE)))
    output <- list('data' = norm_data, 'model' = list('center' = attr(x = norm_data, which = 'scaled:center'), 'scale' = attr(x = norm_data, which = 'scaled:scale')))
  }
  else
  {
    norm_data <- (1 - exp(-scale(x = data, center = model$center, scale = model$scale)))/(1 + exp(-scale(x = data, center = model$center, scale = model$scale)))
    output <- list('data' = norm_data, 'model' = model)
  }

  return(output)
}

wc_norm_softmax <- function(data, model = NULL)
{
  if(is.null(model))
  {
    norm_data <- 1/(1 + exp(-scale(x = data, center = TRUE, scale = TRUE)))
    output <- list('data' = norm_data, 'model' = list('center' = attr(x = norm_data, which = 'scaled:center'), 'scale' = attr(x = norm_data, which = 'scaled:scale')))
  }
  else
  {
    norm_data <- 1/(1 + exp(-scale(x = data, center = model$center, scale = model$scale)))
    output <- list('data' = norm_data, 'model' = model)
  }

  return(output)
}

wc_norm_types <- data.frame()

wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'No', 'Method' = 'wc_norm_no'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Z', 'Method' = 'wc_norm_z'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'L2', 'Method' = 'wc_norm_l2'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'L1', 'Method' = 'wc_norm_l1'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Linf', 'Method' = 'wc_norm_linf'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Max-Min', 'Method' = 'wc_norm_max_min'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Mean', 'Method' = 'wc_norm_mean'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Logarithmic', 'Method' = 'wc_norm_log'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Non-monotonic', 'Method' = 'wc_norm_non_monotonic'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Comprehensive', 'Method' = 'wc_norm_comprehensive'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Decimal Scaling', 'Method' = 'wc_norm_decimal_scaling'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Sigmoid', 'Method' = 'wc_norm_sigmoid'))
wc_norm_types <- rbind.data.frame(wc_norm_types, data.frame('Type' = 'Softmax', 'Method' = 'wc_norm_softmax'))
