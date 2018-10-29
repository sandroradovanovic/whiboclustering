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

wc_norm_no <- function(data)
{
  norm_data <- data
  return(norm_data)
}

wc_norm_z <- function(data)
{
  norm_data <- scale(x = data, center = TRUE, scale = TRUE)
  return(norm_data)
}

wc_norm_l2 <- function(data)
{
  norm_data <- apply(X = data, MARGIN = 2, FUN = function(x) {if(min(x, na.rm = TRUE) > 0) {x/sqrt(sum(x^2, na.rm = TRUE))} else {(x - min(x, na.rm = TRUE))/sqrt(sum((x - min(x, na.rm = TRUE) + 0.01)^2, na.rm = TRUE))} })
  return(norm_data)
}

wc_norm_l1 <- function(data)
{
  norm_data <- apply(X = data, MARGIN = 2, FUN = function(x) {if(min(x, na.rm = TRUE) > 0) {x/sum(x, na.rm = TRUE)} else {(x - min(x, na.rm = TRUE))/(sum(x - min(x, na.rm = TRUE) + 0.01, na.rm = TRUE))}})
  return(norm_data)
}

wc_norm_linf <- function(data)
{
  norm_data <- apply(X = data, MARGIN = 2, FUN = function(x) {if(min(x, na.rm = TRUE) > 0) {x/max(x, na.rm = TRUE)} else {(x - min(x, na.rm = TRUE))/(max(x - min(x, na.rm = TRUE) + 0.01, na.rm = TRUE))}})
  return(norm_data)
}

wc_norm_max_min <- function(data)
{
  norm_data <- apply(X = data, MARGIN = 2, FUN = function(x) {(x - min(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE)) + 0.01})
  return(norm_data)
}

wc_norm_mean <- function(data)
{
  norm_data <- apply(X = data, MARGIN = 2, FUN = function(x) {(x - mean(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE)) + 0.01})
  return(norm_data)
}

wc_norm_log <- function(data)
{
  norm_data <- apply(X = data, MARGIN = 2, FUN = function(x) {if(min(x, na.rm = TRUE) > exp(1)) {log(x = x)/log(x = prod(x, na.rm = TRUE))} else {log(x = (x - min(x, na.rm = TRUE) + 0.01 + exp(1)))/log(x = prod(x - min(x, na.rm = TRUE) + 0.01 + exp(1)))}})
  return(norm_data)
}

wc_norm_non_monotonic <- function(data)
{
  norm_data <- apply(X = data, MARGIN = 2, FUN = function(x) {exp(-(x - max(x, na.rm = TRUE) + 0.01)^2/var(x, na.rm = TRUE))})
  return(norm_data)
}

wc_norm_comprehensive <- function(data)
{
  norm_data <- apply(X = data, MARGIN = 2, FUN = function(x) {1 - exp(abs(x - min(x, na.rm = TRUE))/(min(x, na.rm = TRUE) - max(x, na.rm = TRUE)))})
  return(norm_data)
}

wc_norm_decimal_scaling <- function(data)
{
  norm_data <- scale(x = data, center = FALSE, scale = 10^(ceiling(log10(apply(X = abs(data), MARGIN = 2, FUN = max)))))
  return(norm_data)
}

wc_norm_sigmoid <- function(data)
{
  norm_data <- (1 - exp(-scale(x = data, center = TRUE, scale = TRUE)))/(1 + exp(-scale(x = data, center = TRUE, scale = TRUE)))
  return(norm_data)
}

wc_norm_softmax <- function(data)
{
  norm_data <- 1/(1 + exp(-scale(x = data, center = TRUE, scale = TRUE)))
  return(norm_data)
}

wc_norm_types <<- data.frame()

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
