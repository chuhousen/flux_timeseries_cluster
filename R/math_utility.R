na.median <-
  function(x) {
    ifelse(sum(!is.na(x)) >= 0.25 * min(length(x), 30), 
           median(x, na.rm = T), NA)
  }
upp.bd <-
  function(x) {
    ifelse(sum(!is.na(x)) >= 0.25 * min(length(x), 30),
           quantile(x, probs = 0.75, na.rm = T),
           NA)
  }
low.bd <-
  function(x) {
    ifelse(sum(!is.na(x)) >= 0.25 * min(length(x), 30),
           quantile(x, probs = 0.25, na.rm = T),
           NA)
  }
upp.bd2 <-
  function(x) {
    ifelse(sum(!is.na(x)) >= 0.25 * min(length(x), 30),
           quantile(x, probs = 0.975, na.rm = T),
           NA)
  }
low.bd2 <-
  function(x) {
    ifelse(sum(!is.na(x)) >= 0.25 * min(length(x), 30),
           quantile(x, probs = 0.025, na.rm = T),
           NA)
  }
na.mean <- function(x) {
  ifelse(sum(!is.na(x)) > 0, mean(x, na.rm = T), NA)
}
na.max <- function(x) {
  ifelse(sum(!is.na(x)) > 0, max(x, na.rm = T), NA)
}
na.min <- function(x) {
  ifelse(sum(!is.na(x)) > 0, min(x, na.rm = T), NA)
}
sum.na <- function(x) {
  sum(is.na(x), na.rm = T)
}
sum.notna <- function(x) {
  sum(!is.na(x), na.rm = T)
}
