##  define global variables due to non-standard evaluations
utils::globalVariables(
         c("snp.id", "chr", "pos", "a1", "a2",
           "region.id", "region.start", "region.end",
           "p", "i.p", "chisq", "i.chisq",
           "a1_R", "a2_R", "a1_F", "a2_F", "a1_RF", "a2_RF",
           "patterns", "variable")
       )

stop2 <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

warning2 <- function(...) {
  warning(sprintf(...), call. = FALSE, immediate. = TRUE)
}

message2 <- function(...) {
  message(sprintf(...))
}

has_columns <- function(df, columns) {
  df_label <- deparse(substitute(df))
  df_columns <- names(df)
  diff_columns <- setdiff(columns, df_columns)
  if (length(diff_columns) > 0) {
    stop2("'%s' doesn't have column(s): %s.",
      df_label,
      paste(diff_columns, collapse = ", ")
    )
  }
}

check_class <- function(x, class) {
  x_label <- deparse(substitute(x))
  if (!inherits(x, class)) {
    stop2("%s is not an object of class: '%s'.", x_label, class)
  }
}

is_df <- function(df) {
  df_label <- deparse(substitute(df))
  if (!is.data.frame(df)) {
    stop2("'%s' is not a data frame.", df_label)
  }
}

is_tf <- function(x, arg) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    stop2("'%s' must be TRUE or FALSE.", arg)
  }
}

is_named_list <- function(x) {
  x_label <- deparse(substitute(x))
  x_names <- names(x)
  if (!is.list(x) || is.null(x_names) || anyDuplicated(x_names) > 0L) {
    stop2("'%s' must be a named list with an unique name for each set.",
          x_label)
  }
}

is_nonnegative_number <- function(x, arg) {
  if (!is.numeric(x) || is.na(x) || x < 0L || length(x) != 1L) {
    stop2("'%s' must be a non-negative number of length 1.", arg)
  }
}

is_number_between <- function(x, left, right, arg) {
  if (!is.numeric(x) || is.na(x) || x < left || x > right) {
    stop2("'%s' must be a number of length 1 between %s and %s.",
          arg, left, right)
  }
}

is_integer_vector <- function(x, tol = .Machine$double.eps) {
  ## this is an empirical solution
  is_wholenumber <- function(x, tol = tol) {
    ## function from is.integer help page
    abs(x - round(x)) < tol && !is.na(x)
  }
  tryCatch(all(sapply(x, is_wholenumber, tol = tol)),
           error = function(e) invisible(FALSE))
}

is_positive_definite <- function(ev, tol = 1e-7, symmetric = TRUE) {
  ## test the smallest of eigenvalue is negative or almost zero
  n <- length(ev)
  cutoff <- tol * abs(ev[1])
  invisible(ev[n] > cutoff)
}

get_duplicate_indice <- function(x) {
  which(duplicated(x) | duplicated(x, fromLast = TRUE))
}
