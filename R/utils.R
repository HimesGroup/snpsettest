##  define global variables due to non-standard evaluations
utils::globalVariables(
         c(".", "id", "chr", "pos", "A1", "A2",
           "gene.id", "gene.start", "gene.end",
           "pvalue", "i.pvalue", "chisq", "i.chisq",
           "A1_S", "A2_S", "A1_F", "A2_F", "A1_SF", "A2_SF",
           "patterns", "variable", "maf", "callrate",
           "key_", "swapped_", "flipped_", "tag")
       )

has_columns <- function(df, columns) {
  df_name <- deparse(substitute(df))
  df_columns <- names(df)
  diff_columns <- setdiff(columns, df_columns)
  if (length(diff_columns) > 0) {
    stop("'", df_name, "' doesn't have column(s): ",
         paste(diff_columns, collapse = ", "), ".",
         call. = FALSE
    )
  }
}

assert_class <- function(obj, class) {
  obj_name <- deparse(substitute(obj))
  if (!inherits(obj, class)) {
    stop(obj_name, " is not an object of class: '", class, "'",
         call. = FALSE)
  }
}

is_bed_matrix <- function(obj) {
  obj_name <- deparse(substitute(obj))
  ## Possibly `methods::is()` function would be more appropriate since
  ## bed.matrix is S4 class
  if (!inherits(obj, "bed.matrix")) {
    stop(obj_name, " is not a bed.matrix", call. = FALSE)
  }
}

is_df <- function(df) {
  df_name <- deparse(substitute(df))
  if (!is.data.frame(df)) {
    stop("'", df_name, "' is not a data frame.",
         call. = FALSE)
  }
}

is_tf <- function(x, arg_name) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    stop("'", arg_name, "' must be TRUE or FALSE.", call. = FALSE)
  }
}

is_named_list <- function(x) {
  x_name <- deparse(substitute(x))
  x_colnames <- names(x)
  if (!is.list(x) || is.null(x_colnames) || anyDuplicated(x_colnames) > 0L) {
    stop("'", x_name,
         "' must be a named list with an unique name for each set.",
         call. = FALSE)
  }
}

is_nonnegative_number <- function(x, arg) {
  if (!is.numeric(x) || is.na(x) || x < 0L || length(x) != 1L) {
    stop("'", arg, "' must be a non-negative number of length 1.",
         call. = FALSE)
  }
}

is_number_between <- function(x, left, right, arg) {
  if (!is.numeric(x) || is.na(x) || x < left || x > right) {
    stop("'", arg, "' must be a number of length 1 between ",
         left, " and ", right, ".", call. = FALSE)
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

pretty_num <- function(x, ...) {
  prettyNum(x, big.mark = ",", scientific = FALSE, ...)
}

is_positive_definite <- function(ev, tol = 1e-7) {
  ## test the smallest of eigenvalue is negative or almost zero
  n <- length(ev)
  cutoff <- tol * abs(ev[1])
  invisible(ev[n] > cutoff)
}

get_duplicate_indice <- function(x) {
  which(duplicated(x) | duplicated(x, fromLast = TRUE))
}
