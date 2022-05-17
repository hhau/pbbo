design_from_crs2 <- function(discrep, param_set, n_design, n_crs2_iters) {
  init_val <- ParamHelpers::generateDesign(
    n = 1,
    par.set = param_set,
    fun = lhs::maximinLHS
  ) %>%
  unlist()

  par_names <- names(init_val)
  nlr_opt_wrap <- function(x) {
    l <- x
    names(l) <- par_names
    return(discrep(l))
  }

  futile.logger::flog.info("Starting stage one CRS2 optimiser")
  output_file <- tempfile()

  sink(file = output_file, append = TRUE)
  crs_res <- nloptr::nloptr(
    x0 = init_val,
    eval_f = nlr_opt_wrap,
    lb = ParamHelpers::getLower(param_set),
    ub = ParamHelpers::getUpper(param_set),
    opts = list(
      algorithm = "NLOPT_GN_CRS2_LM",
      print_level = 3,
      maxeval = n_crs2_iters,
      xtol_rel = -1
    )
  )
  sink(file = NULL)

  res <- process_output_file(output_file, param_set) %>%
    dplyr::filter(!is.na(y)) %>%
    resample_by_y(n_rows = n_design)

  file.remove(output_file)
  return(res)
}

begins_with_c <- function(x, c, return_index = FALSE) {
  chars <- strsplit(x, split = "")
  index <- lapply(chars, function(x_str) {
    any(x_str[1] == c)
  }) %>%
    unlist()

  if (return_index) {
    return(index)
  } else {
    x[index]
  }
}

process_output_file <- function(outfile, param_set) {
  output_contents <- readLines(con = outfile) %>%
    begins_with_c(c = c("i", "\t"))

  iteration_index <- output_contents %>%
    begins_with_c(c = "i", return_index = TRUE) %>%
    which()

  full_res <- lapply(iteration_index, function(x) {
    iteration_string <- output_contents[x]
    param_string <- output_contents[x + 1]
    loss_string <- output_contents[x + 2]

    iteration_num <- iteration_string %>%
      stringr::str_match("iteration: (\\d+)") %>%
      magrittr::extract(, 2) %>%
      as.numeric()

    param_vec <- param_string %>%
      stringr::str_match("\\tx = \\((.+)\\)") %>%
      magrittr::extract(, 2) %>%
      stringr::str_split(",") %>%
      unlist() %>%
      as.numeric()

    names(param_vec) <- ParamHelpers::generateDesign(n = 1, param_set) %>%
      unlist() %>%
      names()

    loss_num <- loss_string %>%
      stringr::str_match("\\tf\\(x\\) = ([+-]?([0-9]*[.])?[0-9]+)") %>%
      magrittr::extract(, 2) %>%
      as.numeric()

    res <- as.list(param_vec)
    res$y <- loss_num

    return(tibble::as_tibble(res))
  }) %>%
    dplyr::bind_rows()

  return(full_res)
}

resample_by_y <- function(full_design, n_rows) {
  log_raw_weights <- -exp(full_design$y)
  stopifnot(any(!is.na(log_raw_weights)))

  weights <- exp(log_raw_weights - matrixStats::logSumExp(log_raw_weights))
  indices <- sample(
    x = seq_len(nrow(full_design)),
    size = min(n_rows, nrow(full_design)),
    replace = FALSE,
    prob = weights
  )

  return(full_design[indices, ])
}
