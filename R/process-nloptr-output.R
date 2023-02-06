design_from_crs2 <- function(
  discrep,
  param_set,
  n_design,
  n_design_pad,
  n_crs2_iters,
  crs2_ranseed
) {
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

  crs_res <- run_nlopt_crs2(
    output_file = output_file,
    init_val = init_val,
    nlr_opt_wrap = nlr_opt_wrap,
    param_set = param_set,
    n_crs2_iters = n_crs2_iters,
    crs2_ranseed = crs2_ranseed
  )

  res <- process_output_file(output_file, param_set) %>%
    dplyr::filter(!is.na(y)) %>%
    resample_by_y(n_rows = n_design)

  # sometimes output from crs2 is insufficiently diverse, so add extra points
  # because otherwise everything downstream breaks. 10 is just a magic number
  base_desgin <- ParamHelpers::generateDesign(
    n = n_design_pad,
    par.set = param_set,
    fun = lhs::maximinLHS
  )

  base_desgin$y <- sapply(1 : nrow(base_desgin), function(row_id) {
    base_desgin[row_id,] %>%
      unlist() %>%
      discrep()
  })

  final_design <- dplyr::bind_rows(res, base_desgin)

  file.remove(output_file)
  return(final_design)
}

run_nlopt_crs2 <- function(
  output_file,
  init_val,
  nlr_opt_wrap,
  param_set,
  n_crs2_iters,
  crs2_ranseed
) {
  sink(file = output_file, append = TRUE)
  on.exit(expr = sink(file = NULL), add = TRUE, after = FALSE)

  if (is.null(crs2_ranseed)) {
    crs2_ranseed <- 0
  }

  crs_res <- nloptr::nloptr(
    x0 = init_val,
    eval_f = nlr_opt_wrap,
    lb = ParamHelpers::getLower(param_set),
    ub = ParamHelpers::getUpper(param_set),
    opts = list(
      algorithm = "NLOPT_GN_CRS2_LM",
      print_level = 3,
      maxeval = n_crs2_iters,
      xtol_rel = -1,
      ranseed = crs2_ranseed
    )
  )

  return(crs_res)
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
  # need to be careful of overflow in the exp() step -- filter anything that
  # will overflow. 700 ~= log(.Machine$double.xmax)
  feasible_design <- full_design %>%
    dplyr::filter(y < 700)

  log_raw_weights <- -exp(feasible_design$y)
  stopifnot(any(!is.na(log_raw_weights)))
  weights <- exp(log_raw_weights - matrixStats::logSumExp(log_raw_weights))
  n_feasible <- sum(weights > 0)

  if (n_feasible < n_rows) {
    futile.logger::flog.info(
      paste(
        "Not enough feasible points from CRS2 optimisation path",
        "to satisfy the requested number of design points.",
        "Consider using more CRS2 iterations.",
        "Optimisation will continue using only the feasible points."
      )
    )
  }

  indices <- sample(
    x = seq_len(nrow(feasible_design)),
    size = min(n_rows, n_feasible),
    replace = FALSE,
    prob = weights
  )

  return(feasible_design[indices, ])
}
