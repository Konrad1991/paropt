optimize <- function(ode, lb, ub,
                     npop = 40, ngen = 10000,
                     reltol = 1e-06, abstol = 1e-08,
                     error = 0.0001,
                     states,
                     solvertype = "bdf",
                     own_error_fct,
                     own_spline_fct,
                     own_jac_fct,
                     number_threads,
                     verbose = FALSE) {

  stopifnot(!missing(ode))

  stopifnot(!missing(lb))
  stopifnot(!missing(ub))
  stopifnot(is.data.frame(lb))
  stopifnot(is.data.frame(ub))

  stopifnot(!missing(states))
  stopifnot(is.data.frame(states))

  assert_num_length_one(reltol, s(reltol, 1L))
  stopifnot("abstol has to be numeric" = is.numeric(abstol))
  assert_num_length_one(error, s(error, 1L))
  assert_num_length_one(npop, s(npop, 1L))
  assert_num_length_one(ngen, s(ngen, 1L))

  stopifnot("number of particles has to be > 0" = npop > 0)
  if (npop > 1000) warning("Unusually high number of particles.")
  stopifnot("number of generations has to be > 0" = ngen > 0)
  if (ngen > 10^5) warning("Unusually high number of generations")

  stopifnot("time has to be the first column in lb" = names(lb)[1] == "time")
  stopifnot("time has to be the first column in ub" = names(ub)[1] == "time")
  stopifnot("names have to be the same in ub and lb" = names(lb) == names(ub))
  stopifnot(
    "Found difference in dim() between lower and upper boundary" =
      identical(dim(lb), dim(ub))
  )

  stopifnot("time has to be the first column in states" = names(states)[1] == "time")
  integration_times <- states[[1]]

  assert_logical_length_one(verbose, s(verbose, 1L))

  # threads
  number_threads <- resolve_threads(number_threads)

  fct_ret <- resolve_ode_function(ode, verbose, optimizer = TRUE)
  ecf <- resolve_error_function(own_error_fct, verbose, optimizer = TRUE)
  sf <- resolve_spline_function(own_spline_fct, verbose, optimizer = TRUE)

  # own jac function
  stype <- NULL
  if (solvertype == "bdf") {
    stype <- 1
  } else if (solvertype == "adams") {
    stype <- 2
  }
  jf <- get_mock_jac_fct()
  if (!missing(own_jac_fct)) {
    if (stype == 2) {
      warning("own jacobian function cannot be used by solver adams. The function is ignored")
    } else if (is.function(own_jac_fct)) {
      stype <- 3
      jf <- resolve_jacobian(own_jac_fct, verbose, optimizer = TRUE)
    }
  }

  # boundaries
  par_time <- c()
  par_cut_idx <- c()
  lowb <- c()
  upb <- c()
  for (i in 2:ncol(lb)) {
    temp_lb <- lb[, i]
    temp_ub <- ub[, i]
    idx_lb <- !is.na(temp_lb)
    idx_ub <- !is.na(temp_ub)
    temp_lb <- temp_lb[idx_lb]
    temp_ub <- temp_ub[idx_ub]
    if (length(temp_ub) != length(temp_lb)) {
      message(paste("In column", i, "difference found between lb and ub"))
      stop("Error")
    }

    lowb <- c(lowb, temp_lb)
    upb <- c(upb, temp_ub)
    par_cut_idx <- c(par_cut_idx, length(temp_lb))

    lb_time <- lb[, 1]
    ub_time <- ub[, 1]
    lb_time <- lb_time[idx_lb]
    ub_time <- ub_time[idx_ub]
    stopifnot(
      "Found difference in time column between lower and upper boundary" =
        identical(lb_time, ub_time)
    )
    if (length(temp_ub) != length(lb_time)) {
      stopifnot("Found differences in number of entries between boundary and time. Maybe you have NA entries in your time column.")
    }
    par_time <- c(par_time, lb_time)
  }

  # states
  st <- c()
  for (i in 2:ncol(states)) {
    st <- c(st, states[, i])
  }
  state_idx_cuts <- rep(nrow(states), ncol(states) - 1L)

  # tolerances
  atol <- NULL
  if (missing(abstol)) {
    atol <- rep(1e-08, ncol(states) - 1L)
  } else {
    stopifnot(
      "Wrong number of absolute tolerances" =
        (ncol(states) - 1L) == length(abstol)
    )
    atol <- abstol
  }

  par_time <- as.vector(par_time)
  par_cut_idx <- as.integer(par_cut_idx)
  istate <- unlist(states[1, 2:ncol(states)])

  ret <- wrapper_optimizer(
    init_state = istate,
    par_times = par_time,
    param_idx_cuts = par_cut_idx,
    lb_ = lowb, ub_ = upb,
    state_measured = st, state_idx_cuts = state_idx_cuts,
    integration_times = integration_times,
    reltol, atol, fct_ret, npop, ngen,
    error, stype, ecf, sf, jf, number_threads
  )

  # states
  is_states <- data.frame(states$time, ret[[3]])
  names(is_states) <- names(states)

  # parameter
  indevidual_time_for_params <- list()

  params <- data.frame(matrix(NA, ncol = ncol(lb), nrow = length(lb$time)))
  params[, 1] <- lb$time
  counter <- 2
  increment <- 1
  curr_idx <- 1
  intermediate <- numeric(length(lb$time))
  p <- ret[[2]]
  for (i in seq_along(par_cut_idx)) {
    idx_time_used <- NULL
    counter_r <- 1
    increment <- par_cut_idx[i]
    r <- p[curr_idx:(curr_idx + increment - 1)]
    temp_lb <- lb[, i + 1] # because of time column
    idx_time_used <- !is.na(temp_lb)
    for (j in seq_along(idx_time_used)) {
      if (idx_time_used[j] == TRUE) {
        intermediate[j] <- r[counter_r]
        counter_r <- counter_r + 1
      } else if (idx_time_used[j] == FALSE) {
        intermediate[j] <- NA
      }
    }
    params[, counter] <- intermediate
    counter <- counter + 1
    curr_idx <- curr_idx + increment
  }
  names(params) <- names(lb)

  structure(
    list(
      global_best_error = ret[[1]],
      best_parameter_set = params,
      in_silico_states = is_states,
      original_states = states
    ), class = "OptimResPAROPT"
  )
}
