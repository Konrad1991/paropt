optimize <- function(ode, lb, ub,
                     npop = 40, ngen = 10000,
                     reltol = 1e-06, abstol = 1e-08,
                     error = 0.0001,
                     states,
                     solvertype = "bdf",
                     verbose = FALSE) {

    stopifnot(!missing(ode))
    stopifnot(is.function(ode))
    stopifnot(is.data.frame(lb))
    stopifnot(is.data.frame(ub))
    stopifnot(is.data.frame(states))
    stopifnot(is.numeric(reltol))
    stopifnot(is.numeric(abstol))
    stopifnot(is.numeric(error))
    stopifnot("Error has to be >= 0" = error >= 0)
    stopifnot(is.atomic(npop))
    stopifnot(is.atomic(ngen))
    stopifnot("number of particles has to be > 0" = npop > 0)
    stopifnot("number of generations has to be > 0" = ngen > 0)
    stopifnot("time has to be the first column in lb" = names(lb)[1] == "time")
    stopifnot("time has to be the first column in ub" = names(ub)[1] == "time")
    stopifnot("time has to be the first column in states" = names(states)[1] == "time")
    integration_times <- states[, 1]
    stopifnot("names have to be the same in ub and lb" = names(lb) == names(ub) )
    stopifnot(is.logical(verbose))
    name_f <- as.character(substitute(f))
    fct_ret <- ast2ast::translate(ode, verbose = verbose, output = "XPtr", reference = FALSE,
                                  types_of_args = c("double", rep("sexp", 3)),
                                  return_type = "void")
    stopifnot("Found differnce in dim() between lower and upper boundary" =
                identical(dim(lb), dim(ub)))

    # boundaries
    par_time <- c()
    par_cut_idx <- c()
    lowb <- c()
    upb <- c()
    for(i in 2:dim(lb)[2]) {
      temp_lb <- lb[, i]
      temp_ub <- ub[, i]
      idx_lb <- !is.na(temp_lb)
      idx_ub <- !is.na(temp_ub)
      temp_lb <- temp_lb[idx_lb]
      temp_ub <- temp_ub[idx_ub]
      if(length(temp_ub) != length(temp_lb)) {
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
      stopifnot("Found differnce in time column between lower and upper boundary" =
                  identical(lb_time, ub_time))
      if(length(temp_ub) != length(lb_time)) {
        stopifnot("Found differences in number of entries between boundary and time. Maybe you have NA entries in your time column.")
      }
      par_time <- c(par_time, lb_time)
    }

    # states
    st <- c()
    for(i in 2:dim(states)[2]) {
      st <- c(st, states[, i])
    }
    state_idx_cuts <- rep(dim(states)[1], dim(states)[2] -1)

    # tolerances
    atol <- NULL
    if(missing(abstol)) {
      atol <- rep(1e-08, dim(states)[2] - 1)
    } else {
      stopifnot("Wrong number of absolute tolerances" =
                  (dim(states)[2] - 1) == length(abstol) )
      atol = abstol
    }

    stype <- NULL
    if(solvertype == "bdf") {
      stype <- 1
    } else if(solvertype == "adams") {
      stype <- 2
    }

    par_time <- as.vector(par_time)
    par_cut_idx <- as.integer(par_cut_idx)
    istate <- unlist(states[1, 2:dim(states)[2]])

    ret <- wrapper_optimizer(
                      init_state = istate,
                      par_times = par_time,
                      param_idx_cuts = par_cut_idx,
                      lb_ = lowb, ub_ = upb,
                      state_measured = st, state_idx_cuts = state_idx_cuts,
                      integration_times = integration_times,
                      reltol, atol, fct_ret, npop, ngen,
                      error, stype)

    # states
    is_states <- data.frame(states$time, ret[[3]])
    names(is_states) <- names(states)

    # parameter
    indevidual_time_for_params <- list()

    params <- data.frame(matrix(NA, ncol = dim(lb)[2], nrow = length(lb$time)))
    params[,1] <- lb$time
    counter <- 2
    increment <- 1
    curr_idx <- 1
    intermediate <- numeric(length(lb$time))
    p <- ret[[2]]
    for(i in seq_along(par_cut_idx)) {
      idx_time_used <- NULL
      counter_r <- 1
      increment <- par_cut_idx[i]
      r <- p[curr_idx:(curr_idx + increment - 1)]
      temp_lb <- lb[, i + 1] # because of time column
      idx_time_used <- !is.na(temp_lb)
      for(j in seq_along(idx_time_used)) {
        if(idx_time_used[j] == TRUE) {
          intermediate[j] = r[counter_r]
          counter_r <- counter_r + 1
        } else if(idx_time_used[j] == FALSE) {
          intermediate[j] = NA
        }
      }
      params[, counter] <- intermediate
      counter <- counter + 1
      curr_idx <- curr_idx + increment
    }
    names(params) <- names(lb)

    # return
    res <- list(
      global_best_error = ret[[1]],
      best_parameter_set = params,
      in_silico_states = is_states
    )

    return(res)
}