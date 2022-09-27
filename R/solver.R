solve <- function(ode, parameter,
                  reltol = 1e-06, abstol = 1e-08,
                  states,
                  solvertype = "bdf",
                  verbose = FALSE) {

  stopifnot(!missing(ode))
  stopifnot(is.function(ode))
  stopifnot(is.data.frame(parameter))
  stopifnot(is.data.frame(states))
  stopifnot(is.numeric(reltol))
  stopifnot(is.numeric(abstol))
  stopifnot("time has to be the first column in parameter" = names(parameter)[1] == "time")
  stopifnot("time has to be the first column in states" = names(states)[1] == "time")
  integration_times <- states[, 1]
  stopifnot(is.logical(verbose))
  name_f <- as.character(substitute(f))

  fct_ret <- ast2ast::translate(ode, verbose = verbose, output = "XPtr", reference = FALSE,
                                    types_of_args = c("double", rep("sexp", 3)),
                                    return_type = "void")

  # boundaries
  par_time <- c()
  par_cut_idx <- c()
  parb <- c()
  for(i in 2:dim(parameter)[2]) {
    temp_par <- parameter[, i]
    idx_par <- !is.na(temp_par)
    temp_par <- temp_par[idx_par]
    parb <- c(parb, temp_par)
    par_cut_idx <- c(par_cut_idx, length(temp_par))

    param_time <- parameter[, 1]
    param_time <- param_time[idx_par]
    par_time <- c(par_time, param_time)
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

  ret <- wrapper_solver(
    init_state = istate,
    par_times = par_time,
    param_idx_cuts = par_cut_idx,
    parameter_vec = parb,
    state_measured = st, state_idx_cuts = state_idx_cuts,
    integration_times = integration_times,
    reltol, atol, fct_ret, stype)

  # states
  is_states <- data.frame(states$time, ret[[2]])
  names(is_states) <- names(states)

  res <- list(error = ret[[1]], insilico_states = is_states)
  return(res)
}
