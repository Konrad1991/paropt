resolve_threads <- function(number_threads) {
  if (missing(number_threads)) {
    number_threads <- RcppThread::detectCores()
  } else {
    stopifnot(is.numeric(number_threads))
    stopifnot(number_threads >= 1)
  }
  number_threads
}

s <- function(var, n = 1L) {
  env <- parent.frame(n)
  deparse1(evalq(substitute(var), env))
}

assert_num_length_one <- function(var, var_name) {
  if (!is.numeric(var)) stop(sprintf("%s is not of type numeric", var_name))
  if (length(var) != 1) stop(sprintf("%s, has not length 1", var_name))
}
assert_logical_length_one <- function(var, var_name) {
  if (!is.logical(var)) stop(sprintf("%s is not of type numeric", var_name))
  if (length(var) != 1) stop(sprintf("%s, has not length 1", var_name))
}

check_fct <- function(f, optimizer = TRUE) {
  e <- new.env()
  e$found <- FALSE
  walk_ast <- function(code) {
    if (!is.call(code)) {
      return(code)
    }
    code <- as.list(code)

    if (code[[1]] == as.name("return")) {
      e$found <- TRUE
      e$found_what <- c(e$found_what, code[[2]])
    } else if (deparse(code[[1]]) %in% e$not_thread_safe) {
      stop(paste("The function", deparse(code[[1]]), "is not thread safe."))
    }
    lapply(code, e$walk_ast)
  }
  e$walk_ast <- walk_ast
  if (optimizer == TRUE) {
    e$not_thread_safe <- c("print")
  } else {
    e$not_thread_safe <- NULL
  }

  e$found_what <- NULL
  trash <- e$walk_ast(body(f))
  stopifnot("Found a return statement in the ode-function" = e$found == FALSE)
}

check_fct_ode <- function(f, optimizer = TRUE) {
  stopifnot("ode system has to be a function" = is.function(f))
  check_fct(f, optimizer)
  args <- formalArgs(f)
  stopifnot("Four arguments have to be passed to ode-function!" = length(args) == 4)
}
resolve_ode_function <- function(ode, verbose, optimizer) {
  check_fct_ode(ode, optimizer = TRUE)
  fct_input <- ode
  var_names <- methods::formalArgs(ode)
  var_names[[1]] <- paste0(var_names[[1]], " |> type(double) |> ref()")
  var_names[[2]] <- paste0(var_names[[2]], " |> type(borrow_vec(double)) |> ref()")
  var_names[[3]] <- paste0(var_names[[3]], " |> type(borrow_vec(double)) |> ref()")
  var_names[[4]] <- paste0(var_names[[4]], " |> type(borrow_vec(double)) |> ref()")
  body <- paste0("{\n", paste(var_names, collapse = "\n"), "\n}")
  body(fct_input) <- str2lang(body)
  ast2ast::translate(
    ode, fct_input, verbose = verbose, output = "XPtr"
  )
}

check_fct_error_function <- function(f, optimizer) {
  stopifnot("error calculation function has to be a function" = is.function(f))
  check_fct(f, optimizer)
  args <- formalArgs(f)
  stopifnot("Three arguments have to be passed to the error calculation function!" = length(args) == 3)
}
resolve_error_function <- function(own_error_fct, verbose, optimizer) {
  if (missing(own_error_fct)) {
    get_default_error_fct()
  } else {
    check_fct_error_function(own_error_fct, optimizer)

    fct_input <- own_error_fct
    var_names <- methods::formalArgs(own_error_fct)
    var_names[[1]] <- paste0(var_names[[1]], " |> type(int)")
    var_names[[2]] <- paste0(var_names[[2]], " |> type(double)")
    var_names[[3]] <- paste0(var_names[[3]], " |> type(double)")
    body <- paste0("{\n", paste(var_names, collapse = "\n"), "\n}")
    body(fct_input) <- str2lang(body)

    ast2ast::translate(own_error_fct, fct_input, verbose = verbose, output = "XPtr")
  }
}

check_fct_spline_function <- function(f, optimizer) {
  stopifnot("spline function has to be a function" = is.function(f))
  check_fct(f, optimizer)
  args <- formalArgs(f)
  stopifnot("Three arguments have to be passed to the spline function!" = length(args) == 3)
}
resolve_spline_function <- function(own_spline_fct, verbose, optimizer) {
  if (missing(own_spline_fct)) {
    get_default_spline_fct()
  } else {
    check_fct_spline_function(own_spline_fct, optimizer)

    fct_input <- own_spline_fct
    var_names <- methods::formalArgs(own_spline_fct)
    var_names[[1]] <- paste0(var_names[[1]], " |> type(double) |> ref()")
    var_names[[2]] <- paste0(var_names[[2]], " |> type(vec(double)) |> ref()")
    var_names[[3]] <- paste0(var_names[[3]], " |> type(vec(double)) |> ref()")
    body <- paste0("{\n", paste(var_names, collapse = "\n"), "\n}")
    body(fct_input) <- str2lang(body)

    ast2ast::translate(own_spline_fct, fct_input, verbose)
  }
}

check_fct_jacobian <- function(f, optimizer) {
  stopifnot("jacobian function has to be a function" = is.function(f))
  check_fct(f, optimizer)
  args <- formalArgs(f)
  stopifnot("Five arguments have to be passed to the spline function!" = length(args) == 5)
}
resolve_jacobian <- function(own_jac_fct, verbose, optimizer) {
    check_fct_jacobian(own_jac_fct, optimizer)

    fct_input <- own_spline_fct
    var_names <- methods::formalArgs(own_spline_fct)
    var_names[[1]] <- paste0(var_names[[1]], " |> type(double) |> ref()")
    var_names[[2]] <- paste0(var_names[[2]], " |> type(borrow_vec(double)) |> ref()")
    var_names[[3]] <- paste0(var_names[[3]], " |> type(borrow_vec(double)) |> ref()")
    var_names[[4]] <- paste0(var_names[[3]], " |> type(borrow_mat(double)) |> ref()")
    var_names[[5]] <- paste0(var_names[[3]], " |> type(borrow_vec(double)) |> ref()")
    body <- paste0("{\n", paste(var_names, collapse = "\n"), "\n}")
    body(fct_input) <- str2lang(body)

    ast2ast::translate(own_jac_fct, fct_input, verbose)
}
