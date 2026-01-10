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

summary.OptimResPAROPT <- function(object, ...) {
  stopifnot(inherits(object, "OptimResPAROPT"))
  cat(sprintf("Best error: %s", object$global_best_error), "\n")
  cat("Parameters:\n")
  print(object$best_parameter_set)

  cat("In silico states:\n")
  print(object$in_silico_states[1:6, ])
  if (nrow(object$in_silico_states) > 6) {
    cat("\t.\n\t.\n\t.\n")
  }
  cat("True states:\n")
  print(object$original_states[1:6, ])
  if (nrow(object$original_states) > 6) {
    cat("\t.\n\t.\n\t.\n")
  }
  invisible(NULL)
}
summary.SolverResPAROPT <- function(object, ...) {
  stopifnot(inherits(object, "SolverResPAROPT"))
  cat(sprintf("Error: %s", object$error), "\n")
  cat("States:\n")
  print(object$in_silico_states[1:6, ])
  if (nrow(object$in_silico_states) > 6) {
    cat("\t.\n\t.\n\t.\n")
  }
  cat("True states:\n")
  print(object$original_states[1:6, ])
  if (nrow(object$original_states) > 6) {
    cat("\t.\n\t.\n\t.\n")
  }
  invisible(NULL)
}
