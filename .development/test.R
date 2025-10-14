ode <- function(t, y, ydot, parameter) {
  a <- parameter[[1]]
  b <- parameter[[2]]
  c <- parameter[[3]]
  d <- parameter[[4]]
  predator <- y[[1]]
  prey <- y[[2]]
  ydot[1] <- predator * prey * c - predator * d
  ydot[2] <- prey * a - prey * predator * b
}

path <- system.file("examples", package = "paropt")
states <- read.table(paste(path, "/states_LV.txt", sep = ""), header = TRUE)
lb <- data.frame(time = 0, a = 0.8, b = 0.3, c = 0.09, d = 0.09)
ub <- data.frame(time = 0, a = 1.3, b = 0.7, c = 0.4, d = 0.7)
set.seed(1)

res <- paropt::optimize(ode,
  lb = lb, ub = ub,
  reltol = 1e-06, abstol = c(1e-08, 1e-08),
  error = 0.0001,
  npop = 40, ngen = 1000,
  states = states,
  verbose = TRUE
)
summary <- function(object, ...) {
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
summary(res)
