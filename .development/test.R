ode <- function(t, y, ydot, parameter) {
  a_db <- at(parameter, 1)
  b_db <- at(parameter, 2)
  c_db <- at(parameter, 3)
  d_db <- at(parameter, 4)
  predator_db <- at(y, 1)
  prey_db <- at(y, 2)
  ydot[1] <- predator_db * prey_db * c_db - predator_db * d_db
  ydot[2] <- prey_db * a_db - prey_db * predator_db * b_db
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
  states = states
)

res

insilico <- res[[3]]
par(mfrow = c(2, 1))
plot(states$time, states$n1, type = "l")
points(insilico$time, insilico$n1, type = "l", col = "darkred")

plot(states$time, states$n2, type = "l")
points(insilico$time, insilico$n2, type = "l", col = "darkred")
