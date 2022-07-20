MA2 <- R6::R6Class("MA2",
    inherit = ast2ast:::MA,

    public = list(

      build2 = function(verbose = FALSE, reference = FALSE, name_f, args_f) {
        self$getast()
        self$ast2call()
        self$call2char()

        exporter_signature = "SEXP getXPtr();"
        signature = paste("void",
                          name_f,
                          "(double t, double* params, int size_params, double* y_, double* ydot_, int NEQ) {")
        pointer2sexp = paste("sexp",
                             args_f[2],
                             "(size_params, params, 2);",
                             "sexp",
                             args_f[3], "(NEQ, y_, 2);",
                             "sexp",
                             args_f[4], "(NEQ, ydot_, 2);")
        exporter = paste("SEXP getXPtr() {
                    typedef void (*funcPtr)(double t, double* params, int size_params, double* y_, double* ydot_, int NEQ);",
                    "return XPtr<funcPtr>(new funcPtr(&", name_f, ")); }")

        fct = c(
          exporter_signature,
          signature,
          pointer2sexp,
          self$vars_declaration(self$desired_type),
          self$char,
          '}',
          exporter
        )

        fct = paste( unlist(fct), collapse='')
      }


    )
)








#' Translation of R function into C++ code. 
#'
#' @export
#' @param f The function which should be translated from R to C++.
#' \strong{The function should accepts four arguments and does not return anything.}
#'  \enumerate{
#'    \item t = the independent variable most often the time
#'    \item parameter = the vector which contains the parameter at timepoint t (already splined)
#'    \item the state vector y with the current values
#'    \item the left hand side vector ydot which has to be filled. 
#'  }
#' @param verbose If set to TRUE the output of RcppXPtrUtils::cppXPtr or Rcpp::cppFunction is printed.
#' @return The external pointer to the generated C++ function.
#' @details \strong{The following types are supported: }
#'  \enumerate{
#'    \item numeric vectors
#'    \item numeric matrices
#'  }
#'  Variables can be either numeric vectors or matrices.
#'  Notably, it is possible that the variables change the type within the function.
#'  \strong{It is possible to declare a variable of a scalar numeric data type.
#'          This is done by adding '_db' to the end of the variable. Each time '_db' is found
#'          the variable is declared as a scalar numeric data type. In this case the
#'          object cannot change its type!}
#' @details \strong{The following functions are supported:}
#'  \enumerate{
#'    \item assignment: = and <-
#'    \item allocation: vector and matrix
#'    \item information about objects: length and dim
#'    \item Basic operations: +, -, *, /
#'    \item Indices: [] and at
#'    \item mathematical functions: sin, asin, sinh, cos, acos, cosh, tan, atan, tanh, log, ^ and exp
#'    \item concatenate objects: c
#'    \item control flow: for, if, else if, else
#'    \item comparison: ==, !=, >, <, >= and <=
#'    \item printing: print
#'    \item returning objects: return
#'    \item catmull-rome spline: cmr
#'    \item to get a range of numbers the ':' function can be used
#'  }
#' @details  \strong{Some details about the implemented functions}
#' @details  \itemize{
#'    \item allocation of memory works: Following forms are possible: vector(size_of_elements), vector(value, size_of_elements),
#'              matrix(nrows, ncols), matrix(value, nrows, ncols) and matrix(vector, nrows, ncols). The latter fills the matrix or the vector with the specified 'value'.
#'    \item For indices squared brackets '[]' can be used as common in R. Beyond that the function 'at' exists
#'              which accepts as first argument a variable and as the second argument you pass the desired index.
#'              The caveat of using 'at' is that only one entry can be accessed. The function '[]' can return more then one element.
#'              \strong{The 'at'function returns a reference to the vector entry.
#'                Therefore variable[index] can behave differently then at(variable, index).
#'                The function has to be used carefully when 'at' is used.
#'                Especially if '[]' and 'at' are mixed the function behaviour is difficult to predict.
#'                Please test it before using it in a serious project.}
#'    \item For-loops can be written as common in R
#'            \itemize{
#'                \item Nr.1 \cr
#'                      for(index in variable)\{ \cr
#'                        # do whatever \cr
#'                      \} \cr
#'                \item Nr.2 \cr
#'                      for(index in 1:length(variable)\{ \cr
#'                        # do whatever \cr
#'                      \} \cr
#'    }
#'    \item Be aware that it is not possible to assign the result of a comparison to a variable.
#'    \item The print function accepts either a scalar, vector, matrix, string, bool or nothing (empty line).
#'    \item In order to return an object use the 'return' function (The last object is not returned automatically as in R).
#'    \item In order to interpolate values the 'cmr' function can be used. The function needs three arguments.
#'          \enumerate{
#'            \item the first argument is the point of the independent variable (x) for which the dependent variable should be calculated (y). This has to be a vector of length one.
#'            \item the second argument is a vector defining the points of the independent variable (x). This has to be a vector of at least length four.
#'            \item the third argument is a vector defining the points of the dependent variable (y). This has to be a vector of at least length four.
#'        }
#'  }
#' \strong{Be aware that the R code is translated to ETR an expression template library which tries to mimic R.
#' However, it does not behave exactly like R! Please check your compiled function before using it in a serious project.}
#' @examples #Lotka-Volterra
#' \dontrun{
#' # Translating to R_fct
#'f <- function(time, par, y, ydot) {
#'  
#'  a_db = at(par, 1)
#'  b_db = at(par, 2)
#'  c_db = at(par, 3)
#'  d_db = at(par, 4)
#'  
#'  predator_db = at(y,1)
#'  prey_db = at(y, 2)
#'  
#'  ydot[1] = predator_db*prey_db*c_db - predator_db*d_db
#'  ydot[2] = prey_db*a_db - prey_db*predator_db*b_db
#'}
#' fptr <- paropt::convert(f, verbose = TRUE)
#' }
convert <- function(f, verbose = FALSE) {

    stopifnot(!missing(f))
    stopifnot(is.function(f))
    stopifnot(is.logical(verbose))
    
    # get name, args etc. check whether correct
    name_f <- as.character(substitute(f))
    args_f <- methods::formalArgs(f)
    stopifnot(length(args_f) == 4)
    stopifnot("params cannot be used as name as it is used internally in the C++ Code"=
                !("params") %in% args_f)
    stopifnot("size_params cannot be used as name as it is used internally in the C++ Code"=
                !("size_params") %in% args_f)
    stopifnot("y_ cannot be used as name as it is used internally in the C++ Code"=
                !("y_") %in% args_f)
    stopifnot("ydot_ cannot be used as name as it is used internally in the C++ Code"=
                !("ydot_") %in% args_f)
    stopifnot("NEQ cannot be used as name as it is used internally in the C++ Code"=
            !("NEQ") %in% args_f)

    desired_type = 'sexp'

    a = MA2$new(f, desired_type)
    fct <- a$build2(verbose, reference = TRUE, name_f, args_f)

    fct_ret = NULL

    tryCatch(
      expr = {
        fct_ret = Rcpp::cppFunction(code = fct, plugins = c("cpp17"),
                                         depends = c("ast2ast", "RcppArmadillo"),
                                         includes = "#include <etr.hpp>",
                                         verbose =  verbose)
      },
      error = function(e) {
        print(e)
        print("Sorry compilation failed!")
      }
    )

    return(fct_ret)
}
