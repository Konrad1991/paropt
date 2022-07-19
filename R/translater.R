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








#' Translation of R function into C++ code. See R package 'ast2ast' for details
#' @export
#' @param f The function which should be translated from R to C++.
#' @param verbose If set to TRUE the output of RcppXPtrUtils::cppXPtr is printed.
convert <- function(f, verbose = FALSE) {

    stopifnot(!missing(f))
    stopifnot(is.function(f))
    stopifnot(is.logical(verbose))
    
    # get name, args etc. check whether correct
    name_f <- as.character(substitute(f))
    args_f <- methods::formalArgs(f)
    stopifnot(length(args_f) == 4)

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
