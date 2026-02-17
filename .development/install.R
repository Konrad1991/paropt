system('find -name "*.o" | xargs rm')
system('find -name "*.so" | xargs rm')
Rcpp::compileAttributes(".")
install.packages(".", repos = NULL, type = "source")
