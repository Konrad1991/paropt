AC_INIT(paropt, 0.2.1)

# Find the compiler and compiler flags used by R.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${PKG_LIBS} ${ac_cv_search_deflate}"])


openmp_flag=""


AC_MSG_CHECKING([for macOS])
RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
if test x"${RSysinfoName}" == x"Darwin"; then
   AC_MSG_RESULT([found])
   AC_MSG_WARN([OpenMP unavailable and turned off.])
   openmp_flag="-DARMA_DONT_USE_OPENMP"
else
   AC_MSG_RESULT([not found as on ${RSysinfoName}])
   AC_MSG_CHECKING([for OpenMP])
   allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
   hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)
   if test x"${hasOpenMP}" == x""; then
      AC_MSG_RESULT([missing])
      openmp_flag="-DARMA_DONT_USE_OPENMP"
   else
      AC_MSG_RESULT([found])
      openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
   fi
fi


AC_SUBST([openmp], ["${openmp_flag}"])


AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT




