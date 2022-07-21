AC_INIT([hexo], [2022.14.6], [merhardt@rhrk.uni-kl.de])
AM_INIT_AUTOMAKE([subdir-objects -Wpedantic -Wall -Wextra ])
AC_ARG_ENABLE([debug],
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],[debug=false])
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])
: ${CXXFLAGS="-fopenmp -O3 -march=native -fcx-fortran-rules -flto"}
AC_PROG_CXX
AC_PROG_FC
AX_CHECK_COMPILE_FLAG([-std=c++17], [
                        CXXFLAGS="$CXXFLAGS -std=c++17"])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
