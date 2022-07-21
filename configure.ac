AC_INIT([hexo], [2022.14.6], [merhardt@rhrk.uni-kl.de])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall])
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
#AC_CHECK_LIB(blas, )
#AC_CHECK_LIB(curl, )
#AC_CHECK_LIB(eigen3, )
#AC_CHECK_LIB(f2c, )
# Check for curl
AC_CHECK_LIB(curl, curl_easy_init, curl_preinstalled=true, curl_preinstalled=false)
AC_CHECK_LIB(curl, curl_slist_append,                    , curl_preinstalled=false)
AC_CHECK_LIB(curl, curl_easy_setopt,                     , curl_preinstalled=false)
AC_CHECK_LIB(curl, curl_easy_perform,                    , curl_preinstalled=false)
AC_CHECK_LIB(curl, curl_easy_cleanup,                    , curl_preinstalled=false)
AC_CHECK_LIB(curl, curl_slist_free_all,                  , curl_preinstalled=false)
AM_CONDITIONAL([NO_CURL], [test x$curl_preinstalled = xfalse])

# Check for FFTW3
AC_CHECK_LIB(fftw3, fftw_plan_dft_r2c_1d, fftw_preinstalled=true, fftw_preinstalled=false)
AC_CHECK_LIB(fftw3, fftw_plan_dft_1d,                           , fftw_preinstalled=false)
AC_CHECK_LIB(fftw3, fftw_malloc ,                               , fftw_preinstalled=false)
AC_CHECK_LIB(fftw3, fftw_free ,                                 , fftw_preinstalled=false)
AC_CHECK_LIB(fftw3, fftw_execute ,                              , fftw_preinstalled=false)
AC_CHECK_LIB(fftw3, fftw_destroy_plan,                          , fftw_preinstalled=false)
AM_CONDITIONAL([NO_FFTW], [test x$fftw_preinstalled = xfalse])

# Check for sqlite3
AC_CHECK_LIB(sqlite3, sqlite3_column_double, sqlite_preinstalled=true, sqlite_preinstalled=false)
AC_CHECK_LIB(sqlite3, sqlite3_prepare_v2,                            , sqlite_preinstalled=false)
AC_CHECK_LIB(sqlite3, sqlite3_bind_text,                             , sqlite_preinstalled=false)
AC_CHECK_LIB(sqlite3, sqlite3_bind_double,                           , sqlite_preinstalled=false)
AC_CHECK_LIB(sqlite3, sqlite3_bind_parameter_index,                  , sqlite_preinstalled=false)
AC_CHECK_LIB(sqlite3, sqlite3_step,                                  , sqlite_preinstalled=false)
AC_CHECK_LIB(sqlite3, sqlite3_errmsg,                                , sqlite_preinstalled=false)
AC_CHECK_LIB(sqlite3, sqlite3_finalize,                              , sqlite_preinstalled=false)
AM_CONDITIONAL([NO_SQLITE], [test x$sqlite_preinstalled = xfalse])

#AC_CHECK_HEADERS([eigen3/Eigen/Dense])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
