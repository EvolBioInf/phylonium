AC_INIT([phylonium], [1.1], [kloetzl@volbio.mpg.de])
AC_CONFIG_SRCDIR([src/phylonium.cxx])
AM_INIT_AUTOMAKE([-Wall foreign ])
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Make sure, also the C++ programs are compiled with OpenMP
AC_LANG(C++)
AC_OPENMP

AX_CXX_COMPILE_STDCXX(14)

# Execute all tests using C
AC_LANG(C)
AC_OPENMP
AX_GCC_FUNC_ATTRIBUTE(ifunc)

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm], [], [have_gsl=no])
AC_CHECK_LIB([gsl],[gsl_ran_binomial], [], [have_gsl=no])

AS_IF([test "x$have_gsl" = "xno"],[
	AC_MSG_ERROR([Missing the Gnu Scientific Library.])
])

# The libdivsufsort header contains some Microsoft extension making
# compilation fail on certain systems (i.e. OS X). Add the following
# flag so the build runs smoothly.
CPPFLAGS="$CPPFLAGS -fms-extensions"
AC_CHECK_HEADERS([divsufsort.h],[have_libdivsufsort=yes],[have_libdivsufsort=no])
AC_CHECK_LIB(divsufsort, divsufsort, [], [have_libdivsufsort=no])

AS_IF([test "x$have_libdivsufsort" = "xno"],
	[AC_MSG_ERROR([Missing libdivsufsort.])
])

AC_ARG_ENABLE([x86simd],
	[AS_HELP_STRING([--disable-x86simd],[disable x86 SIMD support @<:@default: enabled@:>@])],
	[x86simd=${enableval}],[x86simd=enabled]
	)

# AC_MSG_WARN(["XXX${x86simd}"])
AS_IF([test "x${x86simd}" = xenabled],[AC_DEFINE([ENABLE_X86_SIMD],[1],[Enable x86 SIMD support])])
AM_CONDITIONAL([ENABLE_X86_SIMD],[test "x${x86simd}" = xenabled])

# older compilers don't support AVX512
AC_ARG_ENABLE([avx512],
	[AS_HELP_STRING([--disable-avx512],[disable AVX512 support @<:@default: enabled@:>@])],
	[avx512=${enableval}],[avx512=enabled]
	)

# AC_MSG_WARN(["XXX${avx512}"])
AS_IF([test "x${avx512}" = xenabled],[AC_DEFINE([ENABLE_AVX512],[1],[Enable AVX512 support])])
AM_CONDITIONAL([ENABLE_AVX512],[test "x${avx512}" = xenabled])

# Check for various headers including those used by libdivsufsort.
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h stdint.h inttypes.h err.h errno.h fcntl.h])

AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T
AC_HEADER_STDBOOL

AC_CHECK_FUNCS([floor pow sqrt strdup strerror])
AC_CHECK_FUNCS([strtoul strtod])
AC_CHECK_FUNCS([reallocarray])

AM_CONDITIONAL([HAVE_REALLOCARRAY], [test "x$ac_cv_func_reallocarray" = xyes])

AC_CONFIG_HEADERS([config.h:config.hin])

AC_CONFIG_FILES([
 libs/Makefile
 Makefile
 man/phylonium.1
 src/Makefile
 test/Makefile
])
AC_OUTPUT

