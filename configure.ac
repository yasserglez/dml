AC_PREREQ([2.67])
LT_PREREQ([2.2.6b])

AC_INIT([dml], [0.1.0])
AM_PROG_AR
AM_INIT_AUTOMAKE([-Wall foreign subdir-objects])
LT_INIT

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99
AM_PROG_CC_C_O
PKG_PROG_PKG_CONFIG([0.25])

# Checks for libraries.
PKG_CHECK_MODULES([GLIB], [glib-2.0 >= 2.28.6])
PKG_CHECK_MODULES([GSL], [gsl >= 1.14])
PKG_CHECK_MODULES([IGRAPH], [igraph >= 0.6.5])

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_C_INLINE
if test "$ac_cv_c_inline" != no ; then
  AC_DEFINE([HAVE_INLINE], [1],
            [Define to 1 if the compiler supports the inline keyword.])
  AC_SUBST([HAVE_INLINE])
fi

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])
AC_FUNC_MALLOC

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile dml.pc])

AC_OUTPUT
