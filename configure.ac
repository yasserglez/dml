# DML - Dependence Modeling Library
# Copyright (C) 2011-2013 Yasser Gonzalez-Fernandez <ygonzalezfernandez@gmail.com>
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.67])
LT_PREREQ([2.2.6b])

AC_INIT([dml], [0.0.1])
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
PKG_CHECK_MODULES([IGRAPH], [igraph >= 0.5.3])

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
