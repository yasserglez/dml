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

AM_CPPFLAGS = @GLIB_CFLAGS@ @GSL_CFLAGS@ @IGRAPH_CFLAGS@
AM_CFLAGS = -Wall -Werror -pedantic
AM_LDFLAGS = @GLIB_LIBS@ @GSL_LIBS@ @IGRAPH_LIBS@

ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dml.pc

lib_LTLIBRARIES = libdml.la

check_PROGRAMS = test
TESTS = test

libdml_la_SOURCES = \
    src/measure/measure.c \
    src/measure/rank.c \
    src/measure/tau.c \
    src/measure/cvm.c \
    src/copula/copula.c \
    src/copula/indep.c \
    src/copula/normal.c \
    src/copula/clayton.c \
    src/copula/select.c \
    src/vine/vine.c \
    src/vine/cvine.c \
    src/vine/dvine.c \
    src/vine/rvine.c    

include_HEADERS = src/dml.h

test_SOURCES = \
    tests/common.h \
    tests/common.c \
    tests/main.c \
    tests/measure.c \
    tests/copula.c \
    tests/vine.c    

test_LDADD = libdml.la

EXTRA_DIST = \
    README.txt \
    LICENSE.txt \
    src/measure/cvm.R \
    tests/data/copula.R \
    tests/data/*.dat
