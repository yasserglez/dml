# DML - Dependence Modeling Library
# Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

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
    src/vine/linkern.c \
    src/vine/linkern.h \
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
    tests/data/copula.R \
    tests/data/*.dat
