####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.0                                             #
#  June 30, 1999                                                  #
####################################################################
#
SHELL = /bin/sh
include ../make.inc
#
#  The machine (platform) identifier to append to the library names
#
#PLAT = _SUN4SOL2
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
#OPTS     = -O4 -u -f -mt
#OPTS     = -u -f -dalign -native -xO5 -xarch=v8plusa
#OPTS     = -u -f -dalign -native -xO5 -xarch=v8plusa
DRVOPTS  = $(OPTS)
NOOPT    = -u -f
#NOOPT    = -u -f -mt
LOADER   = $(CC)
#LOADOPTS = -mt
#LOADOPTS = -f -dalign -native -xO5 -xarch=v8plusa
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
#BLASLIB      = ../../blas$(PLAT).a
#BLASLIB     = -xlic_lib=sunperf_mt
#BLASLIB     = -xlic_lib=sunperf
LAPACKLIB    = liblapack.a
#TMGLIB       = tmglib$(PLAT).a
#EIGSRCLIB    = eigsrc$(PLAT).a
#LINSRCLIB    = linsrc$(PLAT).a



