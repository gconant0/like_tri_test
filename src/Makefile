# Master Makefile for Sequence Utilities for Genome analysis
# Type make to build
include ../make.inc

O = o
SRC_DIR = ./

LIBRARY_DIR = -L./ -L../lapack -L../libf2c
INCLUDE = -I. -I../lapack
CFLAGS = -g                 		                     # compiler switches to be applied to every module
OPTIM_SPEED = -O3             	  	                     # switches that give speed priority over size
OPTIM_SIZE = -O1              	  	                     # switches that give size priority over speed
MATH_LIB = -lm                  	                     # math library 
F2C_LIB =  -lf2c                                             # F2C library need for legacy ranlib and lapack
LAPACK_LIB = -llapack                                        # Linear Algebra library
DNAFUNCS_LIB= -ldnafuncs
LIKELIHOOD_LIB = -llikelihood

OPTIONS = $(CFLAGS) $(INCLUDE)


all:  libdnafuncs liblikelihood like_tri_test 

DNALIBRARY = libdnafuncs.a

DNALIBRARY_OBJS =		\
	read_seq.$(O)	\
	gen_dna_funcs.$(O) \
	score_matrix.$(O) \
	exchange.$(O)   \
	gen_code.$(O)   \
	write_seq.$(O) 

libdnafuncs : $(DNALIBRARY) 

$(DNALIBRARY): $(DNALIBRARY_OBJS) 
	ar  cr $(DNALIBRARY) $(DNALIBRARY_OBJS)
	ranlib $(DNALIBRARY)


LIKELIHOODLIBRARY = liblikelihood.a

LIKELIHOODLIBRARY_OBJS =		\
	maxlike.$(O)	\
	codon_like.$(O) \
	nucleotide_like.$(O)   \
	tree.$(O)          \
	powell.$(O)       \
	lin_alg.$(O)  \
	other_like.$(O) \
	genome_list.$(O)


liblikelihood : $(LIKELIHOODLIBRARY) 

$(LIKELIHOODLIBRARY): $(LIKELIHOODLIBRARY_OBJS)
	ar cr $(LIKELIHOODLIBRARY) $(LIKELIHOODLIBRARY_OBJS) 
	ranlib  $(LIKELIHOODLIBRARY)


LIKE_TRI_TEST_OBJS = nrutil.$(O)	\
		      read_tree.$(O)    \
	              pair_dist.$(O)    \
	              like_tri_test.$(O)


like_tri_test: $(LIKE_TRI_TEST_OBJS) $(LIKELIHOODLIBRARY) $(DNALIBRARY)
	$(CC) $(LINUX_BUILD)  $(LIBRARY_DIR) \
	 -o ../like_tri_test $(OPTIONS) $(LIKE_TRI_TEST_OBJS) $(DNAFUNCS_LIB) $(LIKELIHOOD_LIB) $(DNAFUNCS_LIB) $(F2C_LIB) $(MATH_LIB) $(LAPACK_LIB)  $(F2C_LIB)


%.o: %.cpp
	$(CC)  $(OPTIONS) $(OPTIM_SPEED) -c $<


%.o: %.c
	$(CC)  $(OPTIONS) $(OPTIM_SPEED) -c $<





