#Master Makefile for pair_ks.pl package--makes the associated c++ routines
#G. Conant 10/16/01


default: all


all: ./lapack/liblapack.a ./libf2c/libf2c.a  progs

./lapack/liblapack.a:
	cd lapack; make

./libf2c/libf2c.a:
	cd libf2c; make

progs: ./lapack/liblapack.a ./libf2c/libf2c.a
	cd src; make

clean: 
	rm -f like_tri_test  src/*.o src/*.a lapack/*.o lapack/*.a libf2c/*.o libf2c/*.a
