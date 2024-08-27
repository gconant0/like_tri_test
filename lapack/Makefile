LIN_SRC_DIR = ./
include $(LIN_SRC_DIR)make.inc


 
ALLAUX =  $(LIN_SRC_DIR)xerbla.o \
	 $(LIN_SRC_DIR)lsame.o \
	 $(LIN_SRC_DIR)idamax.o \
         $(LIN_SRC_DIR)dswap.o \
	 $(LIN_SRC_DIR)dcopy.o \
	 $(LIN_SRC_DIR)dgemm.o \
         $(LIN_SRC_DIR)dtrmm.o \
	 $(LIN_SRC_DIR)dgemv.o \
         $(LIN_SRC_DIR)dtrmv.o \
	 $(LIN_SRC_DIR)dger.o \
         $(LIN_SRC_DIR)ddot.o \
         $(LIN_SRC_DIR)daxpy.o \
         $(LIN_SRC_DIR)ilaenv.o \
         $(LIN_SRC_DIR)ieeeck.o \
 
  

DZLAUX = \
    $(LIN_SRC_DIR)dlabad.o  $(LIN_SRC_DIR)dlacpy.o \
    $(LIN_SRC_DIR)dlamch.o  \
    $(LIN_SRC_DIR)dscal.o  \
    $(LIN_SRC_DIR)dlascl.o  \
    $(LIN_SRC_DIR)drot.o  \
    $(LIN_SRC_DIR)dnrm2.o  \
    $(LIN_SRC_DIR)dhseqr.o \
    $(LIN_SRC_DIR)dlapy2.o \
    $(LIN_SRC_DIR)dlahrd.o \
    $(LIN_SRC_DIR)dlarfx.o \
    $(LIN_SRC_DIR)dlahqr.o \
    $(LIN_SRC_DIR)dlarf.o \
    $(LIN_SRC_DIR)dlarfg.o \
    $(LIN_SRC_DIR)dlanhs.o \
    $(LIN_SRC_DIR)dlanv2.o \
    $(LIN_SRC_DIR)dorghr.o \
    $(LIN_SRC_DIR)dtrevc.o \
    $(LIN_SRC_DIR)dlaln2.o \
    $(LIN_SRC_DIR)dladiv.o \
    $(LIN_SRC_DIR)dorgqr.o \
    $(LIN_SRC_DIR)dlassq.o \
    $(LIN_SRC_DIR)dorg2r.o \
    $(LIN_SRC_DIR)dlarft.o \
    $(LIN_SRC_DIR)dlaset.o \
    $(LIN_SRC_DIR)dgetri.o \
    $(LIN_SRC_DIR)dtrtri.o \
    $(LIN_SRC_DIR)dtrsm.o \
    $(LIN_SRC_DIR)dtrti2.o \
    $(LIN_SRC_DIR)dgetrf.o \
    $(LIN_SRC_DIR)dgetf2.o \
    $(LIN_SRC_DIR)dlaswp.o \


DLASRC = \
    $(LIN_SRC_DIR)dgebal.o \
    $(LIN_SRC_DIR)dgeev.o  \
    $(LIN_SRC_DIR)dgehrd.o \
    $(LIN_SRC_DIR)dlange.o \
    $(LIN_SRC_DIR)dlartg.o \
    $(LIN_SRC_DIR)dlarfb.o \
    $(LIN_SRC_DIR)dgebak.o \
    $(LIN_SRC_DIR)dgehd2.o \
  




all: double 


double:  $(DLASRC) $(ALLAUX) $(DZLAUX) 
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB)  $(DLASRC) $(ALLAUX) \
	$(DZLAUX) 
	$(RANLIB) $(LAPACKLIB)	

#$(ALLAUX): $(FRC)
#$(DZLAUX): $(FRC)
#$(DLASRC): $(FRC)

#FRC:
#	@FRC=$(FRC)


#%.o: %.f 
#	$(FORTRAN) -O3 -c $(LIN_SRC_DIR)$<

%.o: %.c
	$(cc) -O1 -c $(LIN_SRC_DIR)$<


