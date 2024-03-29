TARGET = ../bin/citool.x
TARGETCITEST = ../bin/generatecitest.x

#################################################################
# keys:  <default>, citest
#  
# for:  CItool
# created:  6 oct 2010 
# machine:  leela.fisica.unimo.it (linux fedora 8)
# compiler:  ifort 11.0(074)
# libraries:  ARPACK (linked with mkl blas/lapack)
#                                            by Andrea Bertoni
#################################################################

## the command shell
SHELL = /bin/csh

## name of the fortran 90, fortran 77, C compilers
#F90 = pgf90 -Mfreeform -r8 -byteswapio
F90 = ifort
F77 = ifort
CC =

## flags to use for fortran 90, fortran 77, C compilation
# NOTE: with -f77rtl the first column in a namelist is skipped
#       i.e. the & must be in the second one
F90FLAGS = -w -CB -C #-O0  # -static -f77rtl -convert big_endian
F77FLAGS = -w -CB -C #-O0  # -static -f77rtl -convert big_endian
CFLAGS =

## name of the program to link the program units
LNK = $(F90)

## flags to use at link time
LFLAGS = $(F90LAGS)

## LLIBS are libraries to use at link time
#LLIBS = -L /home/commenda/lib -larpack_gel-MKL              \
#        -I /opt/intel/Compiler/11.1/059/mkl/include  \
#        -L /opt/intel/Compiler/11.1/059/mkl/lib/32   \
#        -lmkl_intel -lmkl_intel_thread  \
#        -lmkl_core -liomp5 -lpthread -lm
LLIBS = -L ./lib -larpack_leela_fc8-ifort64v110_mkl              \
        -I /opt/intel/Compiler/11.0/074/mkl/include  \
        -L /opt/intel/Compiler/11.0/074/mkl/lib/em64t  \
        -lmkl_intel_lp64 -lmkl_intel_thread  \
        -lmkl_core -liomp5 -lpthread -lm

#################################################################
#      objects files
#################################################################

OBJS =  main.o                     \
        mod_staticdata.o           \
        mod_indata.o               \
        mod_logga.o                \
        mod_myaux.o                \
        mod_specialf.o             \
        mod_cimakeindex.o          \
        mod_createhspace.o         \
        mod_composehspace.o        \
        mod_blockizehamiltonian.o  \
        mod_createhamiltonian.o    \
        mod_diagonalize.o          \
        mod_indexx.o               \
        mod_writempstates.o

#################################################################
#      rules 
#################################################################

$(TARGET) : $(OBJS)
	$(LNK) -o $(TARGET) $(LFLAGS) $^ $(LLIBS)


%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

%.o : %.f
	$(F77) $(F77FLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)
	rm -f *.mod
	rm -f $(TARGET)

citest: generatecitest.o mod_staticdata.o mod_indata.o
	$(LNK) -o $(TARGETCITEST) $(LFLAGS) $^

#################################################################
#      dependencies
#################################################################

main.o:		mod_indata.o                  \
			mod_logga.o                   \
			mod_myaux.o                   \
			mod_specialf.o                \
			mod_cimakeindex.o             \
			mod_createhspace.o            \
			mod_composehspace.o           \
			mod_blockizehamiltonian.o     \
			mod_createhamiltonian.o       \
			mod_diagonalize.o             \
			mod_writempstates.o

mod_staticdata.o: 

mod_indata.o:	mod_staticdata.o              \
				mod_myaux.o

mod_logga.o:		mod_indata.o

mod_myaux.o:

mod_specialf.o:

mod_cimakeindex.o:	mod_indata.o

mod_createhspace.o:	mod_specialf.o

mod_composehspace.o:

mod_blockizehamiltonian.o:	mod_indata.o         \
							mod_logga.o

mod_createhamiltonian.o:	mod_indata.o

mod_diagonalize.o:	mod_indexx.o

mod_indexx.o:

mod_writempstates.o:	mod_indata.o


# citest

generatecitest.o:	mod_indata.o

