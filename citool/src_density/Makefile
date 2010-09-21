TARGET = ../bin/density4CItool.x

#################################################################
# keys:  <default>, citest
#  
# for:  density4CItool
# created:  17 sep 2010 (in Castello)
# machine:  gel.act.uji.es (linux cluster)
# compiler:  ifort 11.1(059)
# libraries:  
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
LLIBS =
#  ./lib/libfftw3.0.1.a  # MAH..... fatta col pgi

# LLIBS = -L ./lib -larpack_leela_fc8-ifort64v110_mkl              \
#         -I /opt/intel/Compiler/11.0/074/mkl/include  \
#         -L /opt/intel/Compiler/11.0/074/mkl/lib/em64t  \
#         -lmkl_solver_lp64 -lmkl_intel_lp64 -lmkl_intel_thread  \
#         -lmkl_core -liomp5 -lpthread -lm

#################################################################
#      objects files
#################################################################
VPATH = ../src/
##FPATH= ../src/

OBJS =
main.o                     \
mod_staticdata.o           \
mod_indata.o               \
mod_logga.o                \
mod_myaux.o                \
mod_indexx.o               \
mod_inoutrs.o              \
mod_denscalc.o

#################################################################
#      rules 
#################################################################

$(TARGET) : $(OBJS)
	$(LNK) -I $(VPATH) -o $(TARGET) $(LFLAGS) $^ $(LLIBS)


%.o : %.f90
	$(F90) -I $(VPATH) $(F90FLAGS) -c $< -o $@

%.o : %.f
	$(F77) $(F77FLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)
	rm -f *.mod
	rm -f $(TARGET)


#################################################################
#      dependencies
#################################################################

main.o :	mod_indata.o               \
			mod_logga.o                \
			mod_myaux.o                \
			mod_indexx.o               \
			mod_inoutrs.o              \
			mod_denscalc.o
			
mod_staticdata.o: 

mod_indata.o:	mod_staticdata.o              \
				mod_myaux.o

mod_logga.o:	mod_indata.o

mod_myaux.o:

mod_indexx.o:

mod_inoutrs.o:

mod_denscalc.o:
