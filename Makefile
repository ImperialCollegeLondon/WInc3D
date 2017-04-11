#=======================================================================
# Makefile for Imcompact3D
#=======================================================================

# Choose pre-processing options
#   -DSHM	   - enable shared-memory implementation
#   -DDOUBLE_PREC  - use double-precision
OPTIONS = -DDOUBLE_PREC

# Choose an FFT engine, available options are:
#   fftw3      - FFTW version 3.x
#   generic    - A general FFT algorithm (no 3rd-party library needed)
FFT= generic

# Paths to FFTW 3
FFTW3_PATH=   # full path of FFTW installation if using fftw3 engine above
FFTW3_INCLUDE = -I$(FFTW3_PATH)/include
FFTW3_LIB = -L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f

# Specify Fortran and C compiler names and flags here
# Normally, use MPI wrappers rather than compilers themselves 
# Supply a Fortran pre-processing flag together with optimisation level flags
# Some examples are given below:

#FC =  
#OPTFC = 
#CC = 
#CFLAGS = 

# PGI
#FC = ftn
#OPTFC = -fast -O3 -Mpreprocess
#CC = cc
#CFLAGS = -O3

# PathScale
#FC = ftn
#OPTFC = -Ofast -cpp
#CC = cc
#CFLAGS = -O3

# GNU
FC = mpif90
OPTFC = -O3 -funroll-loops -ftree-vectorize -fcray-pointer -cpp -ffree-line-length-0
CC = mpicc
CFLAGS = -O3
LIBS = -llapack -lblas 

# Cray
#FC = ftn
#OPTFC = -e Fm
#CC = cc
#CFLAGS = 

#-----------------------------------------------------------------------
# Normally no need to change anything below

# include PATH 
ifeq ($(FFT),generic)
  INC=
else ifeq ($(FFT),fftw3)
  INC=
endif

# library path
ifeq ($(FFT),generic)
   LIBFFT= $(SPUD_LIB)
else ifeq ($(FFT),fftw3)
   LIBFFT=$(FFTW3_LIB) $(SPUD_LIB)
endif

SRC = decomp_2d.f90 Actuator_Line_Model_Utils.f90 Airfoils.f90 dynstall_legacy.f90 Actuator_Line_Element.f90 Actuator_Line_Turbine.f90 Actuator_Line_Write_Output.f90 Actuator_Line_Model.f90 Actuator_Line_Source.f90 glassman.f90 fft_$(FFT).f90 module_param.f90 io.f90 variables.f90 poisson.f90 schemes.f90 convdiff.f90 incompact3d.f90 navier.f90 filter.f90 derive.f90 parameters.f90 tools.f90 visu.f90 tecplot.f90

#-----------------------------------------------------------------------
# Normally no need to change anything below

ifneq (,$(findstring DSHM,$(OPTIONS)))
SRC := FreeIPC.f90 $(SRC)  
OBJ =	$(SRC:.f90=.o) alloc_shm.o FreeIPC_c.o
else
OBJ =	$(SRC:.f90=.o)
endif	

all: incompact3d

alloc_shm.o: alloc_shm.c
	$(CC) $(CFLAGS) -c $<

FreeIPC_c.o: FreeIPC_c.c
	$(CC) $(CFLAGS) -c $<

incompact3d : $(OBJ)
	$(FC) -O3 -o $@ $(OBJ) $(LIBFFT) $(LIBS)

%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) $(INC) -c $<
	
.PHONY: clean 
clean:
	rm -f *.o *.mod incompact3d

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#
