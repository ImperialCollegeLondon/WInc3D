#=======================================================================
# Makefile for Incompact3D
#=======================================================================
#SET COMPILER AND VERSION BASED ON WHETHER a HPC or local PC is used
HOST=$(shell domainname | sed 's/\.//g')
HOSTCX1=cx1hpcicacuk
HOSTCX2=cx2hpcicacuk

OPTIONS = -DDOUBLE_PREC

ifeq ($(HOST),$(HOSTCX1))
	
FFT = fftw3
FFTW3_INCLUDE = -I$(MKLROOT)/include/fftw
FFTW3_LIB = -mkl -L$(MKLROOT)/interface/fftw3xf -lfftw3xf_intel 

FC = mpiF90
OPTFC = -O3 -xAVX -cpp
CC=cc
CFLAGS= -O3 -xAVX

else
# Choose an FFT engine, available options are:
#   fftw3      - FFTW version 3.x
#   generic    - A general FFT algorithm (no 3rd-party library needed)
FFT= generic

# Paths to FFTW 3
FFTW3_PATH=   # full path of FFTW installation if using fftw3 engine above
FFTW3_INCLUDE = -I$(FFTW3_PATH)/include
FFTW3_LIB = -L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f

# GNU
FC = mpif90
OPTFC = -O3 -funroll-loops -ftree-vectorize -fcray-pointer -cpp -ffree-line-length-0 #-ffpe-trap=invalid,zero
CC = mpicc
CFLAGS = -O3
LIBS = #-llapack -lblas 
DEGUG = #-g -static 
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

endif

SRC = decomp_2d.f90 glassman.f90 fft_$(FFT).f90 module_param.f90 io.f90 variables.f90 poisson.f90 les_models.f90 schemes.f90 convdiff.f90 acl_utils.f90 airfoils.f90 dynstall.f90 acl_elem.f90 acl_turb.f90 acl_out.f90 acl_model.f90 acl_source.f90 incompact3d.f90 navier.f90 filter.f90 derive.f90 parameters.f90 tools.f90 visu.f90 probe.f90 cfl.f90 ABL.f90 

SRCALM = decomp_2d.f90 acl_utils.f90 airfoils.f90 dynstall.f90 acl_elem.f90 acl_turb.f90 acl_out.f90 acl_model.f90 uALM.f90 

ifneq (,$(findstring DSHM,$(OPTIONS)))
SRC := FreeIPC.f90 $(SRC)  
OBJ =	$(SRC:.f90=.o) alloc_shm.o FreeIPC_c.o
OBJALM =	$(SRCALM:.f90=.o) alloc_shm.o FreeIPC_c.o
else
OBJ =	$(SRC:.f90=.o)
OBJALM =	$(SRCALM:.f90=.o) 
endif	

all: incompact3d visualize

alloc_shm.o: alloc_shm.c
	$(CC) $(CFLAGS) -c $<

FreeIPC_c.o: FreeIPC_c.c
	$(CC) $(CFLAGS) -c $<

incompact3d : $(OBJ)
	$(FC) -O3 -o $@ $(OBJ) $(LIBFFT) $(LIBS) $(DEBUG)

uALM : $(OBJALM)
	$(FC) -O3 -o $@ $(OBJALM) $(LIBFFT) $(LIBS) $(DEBUG)


%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) $(INC) $(DEBUG) -c $<
	
visualize :
	mpif90 paraview_incompact3d.f90 -o visualize 
.PHONY: clean 
clean:
	rm -f *.o *.mod incompact3d visualize

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#
