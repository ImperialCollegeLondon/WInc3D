#=======================================================================
# Makefile for incompact3D
#=======================================================================

OPTIONS = -DDOUBLE_PREC

# Choose an FFT engine, available options are:
#   fftw3      - FFTW version 3.x
#   generic    - A general FFT algorithm (no 3rd-party library needed)
FFT= generic

# Paths to xbeam 
xbeam_PATH=../xbeam
xbeam_INCLUDE=-I$(xbeam_PATH)/src -I$(xbeam_PATH)/src/xbeam_base
xbeam_LIB=-L$(xbeam_PATH)/lib -lxbeam 

# Paths to FFTW 3
FFTW3_PATH=   # full path of FFTW installation if using fftw3 engine above
FFTW3_INCLUDE = -I$(FFTW3_PATH)/include
FFTW3_LIB = -L$(FFTW3_PATH)/lib -lfftw3 -lfftw3f

# GNU
FC = mpif90
OPTFC = -O3 -funroll-loops -ftree-vectorize -fcray-pointer -cpp -ffree-line-length-0 -g -fbacktrace -ffpe-trap=invalid,zero
CC = mpicc
CFLAGS = -O3 
LIBS = -fopenmp -llapack -lblas  

# include PATH 
ifeq ($(FFT),generic)
  INC=
else ifeq ($(FFT),fftw3)
  INC=
endif

SRC = decomp_2d.f90 glassman.f90 fft_$(FFT).f90 module_param.f90 io.f90 variables.f90 poisson.f90 les_models.f90 SVV.f90 schemes.f90 convdiff.f90 acl_utils.f90 airfoils.f90 dynstall_legacy.f90 dynstall.f90 acl_elem.f90 acl_beam.f90 acl_controller.f90 acl_turb.f90 acl_out.f90 acl_model.f90 acl_source.f90 adm.f90 incompact3d.f90 navier.f90 filters.f90 derive.f90 parameters.f90 tools.f90 visu.f90 probe.f90 cfl.f90 ABL.f90 

ifneq (,$(findstring DSHM,$(OPTIONS)))
SRC := FreeIPC.f90 $(SRC) $(SRCALM)
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
	$(FC) -O3 -o $@ $(OBJ) $(LIBFFT) $(LIBS) $(DEBUG) $(xbeam_LIB) 

%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) $(INC) $(xbeam_INCLUDE) $(DEBUG) $(LIBS) -c $<
	
visualize :
	mpif90 paraview_incompact3d.f90 -o visualize 
.PHONY: clean 
clean:
	rm -f *.o *.mod incompact3d visualize

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#
