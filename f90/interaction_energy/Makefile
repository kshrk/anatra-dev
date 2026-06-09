HOME=.
LIB=../lib

EXEF=interaction_energy.x
TRGT=interaction_energy

# ==========================================================
# Source codes
# ==========================================================

OBJ=$(XDR_ROOT)/xdr.F90     \
    $(LIB)/util.f90         \
    $(LIB)/const.f90        \
    $(LIB)/grid3d.f90       \
    $(LIB)/input.f90        \
    $(LIB)/output.f90       \
    $(LIB)/netcdf.f90       \
    $(LIB)/dcd.f90          \
    $(LIB)/xtclib.f90       \
    $(LIB)/xyz.f90          \
    $(LIB)/traj.f90         \
    $(LIB)/com.f90          \
    $(LIB)/prmtop.f90       \
    $(LIB)/anaparm.f90      \
    $(LIB)/pot.f90          \
    $(HOME)/ctrl.f90        \
    $(HOME)/setup.f90       \
    $(HOME)/pme_str.f90     \
    $(HOME)/pme.f90         \
    $(HOME)/analyze.f90     \
    $(HOME)/main.f90

# ==========================================================
# External libraries
# ==========================================================

XDR_ROOT      := $(realpath ../lib/external/xdr-interface-fortran)
XDR_LIBDIR    := $(XDR_ROOT)/xdrfile-1.1.4/lib

NETCDF_ROOT   := $(realpath ../lib/external/netcdf/netcdf)
NETCDF_LIBDIR := $(NETCDF_ROOT)/lib
NETCDF_INCDIR := $(NETCDF_ROOT)/include

# ==========================================================
# Select Fortran compiler
# ==========================================================

chk_ifort := $(shell which ifort    >/dev/null 2>&1 && echo 1 || echo 0)
chk_ifx   := $(shell which ifx      >/dev/null 2>&1 && echo 1 || echo 0)
chk_gfort := $(shell which gfortran >/dev/null 2>&1 && echo 1 || echo 0)

# ==========================================================
# Select Math Library
#
# Priority:
#   MKLROOT
#   OPENBLAS_ROOT
# ==========================================================

ifndef MATHLIB
    ifdef MKLROOT
        MATHLIB = mkl
    else ifdef OPENBLAS_ROOT
        MATHLIB = openblas
    else
        $(error Neither OPENBLAS_ROOT nor MKLROOT is defined)
    endif
endif

ifeq ($(MATHLIB),openblas)

    ifndef OPENBLAS_ROOT
        $(error OPENBLAS_ROOT is not defined)
    endif

    FFLAGS += -I$(OPENBLAS_ROOT)/include

    LDLIBS += -L$(OPENBLAS_ROOT)/lib \
              -Wl,-rpath,$(OPENBLAS_ROOT)/lib \
              -lopenblas

else ifeq ($(MATHLIB),mkl)

    ifndef MKLROOT
        $(error MKLROOT is not defined)
    endif

    LDLIBS += -L$(MKLROOT)/lib/intel64 \
              -Wl,-rpath,$(MKLROOT)/lib/intel64 \
              -lmkl_rt

else

    $(error Unknown MATHLIB=$(MATHLIB))

endif

$(info Math library     : $(MATHLIB))

# ==========================================================
# Select FFT code 
#
# Priority:
#   MKL
#   FFTE 
# ==========================================================

ifndef FFTLIB
    ifdef MKLROOT
	FFTLIB = mkl
    else ifdef FFTE
	FFTLIB = ffte
    endif
endif

ifeq ($(FFTLIB),mkl)

    ifndef MKLROOT
        $(error MKLROOT is not defined)
    endif

    ifeq ($(MATHLIB), openblas)
        $(error MKL-FFT can not be used when MATH library is OpenBLAS)
    endif

    OBJFFT = $(LIB)/fftmkl.f90

else ifeq ($(FFTLIB),ffte)

    OBJFFT = $(LIB)/ffte/factor.f   \
             $(LIB)/ffte/kernel.f   \
             $(LIB)/ffte/fft235.f   \
             $(LIB)/ffte/zdfft3d.f  \
             $(LIB)/ffte/dzfft3d.f  \
             $(LIB)/ffte.f90

else

    $(error Unknown FFTLIB=$(FFTLIB))

endif

$(info FFT library      : $(FFTLIB))


# ==========================================================
# NetCDF
# ==========================================================

FFLAGS += -I$(NETCDF_INCDIR)

LDLIBS += -L$(NETCDF_LIBDIR) \
          -Wl,-rpath,$(NETCDF_LIBDIR) \
          -lnetcdff \
          -lnetcdf

# ==========================================================
# XDRFILE
# ==========================================================

LDLIBS += -L$(XDR_LIBDIR) \
          -lxdrfile

# ==========================================================
# Select compiler
# ==========================================================

ifeq ($(FC),ifort)
        FFLAGS += -qopenmp

else ifeq ($(FC),ifx)
        FFLAGS += -qopenmp

else ifeq ($(FC),gfortran)
        FFLAGS += -fopenmp

else
        ifeq ($(chk_ifort),1)
                FC = ifort
                #FFLAGS += -qopenmp

        else ifeq ($(chk_ifx),1)
                FC = ifx
                #FFLAGS += -qopenmp

        else ifeq ($(chk_gfort),1)
                FC = gfortran
                #FFLAGS += -fopenmp

        else
                $(error No Fortran compiler found)
        endif
endif

$(info Fortran compiler : $(FC))

# ==========================================================
# Build
# ==========================================================

$(TRGT) : $(OBJ)
	$(FC) $(FFLAGS) -o $(EXEF) $(OBJFFT) $(OBJ) $(LDLIBS)
	mkdir -p ../bin
	cp $(EXEF) ../bin/
	rm -f *.mod
