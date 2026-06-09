HOME=.
LIB=../lib

EXEF=average_function.x
TRGT=average_function

# ==========================================================
# Source codes
# ==========================================================

OBJ=$(LIB)/util.f90      \
    $(LIB)/const.f90     \
    ${LIB}/mtwist.f90    \
    ${LIB}/random.f90    \
    $(LIB)/input.f90     \
    $(LIB)/output.f90    \
    $(LIB)/cv.f90        \
    $(LIB)/bootstrap.f90 \
    $(HOME)/ctrl.f90     \
    $(HOME)/analyze.f90  \
    $(HOME)/main.f90

# Select Fortran compiler
#
chk_ifort := $(shell which ifort    >/dev/null 2>&1 && echo 1 || echo 0)
chk_ifx   := $(shell which ifx      >/dev/null 2>&1 && echo 1 || echo 0)
chk_gfort := $(shell which gfortran >/dev/null 2>&1 && echo 1 || echo 0)

# ==========================================================
# Select Math Library
#
# Usage:
#
#   make
#
#       MKLROOT                 -> MKL
#       otherwise OPENBLAS_ROOT -> OpenBLAS
#
#   make MATHLIB=openblas \
#        OPENBLAS_ROOT=/path/to/openblas
#
#   make MATHLIB=mkl
#
# ==========================================================

#ifndef MATHLIB
#    ifdef MKLROOT
#        MATHLIB = mkl
#    else ifdef OPENBLAS_ROOT
#        MATHLIB = openblas
#    else
#        $(error Neither OPENBLAS_ROOT nor MKLROOT is defined)
#    endif
#endif
#
#ifeq ($(MATHLIB),openblas)
#
#    ifndef OPENBLAS_ROOT
#        $(error OPENBLAS_ROOT is not defined)
#    endif
#
#    FFLAGS += -I$(OPENBLAS_ROOT)/include
#
#    LDLIBS += -L$(OPENBLAS_ROOT)/lib \
#              -Wl,-rpath,$(OPENBLAS_ROOT)/lib \
#              -lopenblas
#
#else ifeq ($(MATHLIB),mkl)
#
#    ifndef MKLROOT
#        $(error MKLROOT is not defined)
#    endif
#
#    LDLIBS += -L$(MKLROOT)/lib/intel64 \
#              -Wl,-rpath,$(MKLROOT)/lib/intel64 \
#              -lmkl_rt
#
#else
#
#    $(error Unknown MATHLIB=$(MATHLIB))
#
#endif
#
#$(info Math library     : $(MATHLIB))

# ==========================================================
# Select Compiler
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
                FFLAGS += -qopenmp

        else ifeq ($(chk_ifx),1)
                FC = ifx
                FFLAGS += -qopenmp

        else ifeq ($(chk_gfort),1)
                FC = gfortran
                FFLAGS += -fopenmp

        else
                $(error No Fortran compiler found)
        endif
endif

$(info Fortran compiler : $(FC))

# ==========================================================
# Build
# ==========================================================

$(TRGT) : $(OBJ)
	$(FC) $(FFLAGS) -o $(EXEF) $(OBJ) $(LDLIBS)
	mkdir -p ../bin
	cp $(EXEF) ../bin/
	rm -f *.mod

