
BINDIR := ./bin
SRCDIR := ./src
OBJDIR := ./obj

EXE := exe_solver

# WEAK or STRONG form, default is STRONG, WEAK only works for PSV+STRESS case
# PSV or SH, must select one
# STRIAN or STRESS form, must select one
# for SH case: SH+STRIAN
BLAS := /opt/homebrew/opt/openblas
NETCDF := /opt/homebrew/opt/netcdf
NETCDFF := /opt/homebrew/opt/netcdf-fortran

FC := mpif90 -O3 -cpp -DpOrder=4 -DPSV -DSTRAIN -DTPV5 # -DVERSION1  # -DSYM # -DFD
#FC := gfortran -Wall -O3 -cpp -DSTRESS -DPSV

#-fopenmp
#LDFLAGS := -lblas -llapack
#LDFLAGS := -L. -lblas -L. -llapack
#LDFLAGS := -Bstatic -L. -lblas -L. -llapack
#LDFLAGS := -L. -lblas -L. -llapack
#LDFLAGS := libblas.a liblapack.a
#LDFLAGS := -Llapack-3.9.1 -llapack -lrefblas
LDFLAGS := -L${BLAS}/lib -lopenblas -L${BLAS}/lib -llapack -L${NETCDF}/lib -lnetcdf -L${NETCDFF}/lib -lnetcdff
INC := -I${BLAS}/include -I$(NETCDFF)/include
#LDFLAGS := -mkl

SRC := mod_para.f90 \
       mod_string.f90 \
       mod_gll.f90 \
       mod_fd.f90 \
       mod_types.f90 \
       mod_mesh.f90 \
       mod_geometry.f90 \
       mod_funcs.f90 \
       mod_numflux_v2.f90 \
       mod_wave.f90 \
       mod_fault.f90 \
       mod_mpi.f90 \
       mod_exchange.f90 \
       mod_damp.f90 \
       mod_pml.f90 \
       mod_smooth.f90 \
       mod_source.f90 \
       mod_io_fault.f90 \
       mod_io_free.f90 \
       mod_io_inter_s.f90 \
       mod_io_inter_f.f90 \
       mod_io_wave.f90 \
       main.f90

OBJS := $(foreach file,$(SRC),$(OBJDIR)/$(file:.f90=.o))
SRC := $(foreach file,$(SRC),$(SRCDIR)/$(file))

# SUFFIXES rules
.SUFFIXES:
.SUFFIXES: .f90 .o

.PHONY: all

all: solver
EXE := $(BINDIR)/$(EXE)

solver: $(EXE)

dir:
	mkdir -p bin obj
# Target
#all: $(EXE)
all: dir solver

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) -o $@ -c $^ -J$(OBJDIR) $(INC)


$(EXE) : $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -rf $(EXE)
	rm -rf $(OBJDIR)
	mkdir $(OBJDIR)
