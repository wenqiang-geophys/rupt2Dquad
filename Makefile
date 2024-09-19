
BINDIR := ./bin
SRCDIR := ./src
OBJDIR := ./obj

EXE := exe_solver

# Order of accuracy O = pOrder+1
pOrder := 4

# compilation on Sherlock, comment it to test on your laptop
# module load netcdf-fortran
# module load flexiblas
Sherlock := ON

BLAS := /opt/homebrew/opt/flexiblas
NETCDF := /opt/homebrew/opt/netcdf
NETCDFF := /opt/homebrew/opt/netcdf-fortran

ifeq "$(Sherlock)" "ON"
BLAS := /share/software/user/open/flexiblas/3.1.3
NETCDF := /share/software/user/open/netcdf-c/4.9.0
NETCDFF := /share/software/user/open/netcdf-fortran/4.5.4
endif

FC := mpif90 -O3 -cpp -DpOrder=$(pOrder) -DPSV -DSTRAIN # -DVERSION1  # -DSYM # -DFD
LDFLAGS := -lnetcdf -L${NETCDFF}/lib -lnetcdff
LDFLAGS := -L${BLAS}/lib -lflexiblas -L${NETCDF}/lib -lnetcdf -L${NETCDFF}/lib -lnetcdff
INC := -I${BLAS}/include -I$(NETCDFF)/include
#LDFLAGS := -mkl

SRC := yaml_types.f90       \
       yaml.f90             \
       yaml_settings.f90    \
       mod_para.f90         \
       mod_string.f90       \
       mod_read.f90         \
       mod_gll.f90          \
       mod_fd.f90           \
       mod_types.f90        \
       mod_check.f90        \
       mod_mesh.f90         \
       mod_geometry.f90     \
       mod_funcs.f90        \
       mod_numflux_pml.f90  \
       mod_numflux_v2.f90   \
       mod_wave.f90         \
       mod_init_fault.f90   \
       mod_fault.f90        \
       mod_mpi.f90          \
       mod_exchange.f90     \
       mod_damp.f90         \
       mod_pml.f90          \
       mod_smooth.f90       \
       mod_source.f90       \
       mod_io_fault.f90     \
       mod_io_free.f90      \
       mod_io_inter_s.f90   \
       mod_io_inter_f.f90   \
       mod_io_wave.f90      \
       main.f90

OBJS := $(foreach file,$(SRC),$(OBJDIR)/$(file:.f90=.o))
SRC := $(foreach file,$(SRC),$(SRCDIR)/$(file))

# SUFFIXES rules
.SUFFIXES:
.SUFFIXES: .f90 .o
.SUFFIXES: .F90 .o

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
