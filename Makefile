PROG =	dgv0_mpi.a

SRCS =	collision_mod.f90 commvar.f90 DGVblzmPL.f90 dgvtools_mod.f90 \
	distributions_mod.f90 gaussian_mod.f90 gaussquad.f90 miscmpiset.f90 \
	miscset.f90 nrroutines_mod.f90 nrtype.f90 nrutil.f90 readwrite.f90 \
	sf02.f90 spat_oper_mod.f90 time_integr_mod.f90 algama.f

OBJS =	collision_mod.o commvar.o DGVblzmPL.o dgvtools_mod.o \
	distributions_mod.o gaussian_mod.o gaussquad.o miscmpiset.o miscset.o \
	nrroutines_mod.o nrtype.o nrutil.o readwrite.o sf02.o spat_oper_mod.o \
	time_integr_mod.o algama.o

LIBS = 	

CC = cc
CFLAGS = -O
FC = mpiifort
FFLAGS = -O
F90 = mpiifort
F90FLAGS = -O3
LDFLAGS = -s 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

collision_mod.o: commvar.o nrtype.o
commvar.o: nrtype.o
DGVblzmPL.o: collision_mod.o commvar.o dgvtools_mod.o distributions_mod.o \
	miscmpiset.o miscset.o nrtype.o readwrite.o sf02.o time_integr_mod.o
dgvtools_mod.o: commvar.o distributions_mod.o gaussian_mod.o nrtype.o
distributions_mod.o: nrtype.o
gaussian_mod.o: nrtype.o nrutil.o
miscmpiset.o: commvar.o dgvtools_mod.o miscset.o nrtype.o readwrite.o
miscset.o: commvar.o gaussian_mod.o gaussquad.o nrtype.o
nrroutines_mod.o: 
nrutil.o: nrtype.o
readwrite.o: commvar.o dgvtools_mod.o nrtype.o
sf02.o: commvar.o nrtype.o
spat_oper_mod.o: collision_mod.o commvar.o dgvtools_mod.o distributions_mod.o \
	nrtype.o readwrite.o
time_integr_mod.o: commvar.o nrtype.o spat_oper_mod.o
