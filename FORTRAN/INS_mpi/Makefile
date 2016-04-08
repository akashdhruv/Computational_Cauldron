MPI_PATH = /usr/local
#MPI_PATH = /usr
#MPI_PATH = /opt/mpich2-1.4.1p1/
#PGI_PATH = /opt/pgi/osx86-64/2016

MPIFF = $(MPI_PATH)/bin/mpif90

FF = gfortran
#FF = $(PGI_PATH)/bin/pgf90

FFLAGS = -c
#FFLAGS = -c -Mmpi=mpich

EXEFLAGS = -o
#EXEFLAGS = -Mmpi=mpich -o

EXE = Solver

%.o %.mod: %.F90
	$(MPIFF) $(FFLAGS) $<

Headers += Solver.h

Modules += Grid_interface.mod Grid_data.mod IO_interface.mod IncompNS_interface.mod IncompNS_data.mod Poisson_interface.mod MPI_interface.mod MPI_data.mod

Objects += Grid_data.o Grid_init.o IncompNS_init.o IncompNS_data.o IncompNS_rk3.o Poisson_solver.o IO_write.o IO_display.o MPIsolver_init.o MPIsolver_finalize.o \
          MPI_data.o MPI_applyBC.o MPI_physicalBC_vel.o MPI_physicalBC_pres.o MPI_CollectResiduals.o Solver.o

#MPI_objects += MPIsolver_init.o MPIsolver_finalize.o MPI_applyBC.o MPI_physicalBC_vel.o MPI_physicalBC_pres.o MPI_CollectResiduals.o

ALL_OBJS = $(Headers) $(Modules) $(Objects)

REBUILDABLES = $(filter-out $(Headers),$(ALL_OBJS)) $(EXE)

#Fortran_objects += $(filter-out $(MPI_objects),$(Objects))

$(EXE): $(ALL_OBJS) 
	$(MPIFF) $(EXEFLAGS) $(EXE) $(Objects)

include Makefile.Depend

clean:
	rm *.dat *.mod *.o $(EXE)
