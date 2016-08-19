
#FF=gfortran
#FFLAGS = -fopenmp

FF=pgfortran
FFLAGS=-Mcuda -mp
#FFLAGS = -mp 

Solver : Solver.o Grid_init.o IO_display.o Grid_data.o IncompNS_init.o IncompNS_data.o IncompNS_rk3.o Poisson_solver.o IO_write.o CUDA_poisson.o
				$(FF) $(FFLAGS) -o Solver Solver.o Grid_init.o IO_display.o Grid_data.o IncompNS_init.o IncompNS_data.o IncompNS_rk3.o Poisson_solver.o IO_write.o CUDA_poisson.o

Grid_interface.o Grid_interface.mod: Grid_interface.F90
				$(FF) $(FFLAGS) -c Grid_interface.F90

Grid_data.o Grid_data.mod: Grid_data.F90 Solver.h
				$(FF) $(FFLAGS) -c Grid_data.F90

Grid_init.o : Grid_init.F90 Grid_data.mod Solver.h
				$(FF) $(FFLAGS) -c Grid_init.F90

IO_interface.o IO_interface.mod: IO_interface.F90 Solver.h
				$(FF) $(FFLAGS) -c IO_interface.F90

IO_display.o: IO_display.F90
				$(FF) $(FFLAGS) -c IO_display.F90

IO_write.o: IO_write.F90 Solver.h
				$(FF) $(FFLAGS) -c IO_write.F90

IncompNS_interface.o IncompNS_interface.mod: IncompNS_interface.F90
				$(FF) $(FFLAGS) -c IncompNS_interface.F90

IncompNS_data.o IncompNS_data.mod: IncompNS_data.F90 Solver.h
				$(FF) $(FFLAGS) -c IncompNS_data.F90 

IncompNS_init.o: IncompNS_init.F90 IncompNS_data.mod IncompNS_interface.mod Grid_data.mod Solver.h
				$(FF) $(FFLAGS) -c IncompNS_init.F90

IncompNS_rk3.o: IncompNS_rk3.F90 Grid_data.mod IncompNS_data.mod IO_interface.mod Poisson_interface.mod Solver.h
				$(FF) $(FFLAGS) -c IncompNS_rk3.F90

CUDA_interface.o CUDA_interface.mod: CUDA_interface.F90
				$(FF) $(FFLAGS) -c CUDA_interface.F90

CUDA_poisson.o: CUDA_poisson.F90
				$(FF) $(FFLAGS) -c CUDA_poisson.F90

Poisson_interface.o Poisson_interface.mod: Poisson_interface.F90 Solver.h
				$(FF) $(FFLAGS) -c Poisson_interface.F90

Poisson_solver.o: Poisson_solver.F90 IncompNS_data.mod Grid_data.mod Solver.h CUDA_interface.mod
				$(FF) $(FFLAGS) -c Poisson_solver.F90

Solver.o : Solver.F90 Grid_interface.mod IO_interface.mod IncompNS_interface.mod
				$(FF) $(FFLAGS) -c Solver.F90


clean:
	rm *.o *.mod *.dat
