Solver : Solver.o Grid_init.o IO_display.o Grid_data.o
				 gfortran -o Solver Solver.o Grid_init.o IO_display.o Grid_data.o

Grid_interface.o Grid_interface.mod: Grid_interface.F90
				gfortran -c Grid_interface.F90

Grid_data.o Grid_data.mod: Grid_data.F90 Solver.h
				gfortran -c Grid_data.F90

Grid_init.o : Grid_init.F90 Grid_data.mod Solver.h
				gfortran -c Grid_init.F90

IO_interface.o IO_interface.mod: IO_interface.F90
				gfortran -c IO_interface.F90

IO_display.o: IO_display.F90 Grid_data.mod
				gfortran -c IO_display.F90

RuntimeParameters_read.o : RuntimeParameters_read.F90
				gfortran -c RuntimeParameters_read.F90

Solver.o : Solver.F90 Grid_interface.mod IO_interface.mod
				gfortran -c Solver.F90


clean:
	rm *.o *.mod
