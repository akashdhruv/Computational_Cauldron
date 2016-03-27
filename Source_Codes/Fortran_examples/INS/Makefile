Solver : Solver.o Grid_init.o IO_display.o Grid_data.o IncompNS_init.o IncompNS_data.o IncompNS_rk3.o
				 gfortran -o Solver Solver.o Grid_init.o IO_display.o Grid_data.o IncompNS_init.o IncompNS_data.o IncompNS_rk3.o

Grid_interface.o Grid_interface.mod: Grid_interface.F90
				gfortran -c Grid_interface.F90

Grid_data.o Grid_data.mod: Grid_data.F90 Solver.h
				gfortran -c Grid_data.F90

Grid_init.o : Grid_init.F90 Grid_data.mod Solver.h
				gfortran -c Grid_init.F90

IO_interface.o IO_interface.mod: IO_interface.F90
				gfortran -c IO_interface.F90

IO_display.o: IO_display.F90 Grid_data.mod IncompNS_data.mod
				gfortran -c IO_display.F90

IncompNS_interface.o IncompNS_interface.mod: IncompNS_interface.F90
				gfortran -c IncompNS_interface.F90

IncompNS_data.o IncompNS_data.mod: IncompNS_data.F90 Solver.h
				gfortran -c IncompNS_data.F90 

IncompNS_init.o: IncompNS_init.F90 Grid_data.mod IncompNS_data.mod IncompNS_interface.mod Solver.h
				gfortran -c IncompNS_init.F90

IncompNS_rk3.o: IncompNS_rk3.F90 Grid_data.mod IncompNS_data.mod Solver.h
				gfortran -c IncompNS_rk3.F90

Solver.o : Solver.F90 Grid_interface.mod IO_interface.mod IncompNS_interface.mod
				gfortran -c Solver.F90


clean:
	rm *.o *.mod
