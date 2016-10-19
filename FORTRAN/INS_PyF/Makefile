FF=f2py
FFLAGS = -c --fcompiler=gfortran -m

all: INS.so POISSON.so HEAT.so

POISSON.so: POISSON.F90
	$(FF) $(FFLAGS) POISSON POISSON.F90

INS.so: INS.F90
	$(FF) $(FFLAGS) INS INS.F90

HEAT.so: HEAT.F90
	$(FF) $(FFLAGS) HEAT HEAT.F90

cleanall:
	rm -f *.so *.dat

clean:
	rm -f *.dat

run: Solver.py
	python Solver.py
