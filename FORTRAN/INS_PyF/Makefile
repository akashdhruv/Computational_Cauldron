FF=f2py
FFLAGS = -c -m

all: INS.so POISSON.so HEAT.so

POISSON.so: POISSON.F90
	$(FF) $(FFLAGS) POISSON POISSON.F90

INS.so: INS.F90
	$(FF) $(FFLAGS) INS INS.F90

HEAT.so: HEAT.F90
	$(FF) $(FFLAGS) HEAT HEAT.F90

clean:
	rm -f *.so

run: Solver.py
	python Solver.py
