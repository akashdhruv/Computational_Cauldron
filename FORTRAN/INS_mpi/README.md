## Fortran based Lid Dirven Cavity Solver for an interconnected mesh network using make utility and MPI 
## Needs Debugging
Execution instructions

  1. Download the source code 
  2. Make sure you have make utility and the latest version of GNU installed
  3. Type

     ~~~terminal 
        make FF=[your_fortran_compiler] MPIFF=[your_mpif90_compiler]
        mpirun -n [number_of_procs] ./Solver 
     ~~~

  4. To plot results, edit the python file (plot.py) to match your grid size and simply type

     ~~~terminal
        python plot.py
     ~~~ 

Author - Akash Dhruv

