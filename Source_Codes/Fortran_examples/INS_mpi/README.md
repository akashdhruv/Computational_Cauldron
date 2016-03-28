
## Fortran based Lid Dirven Cavity Solver using make utility and MPI 

Presently in Debug mode

Execution instructions

  1. Download the source code 
  2. Make sure you have make utility and the latest version of GNU installed
  3. Type

     ~~~terminal 
        make
        mpirun -n [number_of_procs] ./Solver 
     ~~~

  4. To plot results, edit the python file (plot.py) to match your grid size and simply type

     ~~~terminal
        python plot.py
     ~~~ 

Author - Akash Dhruv

