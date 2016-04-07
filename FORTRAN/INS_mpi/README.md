## Fortran based Lid Dirven Cavity Solver for an interconnected mesh network using make utility and MPI 

Execution instructions

  1. Download the source code 
  2. Make sure you have make utility and the latest version of GNU and MPI installed
  3. Edit the Makefile to include your MPI path.

     ~~~terminal 
        make
        mpirun -n [number_of_procs] ./Solver 
     ~~~

  5. Note that the number of processors in each direction is HK defined in the header file Solver.h. And for 2D cases HD = 2.      Therefore, the total number of processors will be HK^HD.

  4. To plot results, edit the python file (plot.py) to match your grid size and simply type (make sure k = HK and d = HD)

     ~~~terminal
        python plot.py
     ~~~ 

Author - Akash Dhruv

