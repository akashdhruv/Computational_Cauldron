subroutine Grid_init()
           
             use Grid_data

#include "Solver.h"

            implicit none

#ifdef SOLVER_GRID_UG
 
             Lx=1.0
             Ly=1.0

             Nxb = 20
             Nyb = 20

             iProcs = 1
             jProcs = 1

             Nx = Nxb * iProcs
             Ny = Nxb * jProcs

             dx=Lx/Nx
             dy=Ly/Ny

             inRe = .01
             omega = 1.1

#endif

#ifdef SOLVER_GRID_NON_UG

             Lx=2.0

#endif

#ifdef POISSON_SOLVER_JACOBI
             
             Ly=2.0

#endif

#ifdef POISSON_SOLVER_GS
             Ly=3.0

#endif

#ifdef POISSON_SOLVER_GSOR

             Ly=4.0
#endif

end subroutine Grid_init
