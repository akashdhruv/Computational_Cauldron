subroutine Grid_init()
           
             use Grid_data

#include "Solver.h"

            implicit none

            integer :: I

#ifdef SOLVER_GRID_UG
 
             iProcs = 1
             jProcs = 1  

             Lx = 1/iProcs
             Ly = 1/jProcs

             dx=Lx/Nxb
             dy=Ly/Nyb

             inRe = .0003125

             do i=1,Nyb+1
                 x(:,i)=dx*(/(I,I=0,Nxb)/)
             enddo

             do i=1,Nxb+1
                 y(i,:)=dy*(/(I,I=0,Nyb)/)
             enddo

             t = 30.0

             !dt = 0.000001
             dt = .05*min(dx,dy)
             !dt = (0.5*(dx**2)*(dy**2))/(inRe*((dx**2)+(dy**2)))

             nt = t/dt

#endif

#ifdef SOLVER_GRID_NON_UG
#endif

#ifdef POISSON_SOLVER_JACOBI
#endif

#ifdef POISSON_SOLVER_GS
#endif

#ifdef POISSON_SOLVER_GSOR
#endif

end subroutine Grid_init
