subroutine Grid_init()
           
             use Grid_data

#include "Solver.h"

            implicit none

            integer :: I

#ifdef SOLVER_GRID_UG
 
             Lx=1.0
             Ly=1.0

             iProcs = 1
             jProcs = 1

             Nx = Nxb * iProcs
             Ny = Nxb * jProcs

             dx=Lx/Nx
             dy=Ly/Ny

             inRe = .01

             do i=1,Nyb+1
                 x(:,i)=dx*(/(I,I=0,Nxb)/)
             enddo

             do i=1,Nxb+1
                 y(i,:)=dy*(/(I,I=0,Nyb)/)
             enddo

             t = 20.0

             !dt = 0.5*min(dx,dy)
             dt = 0.5*(dx**2)*(dy**2)/(inRe*((dx**2)+(dy**2)))

             nt = t/dt

#endif

#ifdef SOLVER_GRID_NON_UG
#endif

#ifdef POISSON_SOLVER_JACOBI
             omega = 0
#endif

#ifdef POISSON_SOLVER_GS
             omega = 0
#endif

#ifdef POISSON_SOLVER_GSOR

            omega = 1.1
#endif

end subroutine Grid_init
