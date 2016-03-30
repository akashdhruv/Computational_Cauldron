subroutine Grid_init()
           
             use Grid_data
             use MPI_data

#include "Solver.h"

            implicit none

            integer :: I

#ifdef SOLVER_GRID_UG

             Lx = 1
             Ly = 1

             Lx = Lx/HK
             Ly = Ly/HK

             dx=Lx/Nxb
             dy=Ly/Nyb

             inRe = .01

             do i=1,Nyb+1
                 x(:,i)=dx*(/(I,I=0,Nxb)/)
             enddo

             do i=1,Nxb+1
                 y(i,:)=dy*(/(I,I=0,Nyb)/)
             enddo

             y=y+(myid/HK)*Ly
             x=x+mod(myid,HK)*Lx                 

             t = 30.0

             dt = .05*min(dx,dy)

             !nt = 1
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
