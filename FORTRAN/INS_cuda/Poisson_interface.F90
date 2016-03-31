module Poisson_interface

#include "Solver.h"
 
       implicit none

       interface
             subroutine Poisson_solver(ut,vt,p_res,p_counter)
                implicit none
                real, dimension(Nxb+1,Nyb+2), intent(inout) :: ut
                real, dimension(Nxb+2,Nyb+1), intent(inout) :: vt
                real, intent(out) :: p_res
                integer, intent(out) :: p_counter
             end subroutine Poisson_solver
       end interface


!       interface
!             attributes(global) subroutine Poisson_kernel(p,u,v,Nx,Ny,dx,dy,dt)
!              implicit none
!              integer, value :: Nx,Ny
!              real, value :: dx, dy, dt
!              real, device :: p(Nx+2,Ny+2), u(Nx+2,Ny+2), v(Nx+2,Ny+2)         

!             end subroutine Poisson_kernel
!       end interface

end module Poisson_interface
