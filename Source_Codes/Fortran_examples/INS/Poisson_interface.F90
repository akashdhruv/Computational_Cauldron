module Poisson_interface

#include "Solver.h"
 
       implicit none

       interface
             subroutine Poisson_solver(u,v,p,dx,dy,dt,p_res,p_counter)
                implicit none
                real, dimension(Nxb+2,Nyb+2), intent(inout) :: p
                real, dimension(Nxb+1,Nyb+2), intent(inout) :: u
                real, dimension(Nxb+2,Nyb+1), intent(inout) :: v
                real, intent(out) :: p_res
                integer, intent(out) :: p_counter
                real, intent(in) :: dx,dy,dt
             end subroutine Poisson_solver

       end interface


end module Poisson_interface
