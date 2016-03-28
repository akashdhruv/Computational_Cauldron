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


end module Poisson_interface
