subroutine IO_display(u_res,v_res,p_res,p_counter)

#include "Solver.h"

         implicit none

         real, intent(in) :: u_res,v_res,p_res
         integer, intent(in) :: p_counter

         print *,"*****************PARAMETER DISPLAY**********************"
         print *,"U Velocity Residual: ",u_res
         print *,"V Velocity Residual: ",v_res
         print *,"Pressure Residual:   ",p_res
         print *,"Poisson Counter:     ",p_counter

end subroutine IO_display
