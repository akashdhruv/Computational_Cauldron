subroutine IncompNS_init()

      use Grid_data
      use IncompNS_data
      use IncompNS_interface, ONLY: IncompNS_rk3

#include "Solver.h"
   
      implicit none

      call IncompNS_rk3()
      
end subroutine IncompNS_init
