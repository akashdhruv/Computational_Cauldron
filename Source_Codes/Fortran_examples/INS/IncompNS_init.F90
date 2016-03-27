subroutine IncompNS_init()

      use Grid_data
      use IncompNS_data
      use IncompNS_interface, ONLY: IncompNS_rk3

#include "Solver.h"
   
      implicit none

      u = 0
      v = 0
      p = 0

      call IncompNS_rk3()
      
end subroutine IncompNS_init
