subroutine IncompNS_rk3()
            
       use IncompNS_data
       use Grid_data

#include "Solver.h"

       implicit none

       call Convective_u()
       call Diffusive()


end subroutine IncompNS_rk3


subroutine Convective_u()
       
      implicit none

      print *,"Convective U"


end subroutine Convective_u


subroutine Diffusive()

      implicit none

      print *,"Diffusive"

end subroutine Diffusive
