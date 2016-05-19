subroutine IncompNS_init()

      use Grid_data

#include "Solver.h"
   
      implicit none
      
      real,pointer,dimension(:,:) :: u,v,p

      p => gr_center(PRES_VAR,:,:)
      u => gr_facex(VELC_VAR,:,:)
      v => gr_facey(VELC_VAR,:,:)


      p = 0.0
      u = 0.0
      v = 0.0

      nullify(p)
      nullify(u)
      nullify(v)
      
end subroutine IncompNS_init
