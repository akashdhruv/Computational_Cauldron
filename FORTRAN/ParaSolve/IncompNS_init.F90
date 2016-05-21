subroutine IncompNS_init()

      use physicaldata
      use IncompNS_data

#include "Solver.h"
   
      implicit none
      
      real,pointer,dimension(:,:) :: u,v,p

      p => ph_center(PRES_VAR,:,:)
      u => ph_facex(VELC_VAR,:,:)
      v => ph_facey(VELC_VAR,:,:)


      p = 0.0
      u = 0.0
      v = 0.0

      nullify(p)
      nullify(u)
      nullify(v)
      
      ins_inRe = 0.001

end subroutine IncompNS_init
