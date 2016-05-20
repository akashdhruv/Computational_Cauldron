subroutine HeatAD_init()

#include "Solver.h"

   use HeatAD_data
   use IncompNS_data
   use physicaldata

   implicit none

   real,pointer,dimension(:,:) :: T

   ht_Pr = 0.7
   ht_Nu = 0.332*(ht_Pr**0.33)/(ins_inRe**0.5)

   T => ph_center(TEMP_VAR,:,:)
  
   T = 313.0

   nullify(T)


end subroutine HeatAD_init
