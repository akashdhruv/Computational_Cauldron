module IncompNS_data

        use Grid_data

#include "Solver.h"

     
        implicit none

        real, save, dimension(Nxb+2,Nyb+2) :: p 

        real, save, dimension(Nxb+1,Nyb+2) :: u

        real, save, dimension(Nxb+2,Nyb+1) :: v

end module IncompNS_data
