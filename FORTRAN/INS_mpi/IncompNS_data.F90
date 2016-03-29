module IncompNS_data

#include "Solver.h"

     
        implicit none

        real, save, dimension(Nxb+2,Nyb+2) :: p 

        real, save, dimension(Nxb+2,Nyb+2) :: u

        real, save, dimension(Nxb+2,Nyb+2) :: v

end module IncompNS_data
