subroutine IO_display()

         use Grid_data
         use IncompNS_data

#include "Solver.h"

         implicit none

         print *,"Size of Domain ",Lx," X ",Ly
         print *,"Relaxation Factor ",omega
         print *,"Time Steps ",nt
         print *,"Poisson Iterations ",MaxIt

end subroutine IO_display
