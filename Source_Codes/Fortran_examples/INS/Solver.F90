program Solver

use Grid_interface, ONLY: Grid_init
use IO_interface, ONLY: IO_display
use IncompNS_interface, ONLY: IncompNS_init

   implicit none

    call Grid_init()
    call IncompNS_init()

end program Solver


