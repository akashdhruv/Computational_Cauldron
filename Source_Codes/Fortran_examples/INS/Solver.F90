program Solver

use Grid_interface, ONLY: Grid_init
use IO_interface, ONLY: IO_display

   implicit none

    call Grid_init()
    call IO_display()

end program Solver


