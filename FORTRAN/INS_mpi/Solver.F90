program Solver

use Grid_interface, ONLY: Grid_init
use IO_interface, ONLY: IO_display
use IncompNS_interface, ONLY: IncompNS_init
use MPI_interface, ONLY: MPIsolver_init, MPIsolver_finalize

   implicit none

    call MPIsolver_init()
    call Grid_init()
    call IncompNS_init()
    call MPIsolver_finalize()

end program Solver


