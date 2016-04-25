program Solver

use Grid_interface, ONLY: Grid_init
use IO_interface, ONLY: IO_display
use IncompNS_interface, ONLY: IncompNS_init

   implicit none

    real :: start,finish,exec_time

    call cpu_time(start)
    call Grid_init()
    call IncompNS_init()
    call cpu_time(finish)

    exec_time = finish-start

    print *,"Execution time = ",exec_time," seconds"

end program Solver


