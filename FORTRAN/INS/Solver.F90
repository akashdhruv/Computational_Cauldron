program Solver

use omp_lib
use Grid_interface, ONLY: Grid_init
use IO_interface, ONLY: IO_display
use IncompNS_interface, ONLY: IncompNS_init

   implicit none

    double precision :: start,finish,exec_time

    start = omp_get_wtime()

    call Grid_init()
    call IncompNS_init()

    finish = omp_get_wtime()

    exec_time = finish-start

    !print *,"Execution time = ",exec_time," seconds"
    print '("Execution time: ",f10.7," seconds")',exec_time

end program Solver


