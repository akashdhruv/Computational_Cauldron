subroutine MPIsolver_init()

    use MPI_data
    implicit none

    include "mpif.h"

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, procs, ierr)

end subroutine MPIsolver_init 
