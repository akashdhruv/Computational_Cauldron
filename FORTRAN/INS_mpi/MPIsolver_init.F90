subroutine MPIsolver_init()

    use MPI_data
    implicit none

    include "mpif.h"

    solver_comm = MPI_COMM_WORLD
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(solver_comm, myid, ierr)
    call MPI_COMM_SIZE(solver_comm, procs, ierr)

end subroutine MPIsolver_init 
