subroutine Solver_init

    use MPI_interface, only: MPIsolver_init
    use Grid_interface, only: Grid_init
    use IncompNS_interface, only: IncompNS_init

    implicit none

    call MPIsolver_init()
    call Grid_init()
    call IncompNS_init()
    call HeatAD_init()
    call Driver_init()

end subroutine Solver_init
