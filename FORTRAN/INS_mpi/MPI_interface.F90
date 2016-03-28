module MPI_interface

    interface
          subroutine MPIsolver_init()
          end subroutine MPIsolver_init
    end interface

    interface 
          subroutine MPIsolver_finalize()
          end subroutine MPIsolver_finalize
    end interface

!    interface
!          subroutine MPI_applyBCu()
!          end subroutine
!    end

!    interface
!          subroutine MPI_applyBCv()
!          end subroutine
!    end

!    interface 
!          subroutine MPI_applyBCp()
!          end subroutine
!    end

!    interface 
!       subroutine MPI_CollectResiduals()
!       end subroutine
!    end interface

end module MPI_interface
