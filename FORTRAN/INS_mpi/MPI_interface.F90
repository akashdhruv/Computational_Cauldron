module MPI_interface
#include "Solver.h"

    implicit none

    interface
          subroutine MPIsolver_init()
          implicit none
          end subroutine MPIsolver_init
    end interface

    interface 
          subroutine MPIsolver_finalize()
          implicit none
          end subroutine MPIsolver_finalize
    end interface

    interface
          subroutine MPI_applyBC(u_ex)
          implicit none
          real, dimension(Nxb+2,Nyb+2),intent(inout) :: u_ex
          end subroutine MPI_applyBC
    end interface


    interface 
       subroutine MPI_CollectResiduals(res)
       implicit none
       real, intent(inout) :: res
       end subroutine MPI_CollectResiduals
    end interface

end module MPI_interface
