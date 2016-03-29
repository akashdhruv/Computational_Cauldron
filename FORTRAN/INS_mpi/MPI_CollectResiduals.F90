subroutine MPI_CollectResiduals(res,res1)

      use MPI_data

      implicit none

      include "mpif.h"
      
      real, intent(inout) :: res,res1
      integer :: status(MPI_STATUS_SIZE)
      
      call MPI_ALLREDUCE(res,res1,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

end subroutine MPI_CollectResiduals
