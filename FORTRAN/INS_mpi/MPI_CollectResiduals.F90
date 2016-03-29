subroutine MPI_CollectResiduals(res)

      use MPI_data

      implicit none

      include "mpif.h"
      
      real, intent(inout) :: res
      integer :: status(MPI_STATUS_SIZE)
      
      call MPI_ALLREDUCE(res,res,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

end subroutine MPI_CollectResiduals
