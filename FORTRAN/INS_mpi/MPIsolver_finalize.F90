subroutine MPIsolver_finalize()

      use MPI_data

      implicit none

      include "mpif.h"

      call MPI_FINALIZE(ierr)

end subroutine MPIsolver_finalize
