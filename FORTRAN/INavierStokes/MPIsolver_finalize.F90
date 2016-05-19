subroutine MPIsolver_finalize()

      use MPI_data

      implicit none

      include "mpif.h"

      call MPI_COMM_FREE(x_comm,ierr)
      call MPI_COMM_FREE(y_comm,ierr)

   
      call cpu_time(finish)
      exec_time = (finish - start) 
      print *,"Execution Time = ",exec_time," seconds"

      call MPI_FINALIZE(ierr)

end subroutine MPIsolver_finalize
