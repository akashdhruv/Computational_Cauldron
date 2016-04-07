subroutine MPI_applyBC(u_ex)

#include "Solver.h"

       use MPI_data

       implicit none

       include "mpif.h"

       real, dimension(Nxb+2,Nyb+2), intent(inout) :: u_ex
       integer :: status(MPI_STATUS_SIZE)

       call MPI_BARRIER(x_comm, ierr)

       if(mod(x_id,2) == 0) then
           
             if(x_id == 0) then

                  call MPI_SEND(u_ex(Nxb+1,:), Nyb+2, MPI_REAL, x_id + 1, 1,x_comm, ierr)
                  call MPI_RECV(u_ex(Nxb+2,:), Nyb+2, MPI_REAL, x_id + 1, 2,x_comm, status, ierr) 

             else if(x_id == HK-1) then
            
                  call MPI_SEND(u_ex(2,:), Nyb+2, MPI_REAL, x_id-1, 3,x_comm, ierr)
                  call MPI_RECV(u_ex(1,:), Nyb+2, MPI_REAL, x_id-1, 4,x_comm,status, ierr)

             else
                  call MPI_SEND(u_ex(Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 1,x_comm, ierr)
                  call MPI_RECV(u_ex(Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 2,x_comm,status, ierr) 
                  
                  call MPI_SEND(u_ex(2,:), Nyb+2, MPI_REAL, x_id-1, 3,x_comm, ierr)
                  call MPI_RECV(u_ex(1,:), Nyb+2, MPI_REAL, x_id-1, 4,x_comm,status, ierr)
                                   

             end if

       else if (mod(x_id,2) == 1) then

             if(x_id == HK-1) then
            
                  call MPI_RECV(u_ex(1,:), Nyb+2, MPI_REAL,x_id-1, 1,x_comm,status, ierr)
                  call MPI_SEND(u_ex(2,:), Nyb+2, MPI_REAL,x_id-1, 2,x_comm, ierr)

             else
                  call MPI_RECV(u_ex(1,:), Nyb+2, MPI_REAL,x_id-1, 1,x_comm,status, ierr)
                  call MPI_SEND(u_ex(2,:), Nyb+2, MPI_REAL,x_id-1, 2,x_comm,ierr) 

                  call MPI_RECV(u_ex(Nxb+2,:), Nyb+2, MPI_REAL,x_id+1, 3,x_comm,status,ierr)
                  call MPI_SEND(u_ex(Nyb+1,:), Nyb+2, MPI_REAL,x_id+1, 4,x_comm,ierr)
                  

             end if

       end if

       call MPI_BARRIER(x_comm,ierr)
  
       call MPI_BARRIER(y_comm,ierr)
 
       !! Second dimension !!

       if(mod(y_id,2) == 0) then

             if(y_id == 0) then

                  call MPI_SEND(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 5,y_comm,ierr)
                  call MPI_RECV(u_ex(:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 6,y_comm,status, ierr)

             else if(y_id == HK-1) then

                  call MPI_SEND(u_ex(:,2), Nxb+2, MPI_REAL, y_id-1, 7,y_comm, ierr)
                  call MPI_RECV(u_ex(:,1), Nxb+2, MPI_REAL, y_id-1, 8,y_comm,status,ierr)

             else 
                  call MPI_SEND(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 5,y_comm, ierr)
                  call MPI_RECV(u_ex(:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 6,y_comm,status,ierr)

                  call MPI_SEND(u_ex(:,2), Nxb+2, MPI_REAL, y_id-1, 7,y_comm, ierr)
                  call MPI_RECV(u_ex(:,1), Nxb+2, MPI_REAL, y_id-1, 8,y_comm,status,ierr)
                  

             end if

       else if (mod(y_id,2) == 1) then

             if(y_id == HK-1) then

                  call MPI_RECV(u_ex(:,1), Nxb+2, MPI_REAL, y_id-1, 5,y_comm,status, ierr)
                  call MPI_SEND(u_ex(:,2), Nxb+2, MPI_REAL, y_id-1, 6,y_comm, ierr)

             else
                  call MPI_RECV(u_ex(:,1), Nxb+2, MPI_REAL, y_id-1, 5,y_comm,status, ierr)
                  call MPI_SEND(u_ex(:,2), Nxb+2, MPI_REAL, y_id-1, 6,y_comm,ierr)

                  call MPI_RECV(u_ex(:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 7,y_comm,status, ierr)
                  call MPI_SEND(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 8,y_comm,ierr)
                  

             end if

       end if            

       call MPI_BARRIER(y_comm,ierr)

end subroutine

