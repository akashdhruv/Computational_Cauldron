subroutine MPI_applyBC(u_ex)

#include "Solver.h"

       use MPI_data

       implicit none

       include "mpif.h"

       real, dimension(Nxb+2,Nyb+2), intent(inout) :: u_ex
       integer :: status(MPI_STATUS_SIZE)

       if(mod(mod(myid,HK),2) == 0) then
           
             if(mod(myid,HK) == 0) then

                  call MPI_SEND(u_ex(Nxb+1,:), Nyb+2, MPI_REAL, myid+(HK**0), myid+(HK**0),MPI_COMM_WORLD, ierr)
                  call MPI_RECV(u_ex(Nxb+2,:), Nyb+2, MPI_REAL, myid+(HK**0), myid, MPI_COMM_WORLD, status, ierr) 

             else if(mod(myid,HK) == HK-1) then
            
                  call MPI_SEND(u_ex(2,:), Nyb+2, MPI_REAL, myid-(HK**0), myid-(HK**0),MPI_COMM_WORLD, ierr)
                  call MPI_RECV(u_ex(1,:), Nyb+2, MPI_REAL, myid-(HK**0), myid,MPI_COMM_WORLD,status, ierr)

             else
                  call MPI_SEND(u_ex(2,:), Nyb+2, MPI_REAL, myid-(HK**0), myid-(HK**0),MPI_COMM_WORLD, ierr)
                  call MPI_RECV(u_ex(1,:), Nyb+2, MPI_REAL, myid-(HK**0), myid,MPI_COMM_WORLD,status, ierr)
       
                  call MPI_SEND(u_ex(Nxb+1,:), Nyb+2, MPI_REAL, myid+(HK**0), myid+(HK**0),MPI_COMM_WORLD, ierr)
                  call MPI_RECV(u_ex(Nxb+2,:), Nyb+2, MPI_REAL, myid+(HK**0), myid,MPI_COMM_WORLD,status, ierr)                  

             end if

       else if (mod(mod(myid,HK),2) == 1) then

             if(mod(myid,HK) == HK-1) then
            
                  call MPI_RECV(u_ex(1,:), Nyb+2, MPI_REAL, myid-(HK**0), myid,MPI_COMM_WORLD,status, ierr)
                  call MPI_SEND(u_ex(2,:), Nyb+2, MPI_REAL, myid-(HK**0), myid-(HK**0),MPI_COMM_WORLD, ierr)

             else
                  call MPI_RECV(u_ex(Nxb+2,:), Nyb+2, MPI_REAL, myid+(HK**0), myid,MPI_COMM_WORLD,status,ierr)
                  call MPI_SEND(u_ex(Nyb+1,:), Nyb+2, MPI_REAL, myid+(HK**0), myid+(HK**0), MPI_COMM_WORLD,ierr)
       
                  call MPI_RECV(u_ex(1,:), Nyb+2, MPI_REAL, myid-(HK**0), myid,MPI_COMM_WORLD,status, ierr)
                  call MPI_SEND(u_ex(2,:), Nyb+2, MPI_REAL, myid-(HK**0), myid+(HK**0),MPI_COMM_WORLD,ierr) 

             end if

       end if

!       call MPI_BARRIER(MPI_COMM_WORLD)

       !! Second dimension !!

       if(mod(myid/HK,2) == 0) then

             if(myid/HK == 0) then

                  call MPI_SEND(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, myid+(HK**1), myid+(HK**1), MPI_COMM_WORLD,ierr)
                  call MPI_RECV(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, myid+(HK**1), myid,MPI_COMM_WORLD,status, ierr)

             else if(myid/HK == HK-1) then

                  call MPI_SEND(u_ex(:,2), Nxb+2, MPI_REAL, myid-(HK**1), myid-(HK**1),MPI_COMM_WORLD, ierr)
                  call MPI_RECV(u_ex(:,1), Nxb+2, MPI_REAL, myid-(HK**1), myid, MPI_COMM_WORLD,status,ierr)

             else
                  call MPI_SEND(u_ex(:,2), Nxb+2, MPI_REAL, myid-(HK**1), myid-(HK**1),MPI_COMM_WORLD, ierr)
                  call MPI_RECV(u_ex(:,1), Nxb+2, MPI_REAL, myid-(HK**1), myid, MPI_COMM_WORLD,status,ierr)

                  call MPI_SEND(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, myid+(HK**1), myid-(HK**1),MPI_COMM_WORLD, ierr)
                  call MPI_RECV(u_ex(:,Nyb+2), Nxb+2, MPI_REAL, myid+(HK**1), myid, MPI_COMM_WORLD,status,ierr)

             end if

       else if (mod(myid/HK,2) == 1) then

             if(myid/HK == HK-1) then

                  call MPI_RECV(u_ex(:,1), Nxb+2, MPI_REAL, myid-(HK**1), myid,MPI_COMM_WORLD,status, ierr)
                  call MPI_SEND(u_ex(:,2), Nxb+2, MPI_REAL, myid-(HK**1), myid-(HK**1),MPI_COMM_WORLD, ierr)

             else
                  call MPI_RECV(u_ex(:,Nyb+2), Nxb+2, MPI_REAL, myid+(HK**1), myid,MPI_COMM_WORLD,status, ierr)
                  call MPI_SEND(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, myid+(HK**1), myid+(HK**1),MPI_COMM_WORLD,ierr)

                  call MPI_RECV(u_ex(:,1), Nxb+2, MPI_REAL, myid-(HK**1), myid,MPI_COMM_WORLD,status, ierr)
                  call MPI_SEND(u_ex(:,2), Nxb+2, MPI_REAL, myid-(HK**1), myid+(HK**1),MPI_COMM_WORLD,ierr)

             end if

       end if            

!       call MPI_BARRIER(MPI_COMM_WORLD)

end subroutine

