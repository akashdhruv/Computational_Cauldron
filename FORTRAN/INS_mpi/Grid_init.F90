subroutine Grid_init()
           
     use Grid_data
     use MPI_data

#include "Solver.h"

     implicit none

     integer :: I
     real :: pi=4.0*atan(1.0)

     Lx = (D_xmax)-(D_xmin)
     Ly = (D_ymax)-(D_ymin)

     Lx = Lx/HK
     Ly = Ly/HK

     inRe = .001

     dx = Lx/Nxb
     dy = Ly/Nyb


     t = 200.0

#ifdef SOLVER_GRID_UG

     do i=1,Nyb+1
        x(:,i)=D_xmin+mod(myid,HK)*Lx+dx*(/(I,I=0,Nxb)/)
     enddo

     do i=1,Nxb+1
        y(i,:)=D_ymin+(myid/HK)*Ly+dy*(/(I,I=0,Nyb)/)
    enddo
   

#endif

#ifdef SOLVER_GRID_NON_UG

     do i=1,Nyb+1

        if (mod(myid,HK) == 0) then
           x(:,i)=D_xmin+mod(myid,HK)*Lx+Lx*(1-cos((pi/(2*Nxb))*(/(I,I=0,Nxb)/)))
       
        else if (mod(myid,HK) == HK-1) then
           x(:,i)=D_xmin+mod(myid,HK)*Lx+Lx*(cos((pi/(2*Nxb))*(/(I,I=Nxb,0,-1)/)))
       
        else
           !x(:,i)=D_xmin+mod(myid,HK)*Lx+dx*(/(I,I=0,Nxb)/)
          
           x(1:Nxb/2+1,i)=D_xmin+mod(myid,HK)*Lx+Lx*(cos((pi/(Nxb))*(/(I,I=Nxb/2,0,-1)/)))/2
          
           x(Nxb/2+2:Nxb+1,i)=D_xmin+mod(myid,HK)*Lx+Lx/2+&
                              Lx*(1-cos((pi/(Nxb))*(/(I,I=1,Nxb/2)/)))/2

        end if

      enddo

     do i=1,Nxb+1

        if (myid/HK == 0) then
            y(i,:)=D_ymin+(myid/HK)*Ly+Ly*(1-cos((pi/(2*Nyb))*(/(I,I=0,Nyb)/)))
   
        else if (myid/HK == HK-1) then
            y(i,:)=D_ymin+(myid/HK)*Ly+Ly*(cos((pi/(2*Nyb))*(/(I,I=Nyb,0,-1)/)))

        else
          !y(i,:)=D_ymin+(myid/HK)*Ly+dy*(/(I,I=0,Nyb)/)
          
          y(i,1:Nyb/2+1)=D_ymin+(myid/HK)*Ly+Ly*(cos((pi/(Nyb))*(/(I,I=Nyb/2,0,-1)/)))/2
          
          y(i,Nyb/2+2:Nyb+1)=D_ymin+(myid/HK)*Ly+Ly/2+&
                               Ly*(1-cos((pi/(Nyb))*(/(I,I=1,Nyb/2)/)))/2
        end if

     enddo

#endif


dx_nodes(2:Nxb+1,1:Nyb+1) = x(2:Nxb+1,:)-x(1:Nxb,:)
dy_nodes(1:Nxb+1,2:Nyb+1) = y(:,2:Nyb+1)-y(:,1:Nyb)

dx_nodes(1,:) = dx_nodes(2,:)
dy_nodes(:,1) = dy_nodes(:,2)

dx_nodes(:,1) = dx_nodes(:,2)
dy_nodes(1,:) = dy_nodes(2,:)

dx_nodes(Nxb+2,:) = dx_nodes(Nxb+1,:)
dy_nodes(:,Nyb+2) = dy_nodes(:,Nyb+1)

dx_nodes(:,Nyb+2) = dx_nodes(:,Nyb+1)
dy_nodes(Nxb+2,:) = dy_nodes(Nxb+1,:)

dx_centers = 0.5*(dx_nodes(1:Nxb+1,1:Nyb+1)+dx_nodes(2:Nxb+2,1:Nyb+1))
dy_centers = 0.5*(dy_nodes(1:Nxb+1,1:Nyb+1)+dy_nodes(1:Nxb+1,2:Nyb+2)) 

dt = 0.03*min(minval(dx_nodes),minval(dy_nodes))

nt = t/dt
!nt = 1

end subroutine Grid_init
