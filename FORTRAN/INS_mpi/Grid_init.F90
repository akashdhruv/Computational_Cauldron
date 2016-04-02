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


     t = 50.0

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
           x(:,i)=D_xmin+mod(myid,HK)*Lx+dx*(/(I,I=0,Nxb)/)
        end if

      enddo

     do i=1,Nxb+1

        if (myid/HK == 0) then
            y(i,:)=D_ymin+(myid/HK)*Ly+Ly*(1-cos((pi/(2*Nyb))*(/(I,I=0,Nyb)/)))
   
        else if (myid/HK == HK-1) then
            y(i,:)=D_ymin+(myid/HK)*Ly+Ly*(cos((pi/(2*Nyb))*(/(I,I=Nyb,0,-1)/)))

        else
           y(i,:)=D_ymin+(myid/HK)*Ly+dy*(/(I,I=0,Nyb)/)
        end if

     enddo

#endif


dx_n(1:Nxb,:) = x(2:Nxb+1,:)-x(1:Nxb,:)
dy_n(:,1:Nyb) = y(:,2:Nyb+1)-y(:,1:Nyb)

dx_n(Nxb+1,:) = dx_n(Nxb,:)
dy_n(:,Nyb+1) = dy_n(:,Nyb)

dx_a = 0.5*(dx_n(1:Nxb,1:Nyb)+dx_n(2:Nxb+1,1:Nyb))
dy_a = dy_n(1:Nxb,1:Nyb)

dx_b = dx_n(1:Nxb,1:Nyb)
dy_b = 0.5*(dy_n(1:Nxb,1:Nyb)+dy_n(1:Nxb,2:Nyb+1))

dx_c(2:Nxb+1,:) = 0.5*(dx_n(1:Nxb,:)+dx_n(2:Nxb+1,:))
dx_c(1,:) = dx_c(2,:)

dy_c(:,2:Nyb+1) = 0.5*(dy_n(:,1:Nyb)+dy_n(:,2:Nyb+1))
dy_c(:,1) = dy_c(:,2)

!dx_u(1:Nxb,1:Nyb) = dx
!dx_u(Nxb+1,1:Nyb) = dx(Nxb,:)
!dx_u(:,Nyb+1) = dx_u(:,Nyb)

!dy_u(1:Nxb,1:Nyb) = 0.5*(dy(1:Nxb,1:Nyb)+dy(1:Nxb,2:Nyb+1))
!dy_u(1:Nxb,Nyb+1) = dy(:,Nyb)
!dy_u(Nxb+1,:)=dy_u(Nxb,:)

!dx_v(1:Nxb,1:Nyb) = 0.5*(dx(1:Nxb,1:Nyb)+dx(2:Nxb+1,1:Nyb))
!dx_v(Nxb+1,1:Nyb) = dx(Nxb,:)
!dx_v(:,Nyb+1) = dx_v(:,Nyb)

!dy_v(1:Nxb,1:Nyb) = dy
!dy_v(1:Nxb,Nyb+1) = dy(:,Nyb)
!dy_v(Nxb+1,:) = dy_v(Nxb,:)


!dt = 0.05*min(minval(dx),minval(dy))
dt=0.05*min(dx,dy)
!dt = .0000001

nt = t/dt
!nt = 1
end subroutine Grid_init
