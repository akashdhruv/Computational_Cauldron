subroutine IncompNS_rk3()

       use IO_interface, ONLY: IO_display, IO_write
       use Poisson_interface, ONLY: Poisson_solver            
       use IncompNS_data
       use Grid_data
       use MPI_data
       use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals

#include "Solver.h"

       implicit none
       
       real, dimension(Nxb+2,Nyb+2) :: ut
       real, dimension(Nxb+2,Nyb+2) :: vt

       real, dimension(Nxb+1,Nyb+1) :: uu
       real, dimension(Nxb+1,Nyb+1) :: vv

       real, dimension(Nxb+2,Nyb+2) :: u_old
       real, dimension(Nxb+2,Nyb+2) :: v_old

       real, dimension(Nxb,Nyb)   :: C1
       real, dimension(Nxb,Nyb)   :: G1
       real, dimension(Nxb,Nyb)   :: D1
       real, dimension(Nxb,Nyb)   :: G1_old

       real, dimension(Nxb,Nyb)   :: C2
       real, dimension(Nxb,Nyb)   :: G2
       real, dimension(Nxb,Nyb)   :: D2
       real, dimension(Nxb,Nyb)   :: G2_old

       real :: p_res, v_res, u_res, u_res1, v_res1

       integer :: tstep, p_counter, i

       ut=0
       vt=0

       u_old=0
       v_old=0
  
       C1 = 0
       G1 = 0
       D1 = 0
       G1_old = 0

       C2 = 0
       G2 = 0
       D2 = 0
       G2_old = 0

       tstep = 0


     do while (tstep<nt)

       v_res = 0
       u_res = 0

       v_res1 = 0
       u_res1 = 0

       u_old = u
       v_old = v

       ! Predictor Step

       call Convective_U(u,v,dx,dy,C1)
       call Diffusive_U(u,dx,dy,inRe,D1)
       G1 = C1 + D1

       if (tstep == 0) then

              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dt/1)*(G1)
              G1_old = G1
       else

              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dt/2)*(3*G1_old-G1)
              G1_old = G1
       endif


       call Convective_V(u,v,dx,dy,C2)
       call Diffusive_V(v,dx,dy,inRe,D2)
       G2 = C2 + D2

       if (tstep == 0) then

              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dt/1)*(G2)
              G2_old = G2
       else

              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dt/2)*(3*G2_old-G2)
              G2_old = G2
       endif

       ! Boundary Conditions

       if ( mod(myid,HK) == 0) then
 
           vt(1,:)=-vt(2,:)
           ut(1,:)=0
      
       end if

       if ( mod(myid,HK) == HK-1) then

           vt(Nxb+2,:)=-vt(Nxb+1,:)
           ut(Nxb+1,:)=0
           ut(Nxb+2,:)=0

       end if


       if ( myid/HK == 0) then

           vt(:,1)=0
           ut(:,1)=-ut(:,2)

       end if

       if ( myid/HK == HK-1) then
    
           vt(:,Nyb+2)=0
           vt(:,Nyb+1)=0
           ut(:,Nyb+2)=2-ut(:,Nyb+1)

       end if

       call MPI_applyBC(ut)
       call MPI_applyBC(vt)

       ! Poisson Solver

       call Poisson_solver(ut,vt,p_res,p_counter)

       u(2:Nxb+1,2:Nyb+1) = ut(2:Nxb+1,2:Nyb+1) - (dt/dx)*(p(3:Nxb+2,2:Nyb+1)-p(2:Nxb+1,2:Nyb+1))
       v(2:Nxb+1,2:Nyb+1) = vt(2:Nxb+1,2:Nyb+1) - (dt/dy)*(p(2:Nxb+1,3:Nyb+2)-p(2:Nxb+1,2:Nyb+1))

       ! Boundary Conditions

       if ( mod(myid,HK) == 0) then

           v(1,:)=-v(2,:)
           u(1,:)=0

       end if

       if ( mod(myid,HK) == HK-1) then

           v(Nxb+2,:)=-v(Nxb+1,:)
           u(Nxb+1,:)=0
           u(Nxb+2,:)=0

       end if


       if ( myid/HK == 0) then

           v(:,1)=0
           u(:,1)=-u(:,2)

       end if

       if ( myid/HK == HK-1) then

           v(:,Nyb+2)=0
           v(:,Nyb+1)=0
           u(:,Nyb+2)=2-u(:,Nyb+1)

       end if
      
       call MPI_applyBC(u)
       call MPI_applyBC(v)


       do i=1,Nyb+2
          u_res = u_res + sum((u(:,i)-u_old(:,i))**2)
       enddo

       call MPI_CollectResiduals(u_res,u_res1)
  
       u_res = sqrt(u_res1/((HK**HD)*(Nxb+2)*(Nyb+2)))

       do i=1,Nyb+1
          v_res = v_res + sum((v(:,i)-v_old(:,i))**2)
       enddo

       call MPI_CollectResiduals(v_res,v_res1)

       v_res = sqrt(v_res1/((HK**HD)*(Nxb+2)*(Nyb+2)))

       if (mod(tstep,5) == 0 .and. myid == 0) then       

          call IO_display(u_res,v_res,p_res,p_counter,tstep*dt)

       end if

       if( (u_res .lt. 0.000001) .and. (u_res .ne. 0).and. (v_res .lt. 0.000001) .and. (v_res .ne. 0) ) exit

       tstep = tstep +1

     end do

     uu = (u(1:Nxb+1,1:Nyb+1)+u(1:Nxb+1,2:Nyb+2))/2
     vv = (v(1:Nxb+1,1:Nyb+1)+v(2:Nxb+2,1:Nyb+1))/2

     call IO_write(x,y,uu,vv,myid)

end subroutine IncompNS_rk3


!! CONVECTIVE U !!
subroutine Convective_U(ut,vt,dx,dy,C1)

#include "Solver.h"
       
      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, dimension(Nxb,Nyb) :: ue
      real, dimension(Nxb,Nyb) :: uw
      real, dimension(Nxb,Nyb) :: us
      real, dimension(Nxb,Nyb) :: un
      real, dimension(Nxb,Nyb) :: vs
      real, dimension(Nxb,Nyb) :: vn
      real, dimension(Nxb,Nyb), intent(out) :: C1

      ue = (ut(2:Nxb+1,2:Nyb+1)+ut(3:Nxb+2,2:Nyb+1))/2
      uw = (ut(2:Nxb+1,2:Nyb+1)+ut(1:Nxb,2:Nyb+1))/2
      us = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,1:Nyb))/2
      un = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,3:Nyb+2))/2
      vs = (vt(2:Nxb+1,1:Nyb)+vt(3:Nxb+2,1:Nyb))/2
      vn = (vt(2:Nxb+1,2:Nyb+1)+vt(3:Nxb+2,2:Nyb+1))/2

      C1 = -((ue**2)-(uw**2))/dx - ((un*vn)-(us*vs))/dy

end subroutine Convective_U

!! CONVECTIVE V !!
subroutine Convective_V(ut,vt,dx,dy,C2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, dimension(Nxb,Nyb) :: vn, vs, ve, vw, ue, uw
      real, dimension(Nxb,Nyb), intent(out) :: C2

      vs = (vt(2:Nxb+1,2:Nyb+1)+vt(2:Nxb+1,1:Nyb))/2
      vn = (vt(2:Nxb+1,2:Nyb+1)+vt(2:Nxb+1,3:Nyb+2))/2
      ve = (vt(2:Nxb+1,2:Nyb+1)+vt(3:Nxb+2,2:Nyb+1))/2
      vw = (vt(2:Nxb+1,2:Nyb+1)+vt(1:Nxb,2:Nyb+1))/2
      ue = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,3:Nyb+2))/2
      uw = (ut(1:Nxb,2:Nyb+1)+ut(1:Nxb,3:Nyb+2))/2

      C2 = -((ue*ve)-(uw*vw))/dx - ((vn**2)-(vs**2))/dy

end subroutine Convective_V

!! DIFFUSIVE U !!
subroutine Diffusive_U(ut,dx,dy,inRe,D1)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      real, dimension(Nxb,Nyb) :: uP
      real, dimension(Nxb,Nyb) :: uN
      real, dimension(Nxb,Nyb) :: uS
      real, dimension(Nxb,Nyb) :: uE
      real, dimension(Nxb,Nyb) :: uW

      real, dimension(Nxb,Nyb), intent(out) :: D1

      uP = ut(2:Nxb+1,2:Nyb+1)
      uE = ut(3:Nxb+2,2:Nyb+1)
      uW = ut(1:Nxb,2:Nyb+1)
      uN = ut(2:Nxb+1,3:Nyb+2)
      uS = ut(2:Nxb+1,1:Nyb)

      D1 = (inRe/dx)*(((uE-uP)/dx)-((uP-uW)/dx)) + (inRe/dy)*(((uN-uP)/dy)-((uP-uS)/dy))

end subroutine Diffusive_U

!! DIFFUSIVE V !!
subroutine Diffusive_V(vt,dx,dy,inRe,D2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      real, dimension(Nxb,Nyb) :: vP,vE,vW,vN,vS

      real, dimension(Nxb,Nyb), intent(out) :: D2

      vP = vt(2:Nxb+1,2:Nyb+1)
      vE = vt(3:Nxb+2,2:Nyb+1)
      vW = vt(1:Nxb,2:Nyb+1)
      vN = vt(2:Nxb+1,3:Nyb+2)
      vS = vt(2:Nxb+1,1:Nyb)

      D2 = (inRe/dx)*(((vE-vP)/dx)-((vP-vW)/dx)) + (inRe/dy)*(((vN-vP)/dy)-((vP-vS)/dy))

end subroutine Diffusive_V

