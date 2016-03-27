subroutine IncompNS_rk3()

       use IO_interface, ONLY: IO_display
       use Poisson_interface, ONLY: Poisson_solver            
       use IncompNS_data
       use Grid_data

#include "Solver.h"

       implicit none
       
       real, dimension(Nxb+1,Nyb+2) :: ut
       real, dimension(Nxb+2,Nyb+1) :: vt

       real, dimension(Nxb+1,Nyb+2) :: u_old
       real, dimension(Nxb+2,Nyb+1) :: v_old

       real, dimension(Nxb-1,Nyb)   :: C1
       real, dimension(Nxb-1,Nyb)   :: G1
       real, dimension(Nxb-1,Nyb)   :: D1
       real, dimension(Nxb-1,Nyb)   :: G1_old

       real, dimension(Nxb,Nyb-1)   :: C2
       real, dimension(Nxb,Nyb-1)   :: G2
       real, dimension(Nxb,Nyb-1)   :: D2
       real, dimension(Nxb,Nyb-1)   :: G2_old

       real :: p_res, v_res, u_res

       integer :: tstep, p_counter, i

       tstep = 0

     do while (tstep<nt)

       u_old = u
       v_old = v

       ! Predictor Step

       call Convective_U(ut,vt,dx,dy,C1)
       call Diffusive_U(ut,vt,dx,dy,inRe,D1)
       G1 = C1 + D1

       if (tstep == 0) then

              ut(2:Nxb,2:Nyb+1)=ut(2:Nxb,2:Nyb+1)+(dt/1)*(G1)
              G1_old = G1
       else

              ut(2:Nxb,2:Nyb+1)=ut(2:Nxb,2:Nyb+1)+(dt/2)*(3*G1_old-G1)
              G1_old = G1
       endif


       call Convective_V(ut,vt,dx,dy,C2)
       call Diffusive_V(ut,vt,dx,dy,inRe,D2)
       G2 = C2 + D2

       if (tstep == 0) then

              vt(2:Nxb+1,2:Nyb)=vt(2:Nxb+1,2:Nyb)+(dt/1)*(G2)
              G2_old = G2
       else

              vt(2:Nxb+1,2:Nyb)=ut(2:Nxb+1,2:Nyb)+(dt/2)*(3*G2_old-G2)
              G2_old = G2
       endif

       ! Boundary Conditions
       vt(1,:)=-vt(2,:)
       vt(Nxb+2,:)=-vt(Nxb+1,:)

       ut(1,:)=0
       ut(Nxb+1,:)=0

       vt(:,1)=0
       vt(:,Nyb+1)=0

       ut(:,1)=-ut(:,2)
       ut(:,Nyb+2)=2-ut(:,Nyb+1)

       ! Poisson Solver

       call Poisson_solver(ut,vt,p,dx,dy,dt,p_res,p_counter)

       u(2:Nxb,2:Nyb+1) = u(2:Nxb,2:Nyb+1) + dt*(p(3:Nxb+1,2:Nyb+1)-p(2:Nxb,2:Nyb+1))
       v(2:Nxb+1,2:Nyb) = v(2:Nxb+1,2:Nyb) + dt*(p(2:Nxb+1,3:Nyb+1)-p(2:Nxb+1,2:Nyb))


       ! Boundary Conditions

       v(1,:)=-v(2,:)
       v(Nxb+2,:)=-v(Nxb+1,:)

       u(1,:)=0
       u(Nxb+1,:)=0

       v(:,1)=0
       v(:,Nyb+1)=0

       u(:,1)=-u(:,2)
       u(:,Nyb+2)=2-u(:,Nyb+1)

       do i=1,Nyb+2
          u_res = u_res + sum((u(:,i)-u_old(:,i))**2)
       enddo
       u_res = sqrt(u_res/((Nxb+1)*(Nyb+2)))

       do i=1,Nyb+1
          v_res = v_res + sum((v(:,i)-v_old(:,i))**2)
       enddo
       v_res = sqrt(v_res/((Nxb+2)*(Nyb+1)))

       if (mod(tstep,5) == 0) then       

          call IO_display(u_res,v_res,p_res,p_counter)
      
          if( (u_res .lt. 0.0000001) .and. (u_res .ne. 0).and. (v_res .lt. 0.0000001) .and. (v_res .ne. 0) ) exit

       endif

       tstep = tstep +1

     end do

end subroutine IncompNS_rk3


!! CONVECTIVE U !!
subroutine Convective_U(ut,vt,dx,dy,C1)

#include "Solver.h"
       
      implicit none

      real,dimension(Nxb+1,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+1), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, dimension(Nxb-1,Nyb) :: ue
      real, dimension(Nxb-1,Nyb) :: uw
      real, dimension(Nxb-1,Nyb) :: us
      real, dimension(Nxb-1,Nyb) :: un
      real, dimension(Nxb-1,Nyb) :: vs
      real, dimension(Nxb-1,Nyb) :: vn
      real, dimension(Nxb-1,Nyb), intent(out) :: C1

      ue = (ut(2:Nxb,2:Nyb+1)+ut(3:Nxb+1,2:Nyb+1))/2
      uw = (ut(2:Nxb,2:Nyb+1)+ut(1:Nxb-1,2:Nyb+1))/2
      us = (ut(2:Nxb,2:Nyb+1)+ut(2:Nxb,1:Nyb))/2
      un = (ut(2:Nxb,2:Nyb+1)+ut(2:Nxb,3:Nyb+2))/2
      vs = (vt(2:Nxb,1:Nyb)+vt(3:Nxb+1,1:Nyb))/2
      vn = (vt(2:Nxb,2:Nyb+1)+vt(3:Nxb+1,2:Nyb+1))/2

      C1 = -((ue**2)-(uw**2))/dx - ((un*vn)-(us*vs))/dy

end subroutine Convective_U

!! CONVECTIVE V !!
subroutine Convective_V(ut,vt,dx,dy,C2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+1,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+1), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, dimension(Nxb,Nyb-1) :: vn, vs, ve, vw, ue, uw
      real, dimension(Nxb,Nyb-1), intent(out) :: C2

      vn = (vt(2:Nxb+1,2:Nyb)+vt(2:Nxb+1,1:Nyb-1))/2
      vs = (vt(2:Nxb+1,2:Nyb)+vt(2:Nxb+1,3:Nyb+1))/2
      ve = (vt(2:Nxb+1,2:Nyb)+vt(3:Nxb+2,2:Nyb))/2
      vw = (vt(2:Nxb+1,2:Nyb)+vt(1:Nxb,2:Nyb))/2
      ue = (ut(2:Nxb+1,2:Nyb)+ut(2:Nxb+1,3:Nyb+1))/2
      uw = (ut(1:Nxb,2:Nyb)+ut(1:Nxb,3:Nyb+1))/2

      C2 = -((ue*ve)-(uw*vw))/dx - ((vn**2)-(vs**2))/dy

end subroutine Convective_V

!! DIFFUSIVE U !!
subroutine Diffusive_U(ut,vt,dx,dy,inRe,D1)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+1,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+1), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      real, dimension(Nxb-1,Nyb) :: uP
      real, dimension(Nxb-1,Nyb) :: uN
      real, dimension(Nxb-1,Nyb) :: uS
      real, dimension(Nxb-1,Nyb) :: uE
      real, dimension(Nxb-1,Nyb) :: uW

      real, dimension(Nxb-1,Nyb), intent(out) :: D1

      uP = ut(2:Nxb,2:Nyb+1)
      uE = ut(3:Nxb+1,2:Nyb+1)
      uW = ut(1:Nxb-1,2:Nyb+1)
      uN = ut(2:Nxb,3:Nyb+2)
      uS = ut(2:Nxb,1:Nyb)

      D1 = (inRe/dx)*(((uE-uP)/dx)-((uP-uW)/dx)) + (inRe/dy)*(((uN-uP)/dy)-((uP-uS)/dy))

end subroutine Diffusive_U

!! DIFFUSIVE V !!
subroutine Diffusive_V(ut,vt,dx,dy,inRe,D2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+1,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+1), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      real, dimension(Nxb,Nyb-1) :: vP,vE,vW,vN,vS

      real, dimension(Nxb,Nyb-1), intent(out) :: D2

      vP = vt(2:Nxb+1,2:Nyb)
      vE = vt(3:Nxb+2,2:Nyb)
      vW = vt(1:Nxb,2:Nyb)
      vN = vt(2:Nxb+1,3:Nyb+1)
      vS = vt(2:Nxb+1,1:Nyb-1)

      D2 = (inRe/dx)*(((vE-vP)/dx)-((vP-vW)/dx)) + (inRe/dy)*(((vN-vP)/dy)-((vP-vS)/dy))

end subroutine Diffusive_V

