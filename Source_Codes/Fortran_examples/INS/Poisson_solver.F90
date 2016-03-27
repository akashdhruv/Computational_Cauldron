subroutine Poisson_solver(u,v,p,dx,dy,dt,p_res,p_counter)

#include "Solver.h"
                
  implicit none

  real, dimension(Nxb+1,Nyb+2), intent(in) :: u
  real, dimension(Nxb+2,Nyb+1), intent(in) :: v

  real, dimension(Nxb+2,Nyb+2) :: p_old,p1,p2
  real, dimension(Nxb+2,Nyb+2), intent(inout) :: p
  real, intent(out) :: p_res
        
  integer, intent(out) :: p_counter
  integer :: i

  real, intent(in) :: dx,dy,dt

  p_counter = 0
 
  do while(p_counter<1500)
       
     p_old = p

     p1(2:Nxb+1,2:Nyb+1) = ((p_old(2:Nxb+1,3:Nyb+2)+p_old(2:Nxb+1,1:Nyb))/(dy*dy)) + ((p_old(3:Nxb+2,2:Nyb+1)+p_old(1:Nxb,2:Nyb+1))/(dx*dx))
   p2(2:Nxb+1,2:Nyb+1) = p1(2:Nxb+1,2:Nyb+1)-(1/(dy*dt))*(v(2:Nxb+1,2:Nyb+1)-v(2:Nxb+1,1:Nyb))-(1/(dx*dt))*(u(2:Nxb+1,2:Nyb+1)-u(1:Nxb,2:Nyb+1))
     p(2:Nxb+1,2:Nyb+1) = omega*(1/((2/(dx*dx))+(2/(dy*dy))))*p(2:Nxb+1,2:Nyb+1) + (1-omega)*p(2:Nxb+1,2:Nyb+1)

     ! Pressure BC
     p(:,1)=p(:,2)
     p(:,Nyb+2)=p(:,Nyb+1)

     p(1,:)=p(2,:)
     p(Nxb+2,:)=p(Nxb+1,:)


     p_counter = p_counter + 1

     do i=1,Nyb+2
          p_res = p_res + sum((p(:,i)-p_old(:,i))**2)
     enddo
     
     p_res = sqrt(p_res/((Nxb+2)*(Nyb+2)))

   if( (p_res .lt. 0.000001 ) .and. (p_res .ne. 0) ) exit

  end do

end subroutine Poisson_solver
