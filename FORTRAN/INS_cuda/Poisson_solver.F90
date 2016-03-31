subroutine Poisson_solver(ut,vt,p_res,p_counter)

  use IncompNS_data
  use Grid_data 
  !use cudafor

#include "Solver.h"
                
  implicit none

  !integer, parameter :: NG_x = (Nxb+2)/2, NG_y = (Nyb+2/2

  !TYPE(dim3), parameter :: gridSize,blockSize
 
  real, dimension(Nxb+1,Nyb+2), intent(in) :: ut
  real, dimension(Nxb+2,Nyb+1), intent(in) :: vt

  real, dimension(Nxb+2,Nyb+2) :: p_old

  real, intent(out) :: p_res
        
  integer, intent(out) :: p_counter
  integer :: i,j

  p_old = 0
  p_counter = 0
 
  do while(p_counter<MaxIt)

     p_res = 0       
     p_old = p

#ifdef POISSON_SOLVER_JACOBI

     p(2:Nxb+1,2:Nyb+1)= (((p_old(2:Nxb+1,3:Nyb+2)+p_old(2:Nxb+1,1:Nyb))/(dy*dy))& 
                         +((p_old(3:Nxb+2,2:Nyb+1)+p_old(1:Nxb,2:Nyb+1))/(dx*dx))&
                         -((1/(dy*dt))*(vt(2:Nxb+1,2:Nyb+1)-vt(2:Nxb+1,1:Nyb)))&
                         -((1/(dx*dt))*(ut(2:Nxb+1,2:Nyb+1)-ut(1:Nxb,2:Nyb+1))))&
                         *(1/((2/(dx*dx))+(2/(dy*dy))))

#else

     do j=2,Nyb+1
        do i=2,Nxb+1
              
            p(i,j)=(((p_old(i,j+1)+p(i,j-1))/(dy*dy))&
                   +((p_old(i+1,j)+p(i-1,j))/(dx*dx))&
                   -(1/(dy*dt))*(vt(i,j)-vt(i,j-1))&
                   -(1/(dx*dt))*(ut(i,j)-ut(i-1,j)))&
                   *(1/((2/(dx*dx))+(2/(dy*dy))))*omega + (1-omega)*p(i,j)

        end do
     end do

#endif

     ! Pressure BC

     p(:,1)=p(:,2)
     p(:,Nyb+2)=p(:,Nyb+1)

     p(1,:)=p(2,:)
     p(Nxb+2,:)=p(Nxb+1,:)


     p_counter = p_counter + 1

     do i=1,Nxb+2
          p_res = p_res + sum((p(i,:)-p_old(i,:))**2)
     enddo
     
     p_res = sqrt(p_res/((Nxb+2)*(Nyb+2)))

     if( (p_res .lt. 0.000001 ) .and. (p_res .ne. 0) ) exit

  end do

end subroutine Poisson_solver
