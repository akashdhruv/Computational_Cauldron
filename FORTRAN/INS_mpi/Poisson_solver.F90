subroutine Poisson_solver(ut,vt,p_res,p_counter)

  use IncompNS_data
  use Grid_data 
  use MPI_data
  use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_pres

#include "Solver.h"
                
  implicit none

  real, dimension(Nxb+2,Nyb+2), intent(in) :: ut
  real, dimension(Nxb+2,Nyb+2), intent(in) :: vt

  real, dimension(Nxb+2,Nyb+2) :: p_old

  real, intent(out) :: p_res
        
  real :: p_res1

  integer, intent(out) :: p_counter
  integer :: i,j

  p_old = 0
  p_counter = 0
 
  do while(p_counter<MaxIt)

     p_res = 0  
     p_res1 = 0     
     p_old = p

! Warning - Jacobi not working for Non Uniform Grid. Debugging required

#ifdef POISSON_SOLVER_JACOBI

     p(2:Nxb+1,2:Nyb+1)= ((p_old(2:Nxb+1,3:Nyb+2)/(dy_centers(2:Nxb+1,2:Nyb+1)*dy_nodes(2:Nxb+1,2:Nyb+1)))&
                          +(p_old(2:Nxb+1,1:Nyb)/(dy_centers(1:Nxb,1:Nyb)*dy_nodes(2:Nxb+1,2:Nyb+1)))&
                         +(p_old(3:Nxb+2,2:Nyb+1)/(dx_centers(2:Nxb+1,2:Nyb+1)*dy_nodes(2:Nxb+1,2:Nyb+1)))&
                          +(p_old(1:Nxb,2:Nyb+1)/(dx_centers(1:Nxb,1:Nyb)*dx_nodes(2:Nxb+1,2:Nyb+1)))&
                         -((1/(dy_nodes(2:Nxb+1,2:Nyb+1)*dt))*(vt(2:Nxb+1,2:Nyb+1)-vt(2:Nxb+1,1:Nyb)))&
                         -((1/(dx_nodes(2:Nxb+1,2:Nyb+1)*dt))*(ut(2:Nxb+1,2:Nyb+1)-ut(1:Nxb,2:Nyb+1))))&
                         *(1/((1/(dx_nodes(2:Nxb+1,2:Nyb+1)*dx_centers(1:Nxb,1:Nyb)))&
                         +(1/(dy_nodes(2:Nxb+1,2:Nyb+1)*dy_centers(1:Nxb,1:Nyb)))&
                         +(1/(dx_centers(2:Nxb+1,2:Nyb+1)*dx_nodes(2:Nxb+1,2:Nyb+1)))&
                         +(1/(dy_centers(2:Nxb+1,2:Nyb+1)*dy_nodes(2:Nxb+1,2:Nyb+1)))))

#endif

#ifdef POISSON_SOLVER_GS

     do j=2,Nyb+1
        do i=2,Nxb+1

           p(i,j)=((p_old(i,j+1)/(dy_centers(i,j)*dy_nodes(i,j)))+(p(i,j-1)/(dy_nodes(i,j)*dy_centers(i-1,j-1)))&
                  +(p_old(i+1,j)/(dx_centers(i,j)*dx_nodes(i,j)))+(p(i-1,j)/(dx_nodes(i,j)*dx_centers(i-1,j-1)))&
                  -((1/(dy_nodes(i,j)*dt))*(vt(i,j)-vt(i,j-1)))&
                  -((1/(dx_nodes(i,j)*dt))*(ut(i,j)-ut(i-1,j))))&
                  *(1/((1/(dx_nodes(i,j)*dx_centers(i-1,j-1)))+(1/(dy_nodes(i,j)*dy_centers(i-1,j-1)))+&
                   (1/(dx_centers(i,j)*dx_nodes(i,j)))+(1/(dy_centers(i,j)*dy_nodes(i,j)))))

        end do
     end do


#endif

#ifdef POISSON_SOLVER_GSOR

     do j=2,Nyb+1
        do i=2,Nxb+1

           p(i,j)=((p_old(i,j+1)/(dy_centers(i,j)*dy_nodes(i,j)))+(p(i,j-1)/(dy_nodes(i,j)*dy_centers(i-1,j-1)))&
                  +(p_old(i+1,j)/(dx_centers(i,j)*dx_nodes(i,j)))+(p(i-1,j)/(dx_nodes(i,j)*dx_centers(i-1,j-1)))&
                  -((1/(dy_nodes(i,j)*dt))*(vt(i,j)-vt(i,j-1)))&
                  -((1/(dx_nodes(i,j)*dt))*(ut(i,j)-ut(i-1,j))))&
                  *(1/((1/(dx_nodes(i,j)*dx_centers(i-1,j-1)))+(1/(dy_nodes(i,j)*dy_centers(i-1,j-1)))+&
                   (1/(dx_centers(i,j)*dx_nodes(i,j)))+(1/(dy_centers(i,j)*dy_nodes(i,j)))))*omega + (1-omega)*p(i,j)
                  
        end do
     end do

#endif

#ifdef POISSON_SOLVER_FFT

! FFT Solver 

#endif

     ! Pressure BC

     call MPI_applyBC(p)
     call MPI_physicalBC_pres(p)

     p_counter = p_counter + 1

     do i=1,Nxb+2
          p_res = p_res + sum((p(i,:)-p_old(i,:))**2)
     enddo
     
     call MPI_CollectResiduals(p_res,p_res1)

     p_res = sqrt(p_res1/((Nxb+2)*(Nyb+2)*(HK**HD)))

     if( (p_res .lt. 0.000001 ) .and. (p_res .ne. 0) ) exit

  end do

end subroutine Poisson_solver
