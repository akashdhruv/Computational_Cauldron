subroutine Poisson_solver(ps_RHS,ps,ps_res,ps_counter,ps_quant)

  use Grid_data
  use MPI_data
  use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_pres

#include "Solver.h"
                
  implicit none

  real, dimension(Nxb,Nyb), intent(in) :: ps_RHS

  real, dimension(Nxb+2,Nyb+2), intent(inout) :: ps

  integer, intent(in) :: ps_quant

  real, dimension(Nxb+2,Nyb+2) :: ps_old

  real, intent(out) :: ps_res
        
  real :: ps_res1

  integer, intent(out) :: ps_counter
  integer :: i,j

  ps_old = 0
  ps_counter = 0
 
  do while(ps_counter<MaxIt)

     ps_res = 0  
     ps_res1 = 0     
     ps_old = ps

#ifdef POISSON_SOLVER_JACOBI

     ps(2:Nxb+1,2:Nyb+1)= ((ps_old(2:Nxb+1,3:Nyb+2)/(gr_dy_centers(2:Nxb+1,2:Nyb+1)*gr_dy_nodes(2:Nxb+1,2:Nyb+1)))&
                          +(ps_old(2:Nxb+1,1:Nyb)/(gr_dy_centers(1:Nxb,1:Nyb)*gr_dy_nodes(2:Nxb+1,2:Nyb+1)))&
                          +(ps_old(3:Nxb+2,2:Nyb+1)/(gr_dx_centers(2:Nxb+1,2:Nyb+1)*gr_dx_nodes(2:Nxb+1,2:Nyb+1)))&
                          +(ps_old(1:Nxb,2:Nyb+1)/(gr_dx_centers(1:Nxb,1:Nyb)*gr_dx_nodes(2:Nxb+1,2:Nyb+1)))&
                          +ps_RHS)&
                          *(1/((1/(gr_dx_nodes(2:Nxb+1,2:Nyb+1)*gr_dx_centers(1:Nxb,1:Nyb)))&
                          +(1/(gr_dy_nodes(2:Nxb+1,2:Nyb+1)*gr_dy_centers(1:Nxb,1:Nyb)))&
                          +(1/(gr_dx_nodes(2:Nxb+1,2:Nyb+1)*gr_dx_centers(2:Nxb+1,2:Nyb+1)))&
                          +(1/(gr_dy_nodes(2:Nxb+1,2:Nyb+1)*gr_dy_centers(2:Nxb+1,2:Nyb+1)))))

#endif

#ifdef POISSON_SOLVER_GS

     do j=2,Nyb+1
        do i=2,Nxb+1

           ps(i,j)=((ps_old(i,j+1)/(gr_dy_centers(i,j)*gr_dy_nodes(i,j)))+(ps(i,j-1)/(gr_dy_nodes(i,j)*gr_dy_centers(i-1,j-1)))&
                  +(ps_old(i+1,j)/(gr_dx_centers(i,j)*gr_dx_nodes(i,j)))+(ps(i-1,j)/(gr_dx_nodes(i,j)*gr_dx_centers(i-1,j-1)))&
                  +ps_RHS(i-1,j-1))&
                  *(1/((1/(gr_dx_nodes(i,j)*gr_dx_centers(i-1,j-1)))+(1/(gr_dy_nodes(i,j)*gr_dy_centers(i-1,j-1)))+&
                   (1/(gr_dx_nodes(i,j)*gr_dx_centers(i,j)))+(1/(gr_dy_nodes(i,j)*gr_dy_centers(i,j)))))

        end do
     end do


#endif

#ifdef POISSON_SOLVER_GSOR

     do j=2,Nyb+1
        do i=2,Nxb+1

           ps(i,j)=((ps_old(i,j+1)/(gr_dy_centers(i,j)*gr_dy_nodes(i,j)))+(ps(i,j-1)/(gr_dy_nodes(i,j)*gr_dy_centers(i-1,j-1)))&
                  +(ps_old(i+1,j)/(gr_dx_centers(i,j)*gr_dx_nodes(i,j)))+(ps(i-1,j)/(gr_dx_nodes(i,j)*gr_dx_centers(i-1,j-1)))&
                  +ps_RHS(i-1,j-1))&
                  *(1/((1/(gr_dx_nodes(i,j)*gr_dx_centers(i-1,j-1)))+(1/(gr_dy_nodes(i,j)*gr_dy_centers(i-1,j-1)))+&
                   (1/(gr_dx_nodes(i,j)*gr_dx_centers(i,j)))+(1/(gr_dy_nodes(i,j)*gr_dy_centers(i,j)))))*omega + (1-omega)*ps(i,j)
                  
        end do
     end do

#endif

     ! Pressure BC

     call MPI_applyBC(ps)

     if(ps_quant == PRES_VAR) call MPI_physicalBC_pres(ps)
 
     ps_counter = ps_counter + 1

     do i=1,Nxb+2
          ps_res = ps_res + sum((ps(i,:)-ps_old(i,:))**2)
     enddo
     
     call MPI_CollectResiduals(ps_res,ps_res1)

     ps_res = sqrt(ps_res1/((Nxb+2)*(Nyb+2)*(HK**HD)))

     if( (ps_res .lt. 0.000001 ) .and. (ps_res .ne. 0) ) exit

  end do

end subroutine Poisson_solver
