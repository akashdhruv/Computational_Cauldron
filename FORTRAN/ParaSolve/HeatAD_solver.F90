subroutine HeatAD_solver(tstep)

#include "Solver.h"

      use Grid_data
      use IncompNS_data
      use HeatAD_data
      use Driver_data
      use physicaldata
      use MPI_interface, only: MPI_applyBC, MPI_physicalBC_temp, MPI_CollectResiduals

      implicit none
      
      integer, intent(in) :: tstep
      real, pointer, dimension(:,:) :: T,u,v
      real, dimension(Nxb+2,Nyb+2) :: T_old

      integer :: i,j

      real :: u_plus, u_mins, v_plus, v_mins, u_conv, v_conv
      real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins

      real :: T_res1

      ht_T_res = 0.0
      T_res1 = 0.0

      T => ph_center(TEMP_VAR,:,:)
      u => ph_facex(VELC_VAR,:,:)
      v => ph_facey(VELC_VAR,:,:)

      T_old = T

#ifdef TEMP_SOLVER_CENTRAL

   T(2:Nxb+1,2:Nyb+1) = T_old(2:Nxb+1,2:Nyb+1) &
  +((dr_dt*ins_inRe)/(ht_Pr*(gr_dx_centers(2:Nxb+1,2:Nyb+1)**2)))*(T_old(3:Nxb+2,2:Nyb+1)+T_old(1:Nxb,2:Nyb+1)-2*T_old(2:Nxb+1,2:Nyb+1))&
  +((dr_dt*ins_inRe)/(ht_Pr*(gr_dy_centers(2:Nxb+1,2:Nyb+1)**2)))*(T_old(2:Nxb+1,3:Nyb+2)+T_old(2:Nxb+1,1:Nyb)-2*T_old(2:Nxb+1,2:Nyb+1))&
  -((dr_dt*(u(2:Nxb+1,2:Nyb+1) + u(1:Nxb,2:Nyb+1))/2)/(gr_dx_centers(2:Nxb+1,2:Nyb+1)+gr_dx_centers(1:Nxb,2:Nyb+1)))&
   *(T_old(3:Nxb+2,2:Nyb+1)-T_old(1:Nxb,2:Nyb+1))&
  -((dr_dt*(v(2:Nxb+1,2:Nyb+1) + v(2:Nxb+1,1:Nyb))/2)/(gr_dy_centers(2:Nxb+1,2:Nyb+1)+gr_dx_centers(2:Nxb+1,1:Nyb)))&
  *(T_old(2:Nxb+1,3:Nyb+2)-T_old(2:Nxb+1,1:Nyb))

#endif


#ifdef TEMP_SOLVER_UPWIND

  do j=2,Nxb+1
     do i=2,Nyb+1

     u_conv = (u(i,j)+u(i-1,j))/2.
     v_conv = (v(i,j)+v(i,j-1))/2.

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)

     Tx_plus = (T_old(i+1,j)-T_old(i,j))/gr_dx_centers(i,j)
     Tx_mins = (T_old(i,j)-T_old(i-1,j))/gr_dx_centers(i-1,j)

     Ty_plus = (T_old(i,j+1)-T_old(i,j))/gr_dy_centers(i,j)
     Ty_mins = (T_old(i,j)-T_old(i,j-1))/gr_dy_centers(i,j-1)

     T(i,j) = T_old(i,j)+((dr_dt*ins_inRe)/(ht_Pr*gr_dx_centers(i,j)*gr_dx_centers(i,j)))*(T_old(i+1,j)+T_old(i-1,j)-2*T_old(i,j))&
                        +((dr_dt*ins_inRe)/(ht_Pr*gr_dy_centers(i,j)*gr_dy_centers(i,j)))*(T_old(i,j+1)+T_old(i,j-1)-2*T_old(i,j))&
                        -((dr_dt))*(u_plus*Tx_mins + u_mins*Tx_plus)&
                        -((dr_dt))*(v_plus*Ty_mins + v_mins*Ty_plus)

     end do
  end do

#endif

     call MPI_applyBC(T)
     call MPI_physicalBC_temp(T)

     do i=1,Nxb+2
          ht_T_res = ht_T_res + sum((T(i,:)-T_old(i,:))**2)
     enddo

     call MPI_CollectResiduals(ht_T_res,T_res1)

     ht_T_res = sqrt(T_res1/((Nxb+2)*(Nyb+2)*(HK**HD)))

     nullify(T)
     nullify(u)
     nullify(v)

end subroutine HeatAD_solver
