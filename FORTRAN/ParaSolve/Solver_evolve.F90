subroutine Solver_evolve       

#include "Solver.h"

    use IncompNS_interface, only: IncompNS_rk3
    use Grid_data
    use IncompNS_data
    use MPI_data
    use IO_interface, only: IO_display, IO_write

    implicit none

    integer :: tstep,p_counter

    real, dimension(Nxb+1,Nyb+1) :: uu
    real, dimension(Nxb+1,Nyb+1) :: vv
    real, dimension(Nxb+1,Nyb+1) :: pp

    real, pointer, dimension(:,:) :: u,v,p  
  
    tstep = 0

    do while(tstep<gr_nt) 


       call IncompNS_rk3(tstep,p_counter)

       if (mod(tstep,5) == 0 .and. myid == 0) then
          call IO_display(ins_u_res,ins_v_res,ins_p_res,p_counter,tstep*gr_dt)
       end if

       if((ins_u_res .lt. 0.0000001) .and. (ins_u_res .ne. 0).and. (ins_v_res .lt. 0.0000001) .and. (ins_v_res .ne. 0) ) exit

     tstep = tstep +1

    end do

    u => gr_facex(VELC_VAR,:,:)
    v => gr_facey(VELC_VAR,:,:)
    p => gr_center(PRES_VAR,:,:)

    uu = ((u(1:Nxb+1,1:Nyb+1)+u(1:Nxb+1,2:Nyb+2))/2 + (u(2:Nxb+2,1:Nyb+1)+u(2:Nxb+2,2:Nyb+2))/2)/2
    vv = ((v(1:Nxb+1,1:Nyb+1)+v(2:Nxb+2,1:Nyb+1))/2 + (v(1:Nxb+1,2:Nyb+2)+v(2:Nxb+2,2:Nyb+2))/2)/2
    pp = ((p(1:Nxb+1,1:Nyb+1)+p(2:Nxb+2,1:Nyb+1))/2 + (p(1:Nxb+1,2:Nyb+2)+p(2:Nxb+2,2:Nyb+2))/2)/2

    nullify(u)
    nullify(v)
    nullify(p)

    call IO_write(gr_x,gr_y,uu,vv,pp,myid)

  
end subroutine Solver_evolve
