module IncompNS_interface

       implicit none

       interface
           subroutine IncompNS_init()
           end subroutine IncompNS_init
       end interface

       interface
           subroutine IncompNS_rk3(tstep,p_counter)
            implicit none
            integer, intent(in) :: tstep
            integer, intent(out) :: p_counter
           end subroutine IncompNS_rk3
       end interface

end module IncompNS_interface
