module IO_interface

       implicit none

       interface
             subroutine IO_display(u_res,v_res,p_res,p_counter)
             implicit none
             real, intent(in) :: u_res,v_res,p_res
             integer, intent(in) :: p_counter
             end subroutine IO_display
       end interface

end module IO_interface
