subroutine Driver_init()
 
      use Driver_data
      use Grid_data

      implicit none

      dr_t  = 25.0
      dr_dt = 0.03*min(minval(gr_dx_nodes),minval(gr_dy_nodes))
      dr_nt = dr_t/dr_dt

end subroutine Driver_init
