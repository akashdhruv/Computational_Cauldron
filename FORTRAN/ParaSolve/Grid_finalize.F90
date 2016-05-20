subroutine Grid_finalize

   use Grid_data

   implicit none
   
   deallocate(gr_dx_centers)
   deallocate(gr_dy_centers)
   deallocate(gr_dx_nodes)
   deallocate(gr_dy_nodes)
   deallocate(gr_x)
   deallocate(gr_y)
   deallocate(gr_center)
   deallocate(gr_facex)
   deallocate(gr_facey)

end subroutine Grid_finalize 
