subroutine IO_display()

         use Grid_data

         implicit none

         print *,"Size of Domain ",Lx," X ",Ly
         print *,"Relaxation Factor ",omega

end subroutine IO_display
