subroutine IO_write(x,y,uu,vv)

#include "Solver.h"

       implicit none
       real, dimension(Nxb+1, Nyb+1), intent(in) :: x,y,uu,vv

       open(unit = 1, file = "X.dat")
       open(unit = 2, file = "Y.dat")
       open(unit = 3, file = "U.dat")
       open(unit = 4, file = "V.dat")

       write(1,*),x
       write(2,*),y
       write(3,*),uu
       write(4,*),vv

       close(1)
       close(2)
       close(3)
       close(4)

end subroutine IO_write

