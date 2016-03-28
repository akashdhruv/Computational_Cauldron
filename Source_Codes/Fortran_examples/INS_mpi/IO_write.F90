subroutine IO_write(x,y,uu,vv,id)

#include "Solver.h"

       implicit none
       real, dimension(Nxb+1, Nyb+1), intent(in) :: x,y,uu,vv
       integer, intent(in) :: id
       character(len=10) :: f1,f2,f3,f4
     
       write (f1, '( "X", I1, ".dat" )' )id
       write (f2, '( "Y", I1, ".dat" )' )id
       write (f3, '( "U", I1, ".dat" )' )id
       write (f4, '( "V", I1, ".dat" )' )id

       open(unit = 1, file = f1)
       open(unit = 2, file = f2)
       open(unit = 3, file = f3)
       open(unit = 4, file = f4)

       write(1,*),x
       write(2,*),y
       write(3,*),uu
       write(4,*),vv

       close(1)
       close(2)
       close(3)
       close(4)

end subroutine IO_write

