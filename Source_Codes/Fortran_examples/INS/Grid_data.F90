module Grid_data

#include "Solver.h"

       implicit none

       real, save :: Lx
       real, save :: Ly

       integer, save :: iProcs
       integer, save :: jProcs
       integer, save :: tProcs
       
       integer, save :: Nx
       integer, save :: Ny

       real, save :: dx
       real, save :: dy

       real, save :: inRe

       real, save, dimension(Nxb+1,Nyb+1) :: x
       real, save, dimension(Nxb+1,Nyb+1) :: y

       real, save :: t
       real, save :: dt

       integer, save :: nt


end module Grid_data 
