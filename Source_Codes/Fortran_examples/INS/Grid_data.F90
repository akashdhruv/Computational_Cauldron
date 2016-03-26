module Grid_data

#include "Solver.h"

       implicit none

       real, save :: Lx
       real, save :: Ly 

       integer, save :: Nxb
       integer, save :: Nyb

       integer, save :: iProcs
       integer, save :: jProcs
       integer, save :: tProcs
       
       integer, save :: Nx
       integer, save :: Ny

       real, save :: dx
       real, save :: dy

       real, save :: inRe
       real, save :: omega

end module Grid_data 
