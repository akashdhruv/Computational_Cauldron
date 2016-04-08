module Grid_data

#include "Solver.h"

       implicit none

       real, save :: Lx
       real, save :: Ly
       
       integer, save :: Nx
       integer, save :: Ny

       real, save :: dx
       real, save :: dy

       real, save, dimension(Nxb+1,Nyb+1) ::  dx_centers, dy_centers

       real, save, dimension(Nxb+2,Nyb+2) :: dx_nodes, dy_nodes 
 
       real, save :: inRe

       real, save, dimension(Nxb+1,Nyb+1) :: x
       real, save, dimension(Nxb+1,Nyb+1) :: y

       real, save :: t
       real, save :: dt

       integer, save :: nt

end module Grid_data 
