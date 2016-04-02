module Grid_data

#include "Solver.h"

       implicit none

       real, save :: Lx
       real, save :: Ly
       
       integer, save :: Nx
       integer, save :: Ny

       real, save :: dx
       real, save :: dy

       real, save, dimension(Nxb+1,Nyb+1) :: dx_n, dy_n, dx_c, dy_c
      
       real, save, dimension(Nxb,Nyb) :: dx_a, dx_b,  dy_a, dy_b     
      
       real, save :: inRe

       real, save, dimension(Nxb+1,Nyb+1) :: x
       real, save, dimension(Nxb+1,Nyb+1) :: y

       real, save :: t
       real, save :: dt

       integer, save :: nt


end module Grid_data 
