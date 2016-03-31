attributes(global) subroutine Poisson_kernel(p,u,v,Nx,Ny,dx,dy,dt)
             implicit none
             integer, value :: Nx, Ny
             real, value :: dx,dy,dt
             real, device :: p(Nx+2,Ny+2),u(Nx+2,Ny+2),v(Nx+2,Ny+2)
             integer :: id

end subroutine Poisson_kernel
