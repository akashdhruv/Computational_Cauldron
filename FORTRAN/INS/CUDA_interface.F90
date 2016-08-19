module CUDA_interface

  contains
  attributes(global) subroutine CUDA_poisson(p,p_old,ut,vt,dx,dy,dt,omega)
    implicit none
    real, intent(inout) :: p(:,:)
    real, intent(in) :: p_old(:,:), ut(:,:), vt(:,:)
    real, value :: dx, dy, dt, omega
  end subroutine CUDA_poisson

end module CUDA_interface
