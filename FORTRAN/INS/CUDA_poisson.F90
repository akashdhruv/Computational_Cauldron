attributes(global) subroutine CUDA_poisson(p,p_old,ut,vt,dx,dy,dt,omega)

    implicit none
    real, intent(inout) :: p(:,:)
    real, intent(in) :: p_old(:,:), ut(:,:), vt(:,:) 
    real, value :: dx, dy, dt, omega

    integer :: i,j,n(2)

    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1

    n(1) = size(p,1)
    n(2) = size(p,2)

    if (i<n(1) .and. j<n(2)) then

       p(i,j)=(((p_old(i,j+1)+p(i,j-1))/(dy*dy))&
                   +((p_old(i+1,j)+p(i-1,j))/(dx*dx))&
                   -(1/(dy*dt))*(vt(i,j)-vt(i,j-1))&
                   -(1/(dx*dt))*(ut(i,j)-ut(i-1,j)))&
                   *(1/((2/(dx*dx))+(2/(dy*dy))))*omega &
                   + (1-omega)*p(i,j)

    end if

end subroutine CUDA_poisson
