! =====================================================================
!  qinit -- lake-at-rest initial conditions: flat free surface
!  eta = h + z = eta_const over the bump bed aux(1), zero velocity, b = 0.
! =====================================================================
subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    implicit none
    integer, intent(in) :: meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx
    double precision, intent(in) :: aux(maux,1-mbc:mx+mbc)
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc)

    double precision, parameter :: eta_const = 0.35d0
    integer :: i

    do i = 1,mx
        q(1,i) = eta_const - aux(1,i)   ! h, so that h+z = eta_const everywhere
        q(2,i) = 0.0d0                  ! hu
        q(3,i) = aux(1,i)               ! z = rigid floor (no extra sediment layer)
        q(4,i) = 0.0d0                  ! b
    end do

end subroutine qinit
