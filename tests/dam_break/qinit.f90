! =====================================================================
!  qinit -- dam-break initial conditions: depth h0 behind a dam at
!  x = x0, dry downstream, flat bed, zero velocity, b = 0.
!  x0 and h0 must match check_dam_break.py.
! =====================================================================
subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    use grid_module, only: xcell
    implicit none
    integer, intent(in) :: meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx
    double precision, intent(in) :: aux(maux,1-mbc:mx+mbc)
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc)

    double precision, parameter :: x0 = 3.0d0
    double precision, parameter :: h0 = 0.5d0
    integer :: i

    do i = 1,mx
        if (xcell(i) < x0) then
            q(1,i) = h0
        else
            q(1,i) = 0.0d0
        end if
        q(2,i) = 0.0d0        ! hu
        q(3,i) = aux(1,i)     ! z = flat floor
        q(4,i) = 0.0d0        ! b
    end do

end subroutine qinit
