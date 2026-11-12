! =====================================================================
!  bc1     -- ghost-cell boundary conditions
!  b4step1 -- positivity fix applied right after bc1, before each step
!
!  mthbc (setrun.py): 0 = user, 1 = extrapolation, 2 = periodic, 3 = wall
!
!  User conditions:
!   left  -- inflow: h = h_in, u = Q_in/(B h_in), z extrapolated, and
!            b = qs_in/u_p so that rp1's qs = u_p b gives exactly qs_in.
!   right -- outflow: h, z extrapolated, u extrapolated unless it is
!            negative (backflow is blocked: hu = 0), b = 0 (rp1's upwind
!            b flux then takes the interior value).
! =====================================================================
subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)

    use exner_params_module, only: h_in, Q_in, B_channel, qs_in, h_wet
    use sediment_module, only: particle_velocity

    implicit none
    integer, intent(in) :: meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx,t,dt
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc)
    double precision, intent(in) :: aux(maux,1-mbc:mx+mbc)
    integer, intent(in) :: mthbc(2)

    integer :: ibc
    double precision :: u_in, b_in, u_last

    ! ---- left ----
    select case (mthbc(1))
    case (0)
        u_in = Q_in / (B_channel * h_in)
        b_in = qs_in / particle_velocity(h_in, u_in, 0.0d0)
        do ibc = 1,mbc
            q(1,1-ibc) = h_in
            q(2,1-ibc) = h_in * u_in
            q(3,1-ibc) = q(3,1)
            q(4,1-ibc) = b_in
        end do
    case (1)
        do ibc = 1,mbc
            q(:,1-ibc) = q(:,1)
        end do
    case (2)
        do ibc = 1,mbc
            q(:,1-ibc) = q(:,mx+1-ibc)
        end do
    case (3)
        do ibc = 1,mbc
            q(:,1-ibc) = q(:,ibc)
            q(2,1-ibc) = -q(2,ibc)
        end do
    case default
        write(6,*) '*** ERROR *** bc1: invalid mthbc(1) =', mthbc(1)
        stop
    end select

    ! ---- right ----
    select case (mthbc(2))
    case (0)
        u_last = 0.0d0
        if (q(1,mx) > h_wet) u_last = q(2,mx) / q(1,mx)
        do ibc = 1,mbc
            q(1,mx+ibc) = q(1,mx)
            q(2,mx+ibc) = merge(q(1,mx) * u_last, 0.0d0, u_last >= 0.0d0)
            q(3,mx+ibc) = q(3,mx)
            q(4,mx+ibc) = 0.0d0
        end do
    case (1)
        do ibc = 1,mbc
            q(:,mx+ibc) = q(:,mx)
        end do
    case (2)
        do ibc = 1,mbc
            q(:,mx+ibc) = q(:,ibc)
        end do
    case (3)
        do ibc = 1,mbc
            q(:,mx+ibc) = q(:,mx+1-ibc)
            q(2,mx+ibc) = -q(2,mx+1-ibc)
        end do
    case default
        write(6,*) '*** ERROR *** bc1: invalid mthbc(2) =', mthbc(2)
        stop
    end select

end subroutine bc1


! Clamp h >= 0 and zero hu in near-dry cells; z and b are left untouched.
subroutine b4step1(mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    use exner_params_module, only: h_wet

    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    double precision, intent(in) :: xlower,dx,t,dt
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc)

    integer :: i

    do i = 1-mbc, mx+mbc
        if (q(1,i) <= h_wet) then
            if (q(1,i) < 0.0d0) q(1,i) = 0.0d0
            q(2,i) = 0.0d0
        end if
    end do

end subroutine b4step1
