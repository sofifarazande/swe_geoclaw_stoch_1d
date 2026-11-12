! =====================================================================
!  rp1 -- coupled shallow-water / Exner / activity Riemann solver
!
!  State:  q = (h, hu, z, b)   depth, discharge, bed elevation, particle activity
!  aux(1) = z_base (rigid floor, used by src1), aux(2) = capacity (=1)
!
!  One HLL solve with shared speeds SL, SR advances h, hu and z together.
!  With wet states and sediment on, SL/SR include the roots of the
!  coupled characteristic cubic. b is advected by its own upwind flux u_p*b.
!
!  Returned as 2 f-waves (use_fwaves = True in setrun.py):
!     fwave(:,1) = F* - F(qL)  at speed SL,   fwave(:,2) = F(qR) - F*  at speed SR
!  with F = (hu, hu^2 + g h^2/2, qs/(1-porosity), u_p b).
! =====================================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

    use exner_params_module, only: g, h_wet, c_dry, porosity, compute_velocity, sediment_on
    use sediment_module, only: qs_equilibrium, qs_stochastic, particle_velocity

    implicit none

    integer, intent(in) :: maxmx,meqn,mwaves,mbc,mx,maux
    double precision, intent(in),  dimension(meqn, 1-mbc:maxmx+mbc) :: ql,qr
    double precision, intent(in),  dimension(maux, 1-mbc:maxmx+mbc) :: auxl,auxr
    double precision, intent(out) :: s(mwaves, 1-mbc:maxmx+mbc)
    double precision, intent(out) :: fwave(meqn, mwaves, 1-mbc:maxmx+mbc)
    double precision, intent(out), dimension(meqn, 1-mbc:maxmx+mbc) :: amdq,apdq

    ! dx is not an rp1 argument; claw1.f provides it in /comxt/
    double precision :: dtcom,dxcom,tcom
    common /comxt/ dtcom,dxcom,tcom

    integer :: i,m,mw
    logical :: wL,wR,ok
    double precision :: hL,huL,zL,bL,hR,huR,zR,bR
    double precision :: uL,uR,aL,aR,SL,SR
    double precision :: hRoe,uRoe,aRoe,Qu,eps_u,lam(3)
    double precision :: FhL,FhuL,FhR,FhuR,Fh_s,Fhu_s
    double precision :: slope_c(1-mbc:mx+mbc),u_pL,u_pR
    double precision :: qsL,qsR,Fqs_s,FbL,FbR,Fb_s
    double precision, parameter :: inv_solid = 1.0d0/(1.0d0 - porosity)

    ! Cell-centred slope -dz/dx (as in src1), so that each cell has one
    ! sediment flux at both of its faces and the f-wave update stays
    ! conservative in z and b.
    do i = 2-mbc, mx+mbc-1
        slope_c(i) = -(ql(3,i+1) - ql(3,i-1)) / (2.0d0*dxcom)
    end do
    slope_c(1-mbc)  = slope_c(2-mbc)
    slope_c(mx+mbc) = slope_c(mx+mbc-1)

    ! Interface i lies between cells i-1 (left) and i (right).
    do i = 2-mbc, mx+mbc

        hL = qr(1,i-1);  huL = qr(2,i-1);  zL = qr(3,i-1);  bL = qr(4,i-1)
        hR = ql(1,i);    huR = ql(2,i);    zR = ql(3,i);    bR = ql(4,i)

        fwave(:,:,i) = 0.d0
        wL = (hL > h_wet)
        wR = (hR > h_wet)

        if (.not. wL .and. .not. wR) then
            ! both dry: no flux
            s(1,i) = -c_dry
            s(2,i) =  c_dry
        else
            uL = compute_velocity(hL,huL)
            uR = compute_velocity(hR,huR)
            aL = merge(sqrt(g*hL), 0.0d0, wL)
            aR = merge(sqrt(g*hR), 0.0d0, wR)

            ! ---- wave-speed bounds ----
            if (wL .and. .not. wR) then
                SL = uL - aL
                SR = uL + 2.0d0*aL
            else if (.not. wL .and. wR) then
                SL = uR - 2.0d0*aR
                SR = uR + aR
            else
                ! Einfeldt bounds
                hRoe = 0.5d0*(hL+hR)
                uRoe = (sqrt(hL)*uL + sqrt(hR)*uR) / (sqrt(hL) + sqrt(hR) + 1.0d-30)
                aRoe = sqrt(g*hRoe)
                SL = min(uL-aL, uRoe-aRoe)
                SR = max(uR+aR, uRoe+aRoe)

                if (sediment_on) then
                    ! Coupled eigenvalues: lambda^3 - 2u lambda^2 - (gQu + gh - u^2) lambda
                    ! + gQu u = 0, with Qu = d(qs_eq/(1-porosity))/du by central difference.
                    eps_u = max(1.0d-4*max(abs(uRoe),0.01d0), 1.0d-6)
                    Qu = (qs_equilibrium(hRoe,uRoe+eps_u) - qs_equilibrium(hRoe,uRoe-eps_u)) &
                         * inv_solid / (2.0d0*eps_u)
                    call solve_cubic_real(-2.0d0*uRoe, -(g*Qu + g*hRoe - uRoe**2), &
                                          g*Qu*uRoe, lam, ok)
                    if (ok) then
                        SL = min(SL, minval(lam))
                        SR = max(SR, maxval(lam))
                    end if
                end if
            end if

            ! ---- sediment fluxes; also widen [SL,SR] to contain b's speed ----
            if (sediment_on) then
                u_pL = particle_velocity(max(hL,h_wet), uL, slope_c(i-1))
                u_pR = particle_velocity(max(hR,h_wet), uR, slope_c(i))
                qsL  = qs_stochastic(hL, uL, bL, u_pL) * inv_solid
                qsR  = qs_stochastic(hR, uR, bR, u_pR) * inv_solid

                SL = min(SL, sign(u_pL,uL), sign(u_pR,uR))
                SR = max(SR, sign(u_pL,uL), sign(u_pR,uR))

                ! b: upwind on the interface-averaged flow direction
                FbL = u_pL * max(bL,0.0d0)
                FbR = u_pR * max(bR,0.0d0)
                Fb_s = merge(FbL, FbR, 0.5d0*(uL+uR) >= 0.0d0)
            else
                qsL = 0.0d0;  qsR = 0.0d0
            end if

            ! ---- HLL interface flux ----
            FhL  = hL*uL
            FhuL = hL*uL**2 + 0.5d0*g*hL**2
            FhR  = hR*uR
            FhuR = hR*uR**2 + 0.5d0*g*hR**2

            if (SL >= 0.0d0) then
                Fh_s = FhL;  Fhu_s = FhuL;  Fqs_s = qsL
            else if (SR <= 0.0d0) then
                Fh_s = FhR;  Fhu_s = FhuR;  Fqs_s = qsR
            else
                Fh_s  = (SR*FhL  - SL*FhR  + SL*SR*(hR    - hL))    / (SR-SL)
                Fhu_s = (SR*FhuL - SL*FhuR + SL*SR*(hR*uR - hL*uL)) / (SR-SL)
                Fqs_s = (SR*qsL  - SL*qsR  + SL*SR*(zR    - zL))    / (SR-SL)
            end if

            s(1,i) = SL
            s(2,i) = SR

            fwave(1,1,i) = Fh_s  - FhL
            fwave(2,1,i) = Fhu_s - FhuL
            fwave(1,2,i) = FhR   - Fh_s
            fwave(2,2,i) = FhuR  - Fhu_s

            ! with sediment off the bed and b stay frozen (z/b f-waves remain zero)
            if (sediment_on) then
                fwave(3,1,i) = Fqs_s - qsL
                fwave(4,1,i) = Fb_s  - FbL
                fwave(3,2,i) = qsR   - Fqs_s
                fwave(4,2,i) = FbR   - Fb_s
            end if
        end if

        ! ---- split into left-/right-going flux differences ----
        do m = 1,meqn
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw = 1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + fwave(m,mw,i)
                else if (s(mw,i) > 0.d0) then
                    apdq(m,i) = apdq(m,i) + fwave(m,mw,i)
                else
                    amdq(m,i) = amdq(m,i) + 0.5d0*fwave(m,mw,i)
                    apdq(m,i) = apdq(m,i) + 0.5d0*fwave(m,mw,i)
                end if
            end do
        end do

    end do

end subroutine rp1


! =====================================================================
!  Real roots of lambda^3 + b2 lambda^2 + b1 lambda + b0 = 0 (Viete's
!  trigonometric method). ok = .false. unless all three roots are real
!  and distinct; rp1 then keeps the plain shallow-water bounds.
! =====================================================================
subroutine solve_cubic_real(b2,b1,b0,roots,ok)
    use exner_params_module, only: pi
    implicit none
    double precision, intent(in)  :: b2,b1,b0
    double precision, intent(out) :: roots(3)
    logical, intent(out) :: ok
    double precision :: p,q,disc,r,theta,sqrt_term
    integer :: k

    p = b1 - b2**2/3.0d0
    q = 2.0d0*b2**3/27.0d0 - b2*b1/3.0d0 + b0
    disc = 4.0d0*p**3 + 27.0d0*q**2

    if (p >= 0.0d0 .or. disc > 0.0d0) then
        ok = .false.
        roots = 0.0d0
        return
    end if

    sqrt_term = 2.0d0*sqrt(-p/3.0d0)
    r = (3.0d0*q)/(2.0d0*p) * sqrt(-3.0d0/p)
    r = max(-1.0d0, min(1.0d0, r))
    theta = acos(r)

    do k = 0,2
        roots(k+1) = sqrt_term * cos(theta/3.0d0 - 2.0d0*pi*dble(k)/3.0d0) - b2/3.0d0
    end do
    ok = .true.
end subroutine solve_cubic_real
