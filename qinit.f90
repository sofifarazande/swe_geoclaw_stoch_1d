! =====================================================================
!  qinit -- initial conditions: dry channel, b = 0, and a sediment layer
!  of thickness h_bed_init on the rigid floor aux(1) (read from
!  z_base.topo), plus an optional sinusoidal bed perturbation.
! =====================================================================
subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    use exner_params_module, only: PERTURB_ON, A_perturb, lam_perturb, h_bed_init, pi
    use grid_module, only: xcell

    implicit none
    integer, intent(in) :: meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx
    double precision, intent(in) :: aux(maux,1-mbc:mx+mbc)
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc)

    integer :: i

    do i = 1,mx
        q(1,i) = 0.0d0
        q(2,i) = 0.0d0
        q(3,i) = aux(1,i) + h_bed_init
        if (PERTURB_ON) q(3,i) = q(3,i) + A_perturb * sin(2.0d0*pi*xcell(i)/lam_perturb)
        q(4,i) = 0.0d0
    end do

end subroutine qinit
