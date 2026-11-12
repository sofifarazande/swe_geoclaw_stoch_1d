! =====================================================================
!  setprob -- one-time setup, called by claw1ez after set_grid.
!
!  Reads setprob.data (written by setrun.py, fixed order):
!     sediment_on, friction_on, mu_over_sigma, cd_dep, bed_alpha,
!     then optionally mann_n_flow, sgn_on, sgn_bed_sigma (defaults kept
!     if absent, e.g. in the tests' 5-value files)
!  then initialises the sediment closure, which needs mu_over_sigma and
!  cd_dep already set, and dx for the reaction block size.
! =====================================================================
subroutine setprob

    use grid_module, only: xp_edge
    use sediment_module, only: init_stochastic_transport, mu_over_sigma, cd_dep
    use exner_params_module, only: sediment_on, friction_on, bed_alpha, mann_n_flow
    use sgn_module, only: sgn_on, sgn_bed_sigma
    implicit none

    integer, parameter :: iunit = 17
    integer :: sediment_on_flag, friction_on_flag, sgn_on_flag, ios

    call opendatafile(iunit, 'setprob.data')
    read(iunit,*) sediment_on_flag
    read(iunit,*) friction_on_flag
    read(iunit,*) mu_over_sigma
    read(iunit,*) cd_dep
    read(iunit,*) bed_alpha
    read(iunit,*,iostat=ios) mann_n_flow
    if (ios == 0) read(iunit,*,iostat=ios) sgn_on_flag
    if (ios == 0) sgn_on = (sgn_on_flag /= 0)
    if (ios == 0) read(iunit,*,iostat=ios) sgn_bed_sigma
    close(iunit)
    sediment_on = (sediment_on_flag /= 0)
    friction_on = (friction_on_flag /= 0)

    call init_stochastic_transport(xp_edge(2) - xp_edge(1))

end subroutine setprob
