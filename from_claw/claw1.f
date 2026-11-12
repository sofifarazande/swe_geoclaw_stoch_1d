c
c     ==============================================================
      subroutine claw1(meqn,mwaves,maux,mbc,mx,
     &           q,aux,xlower,dx,tstart,tend,dtv,cflv,nv,method,mthlim,
     &           mthbc,work,mwork,use_fwave,info,bc1,rp1,src1,b4step1)
c
c  Advances q from tstart to tend for  capa*q_t + A q_x = psi  with the
c  wave-propagation method (step1, which calls rp1). Each step calls
c  bc1, b4step1, then step1 and src1.
c
c  method(1) = 1  variable dt from the desired Courant number cflv(2);
c                 a step with Courant number > cflv(1) is redone
c  method(2) = 2  second-order corrections, limiters mthlim(1:mwaves)
c  method(4) = 1  print dt and Courant number every step
c  method(5) = 1  Godunov-split source (src1 after step1),
c            = 2  Strang split (src1 for dt/2 before and after)
c  method(6) = index of the capacity function in aux (0 = none)
c
c  On return: tend = time reached, dtv(3:5) = min/max/last dt,
c  cflv(3:4) = max/last Courant number, nv(2) = number of steps,
c  info = 0 on success, otherwise the error code set below.
c
c  Copyright 1994-2002 R. J. LeVeque. Clawpack is distributed under the
c  BSD 3-Clause License, see $CLAW/LICENSE.
c

      use gauges_module, only: num_gauges, update_gauges,
     &                         print_gauges_and_reset_nextLoc
      use topo_module, only: topo_finalized, topo_update, dt_max_dtopo

      implicit double precision (a-h,o-z)
      external bc1,rp1,src1,b4step1
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      dimension work(mwork)
      dimension mthlim(mwaves),method(7),dtv(5),cflv(4),nv(2)
      dimension mthbc(2)
      common /comxt/ dtcom,dxcom,tcom
c
c
      info = 0
      t = tstart
      maxn = nv(1)
      dt = dtv(1)   !# initial dt
      cflmax = 0.d0
      dtmin = dt
      dtmax = dt
      nv(2) = 0
c
c     # check for errors in data:
c
c
      if (method(1) .eq. 0) then
c        # fixed size time steps.  Compute the number of steps:
         if (tend .lt. tstart) then
c             # single step mode
              maxn = 1
           else
              maxn = nint((tend - tstart)/dt)    ! Round to nearest int
              if (dabs(maxn*dt - (tend-tstart)) .gt.
     &                          1d-5*(tend-tstart)) then
c                # dt doesn't divide time interval integer number of times
                 info = 2
                 go to 900
                 endif
           endif
         endif
c
      if (method(1).eq.1 .and. cflv(2).gt.cflv(1)) then
         info = 3
         go to 900
         endif
c
c     # partition work array into pieces for passing into step1:
      i0f = 1
      i0wave = i0f + (mx + 2*mbc) * meqn
      i0s = i0wave + (mx + 2*mbc) * meqn * mwaves
      i0dtdx = i0s + (mx + 2*mbc) * mwaves
      i0qwork = i0dtdx + (mx + 2*mbc)
      i0amdq = i0qwork + (mx + 2*mbc) * meqn
      i0apdq = i0amdq + (mx + 2*mbc) * meqn
      i0dtdx = i0apdq + (mx + 2*mbc) * meqn
      i0end = i0dtdx + (mx + 2*mbc) - 1
c
      if (mwork .lt. i0end) then
         write(6,*) 'mwork must be increased to ',i0end
         info = 4
         go to 900
         endif
c
c     -----------
c     # main loop
c     -----------
c
      if (maxn.eq.0) go to 900
      do 100 n=1,maxn
         told = t   !# time at beginning of time step.

c        # adjust dt to hit tend exactly if we're near end of computation
c        #  (unless tend < tstart, which is a flag to take only a single step)
         if (told+dt.gt.tend .and. tstart.lt.tend) dt = tend - told

         if (method(1).eq.1) then
c           # save old q in case we need to retake step with smaller dt:
            call copyq1(meqn,mbc,mx,q,work(i0qwork))
            endif
c
   40    continue
         dt2 = dt / 2.d0
         thalf = t + dt2  !# midpoint in time for Strang splitting
         t = told + dt    !# time at end of step

c        # store dt and t in the common block comxt in case they are needed
c        # in the Riemann solvers (for variable coefficients)
         tcom = told
         dtcom = dt
         dxcom = dx
c
c        ------------------------------------------------------------------
c        # main steps in algorithm:
c        ------------------------------------------------------------------
c
c        # extend data from grid to bordering boundary cells:
         call bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,told,dt,mthbc)
c
c
c        # call user-supplied routine which might set aux arrays
c        # for this time step, for example.

         call b4step1(mbc,mx,meqn,q,
     &                xlower,dx,told,dt,maux,aux)
c
         call update_gauges(q,aux,xlower,meqn,mx,mbc,maux,dx,told)
c
         if (method(5).eq.2) then
c            # with Strang splitting for source term:
             call src1(meqn,mbc,mx,xlower,dx,q,maux,aux,told,dt2)
             endif
c
c        # take a step on the homogeneous conservation law:
         call step1(meqn,mwaves,mbc,maux,mx,q,aux,dx,dt,
     &             method,mthlim,cfl,work(i0f),work(i0wave),
     &             work(i0s),work(i0amdq),work(i0apdq),work(i0dtdx),
     &             use_fwave,rp1)
c
         if (method(5).eq.2) then
c            # source terms over a second half time step for Strang splitting:
c            # Note it is not so clear what time t should be used here if
c            # the source terms are time-dependent!
             call src1(meqn,mbc,mx,xlower,dx,q,maux,aux,thalf,dt2)
             endif

         if (method(5).eq.1) then
c            # source terms over a full time step:
             call src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
             endif
c

c
c        ------------------------------------------------------------------
c
         if (method(4) .eq. 1) write(6,601) n,cfl,dt,t
  601    format('CLAW1... Step',i4,
     &                   '   Courant number =',f6.3,'  dt =',d12.4,
     &                   '  t =',d12.4)
c
         if (method(1) .eq. 1) then
c           # choose new time step if variable time step
            if (cfl .gt. 0.d0) then
                dt = dmin1(dtv(2), dt * cflv(2)/cfl)
                if (.not. topo_finalized) then
                    dt = dmin1(dt, dt_max_dtopo)
                endif
                dtmin = dmin1(dt,dtmin)
                dtmax = dmax1(dt,dtmax)


              else
                dt = dtv(2)
              endif
            endif
c
c        # check to see if the Courant number was too large:
c
         if (cfl .le. cflv(1)) then
c               # accept this step
                cflmax = dmax1(cfl,cflmax)
              else
c               # reject this step
                t = told
                call copyq1(meqn,mbc,mx,work(i0qwork),q)
c
                if (method(4) .eq. 1) then
                   write(6,602)  cfl
  602              format('CLAW1 rejecting step... ',
     &                         'Courant number too large: ',d18.8)
                   endif
                if (method(1).eq.1) then
c                   # if variable dt, go back and take a smaller step
                    go to 40
                  else
c                   # if fixed dt, give up and return
                    cflmax = dmax1(cfl,cflmax)
                    go to 900
                  endif
               endif

         if (.not. topo_finalized) then
             ! modify topo using dtopo arrays:
             call topo_update(t)
             call setaux(mbc,mx,xlower,dx,maux,aux)
         endif

c
c        # see if we are done:
         nv(2) = nv(2) + 1
         if (t .ge. tend) go to 900
c
  100    continue
c
  900  continue
c
c      # return information
c
       if (method(1).eq.1 .and. t.lt.tend .and. nv(2) .eq. maxn) then
c         # too many timesteps
          info = 11
          endif
c
       if (method(1).eq.0 .and. cflmax .gt. cflv(1)) then
c         # Courant number too large with fixed dt
          info = 12
          endif
       tend = t
       cflv(3) = cflmax
       cflv(4) = cfl
       dtv(3) = dtmin
       dtv(4) = dtmax
       dtv(5) = dt

       ! output gauges:
       do ii = 1, num_gauges
          call print_gauges_and_reset_nextLoc(ii,meqn)
       end do

       return
       end
