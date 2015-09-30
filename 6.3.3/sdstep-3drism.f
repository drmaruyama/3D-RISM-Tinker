c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1998 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sdstep  --  Verlet stochastic dynamics step  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sdstep" performs a single stochastic dynamics time step
c     via the velocity Verlet integration algorithm
c
c     literature references:
c
c     M. P. Allen, "Brownian Dynamics Simulation of a Chemical
c     Reaction in Solution", Molecular Physics, 40, 1073-1087 (1980)
c
c     F. Guarnieri and W. C. Still, "A Rapidly Convergent Simulation
c     Method: Mixed Monte Carlo/Stochastic Dynamics", Journal of
c     Computational Chemistry, 15, 1302-1310 (1994)
c
c
      subroutine sdstep3d (istep,dt,num,gpu)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'freeze.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,istep,num,gpu
      real*8 dt,term
      real*8 epot,etot,e3d
      real*8 eksum
      real*8 temp,pres
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: pfric(:)
      real*8, allocatable :: vfric(:)
      real*8, allocatable :: afric(:)
      real*8, allocatable :: prand(:,:)
      real*8, allocatable :: vrand(:,:)
      real*8, allocatable :: derivs(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (pfric(n))
      allocate (vfric(n))
      allocate (afric(n))
      allocate (prand(3,n))
      allocate (vrand(3,n))
      allocate (derivs(3,n))
c
c     get frictional and random terms for position and velocity
c
      call sdterm (istep,dt,pfric,vfric,afric,prand,vrand)
c
c     store the current atom positions, then find full-step
c     positions and half-step velocities via modified Verlet
c
      do i = 1, n
         if (use(i)) then
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            x(i) = x(i) + v(1,i)*vfric(i) + a(1,i)*afric(i) + prand(1,i)
            y(i) = y(i) + v(2,i)*vfric(i) + a(2,i)*afric(i) + prand(2,i)
            z(i) = z(i) + v(3,i)*vfric(i) + a(3,i)*afric(i) + prand(3,i)
            do j = 1, 3
               v(j,i) = v(j,i)*pfric(i) + 0.5d0*a(j,i)*vfric(i)
            end do
         end if
      end do
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)
c
c     get the potential energy and atomic forces
c
      call gradient3d (epot,e3d,derivs,num,gpu)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using modified Verlet
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + 0.5d0*a(j,i)*vfric(i) + vrand(j,i)
            end do
         end if
      end do
c
c     correct internal virial to account for frictional forces
c
      do i = 1, n
         if (use(i)) then
            term = vfric(i)/dt - 1.0d0
            vxx = term * x(i) * derivs(1,i)
            vyx = 0.5d0 * term * (y(i)*derivs(1,i)+x(i)*derivs(2,i))
            vzx = 0.5d0 * term * (z(i)*derivs(1,i)+x(i)*derivs(3,i))
            vyy = term * y(i) * derivs(2,i)
            vzy = 0.5d0 * term * (z(i)*derivs(2,i)+y(i)*derivs(3,i))
            vzz = term * z(i) * derivs(3,i)
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (pfric)
      deallocate (vfric)
      deallocate (afric)
      deallocate (prand)
      deallocate (vrand)
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     compute and control the temperature and pressure
c
      call kinetic (eksum,ekin)
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot + e3d
c
c     compute statistics and save trajectory for this step
c
      call mdstat3d (istep,dt,etot,epot,e3d,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
      return
      end
