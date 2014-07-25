c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient  --  find energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine gradient3d (energy,sfe,derivs,num,gpu)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inter.i'
      include 'iounit.i'
      include 'potent.i'
      include 'rigid.i'
      include 'vdwpot.i'
      include 'virial.i'
      include 'atmtyp.i'
      include 'charge.i'
      include 'kvdws.i'
      integer i,j
      integer fu,freeunit
      integer count/0/,num,gpu
      real*8 energy,cutoff
      real*8 derivs(3,*)
      real*8 sfe
      real*8 sigfact, qu(maxatm), sigu(maxatm), epsu(maxatm)
      character*80 dummy
      character*4 gn
      save    count
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
      e3d = 0.0d0
c
c     zero out each of the first derivative components
c
      do i = 1, n
         do j = 1, 3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deaa(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            deopd(j,i) = 0.0d0
            deid(j,i) = 0.0d0
            deit(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            debt(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dec(j,i) = 0.0d0
            decd(j,i) = 0.0d0
            ded(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
            der(j,i) = 0.0d0
            des(j,i) = 0.0d0
            delf(j,i) = 0.0d0
            deg(j,i) = 0.0d0
            dex(j,i) = 0.0d0
            de3d(j,i) = 0.0d0
         end do
      end do
c
c     zero out the virial and the intermolecular energy
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = 0.0d0
         end do
      end do
      einter = 0.0d0
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     call the local geometry energy and gradient routines
c
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_tortor)  call etortor1
c
c     call the van der Waals energy and gradient routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
      end if
c
c     call the electrostatic energy and gradient routines
c
      if (use_charge)  call echarge1
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
      if (use_mpole .or. use_polar)  call empole1
      if (use_rxnfld)  call erxnfld1
c
c     call any miscellaneous energy and gradient routines
c
      if (use_solv)  call esolv1
      if (use_metal)  call emetal1
      if (use_geom)  call egeom1
      if (use_extra)  call extra1
c
c     call 3d-rism
c
      fu = freeunit ()
      if (vdwtyp .eq. 'LENNARD-JONES') then
        open (fu, file='xyz.dat')
        sigfact = 2.0d0 ** (5.0d0/6.0d0)
        write(fu, *) n
        do i=1, n
          sigu(i)=rad(class(i))*sigfact
          epsu(i)=eps(class(i))
          qu(i)=pchg(i)
          write(fu,'(f10.4,f10.6,f10.3,3f10.5)')
     .         qu(i), sigu(i), epsu(i), x(i), y(i), z(i)
        end do
        close(fu)
      else
        stop
      end if
      if (count == 0) then
        call system('cat 3drism.dat xyz.dat > 3drism.inp')
        if (gpu == -1) then
           call system('3drism 3drism.inp >> 3drism.out')
        else
           write(gn, '(i4)') gpu
           call system('3drism-cuda 3drism.inp '//gn//' >> 3drism.out')
        endif
      end if
      open (fu, file='3drism.gra',status='old')
      do i = 1, n
        read (fu, *) de3d(1,i), de3d(2,i), de3d(3,i)
      end do
      close(fu)
      count = mod(count + 1, num)
      open (fu, file='3drism.xmu',status='old')
      read (fu, *)  dummy, dummy, dummy, dummy, dummy
      close(fu)
      read(dummy,*) e3d
c
c     sum up to get the total energy and first derivatives
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
     &          + ep + er + es + elf + eg + ex
      energy = esum
      sfe = e3d
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
     &                      + det(j,i) + dept(j,i) + debt(j,i)
     &                      + dett(j,i) + dev(j,i) + dec(j,i)
     &                      + decd(j,i) + ded(j,i) + dem(j,i)
     &                      + dep(j,i) + der(j,i) + des(j,i)
     &                      + delf(j,i) + deg(j,i) + dex(j,i)
     &                      + de3d(j,i)
            derivs(j,i) = desum(j,i)
         end do
      end do
c
c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end
