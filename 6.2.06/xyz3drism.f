c
c     "xyz3drism" takes as input a Cartesian coordinates file,
c     then converts to and writes out a file in 3D-RISM format
c
c
      program xyz3drism
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'vdwpot.i'
      include 'charge.i'
      include 'kvdws.i'
      include 'atmtyp.i'
      integer i,irism,freeunit
      real*8 sigfact,sigu(maxatm),epsu(maxatm),qu(maxatm)
      character*60 rismfile
c
c
c     get the Cartesian coordinates file for the system
c
      call initial
      call getxyz
c
c     get atomic number of each atom and count the molecules
c
c      call field
c      call katom
c      call molecule
c
c     get potential parameter (LJ and charge) of each atom
c
      call mechanic
c
c
c
      if (vdwtyp.eq.'LENNARD-JONES') then
        sigfact=2.0d0**(5.0d0/6.0d0)
        do i=1,n
          sigu(i)=rad(class(i))*sigfact
          epsu(i)=eps(class(i))
          qu(i)=pchg(i)
        end do
      else
        stop
      end if
c
c     write out the 3D-RISM input file
c
      rismfile = filename(1:leng)//'.3d'
      irism = freeunit ()
      open (unit=irism,file=rismfile,status='new')
      write(irism,*) n
      do i=1,n
        write(irism,100) qu(i),sigu(i),epsu(i),x(i),y(i),z(i)
      end do
  100 format(F10.4,F10.6,F10.3,3F10.5)
      close (unit=irism)
c
c     perform any final tasks before program exit
c
      call final
      end
