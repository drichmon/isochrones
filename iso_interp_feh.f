!------------------------------------------------------------------------!
      program iso_interp_feh
!     Aaron Dotter                                                       !
!     Contact aaron.dotter@gmail.com with questions or comments          !
!     2008 orig: includes full range of photometric systems              !
!     2008 July: bug fix in setting [a/Fe] for Y=0.33,0.4 (AS)           !
!     2008 Sept: minor fixes, set default interpolation to cubic (JR)    !
!     2008 Oct:  bug fix in setting -0.5<[Fe/H]<0 with Y=0.33,0.4 (PG)   !
!     2008 Nov:  perform interpolation on [Fe/H] instead of Z/X gives    !
!                smooth results over the entire range of [Fe/H]   (TP)   !
!     2011 Jan:  added support for HST/WFC3 photometric system           !
!     2011 Feb:  fixed quiet mode bugs for Y=0.33,0.40 (MG)              !
!                also fixed bug in [Fe/H] + [alpha/Fe] choices           !
!     2011 Oct:  overhaul and update list of photometric systems         !
!                                                                        !
!     Program interpolates amongst existing isochrones with any number   !
!     of ages (<50) and writes the new isochrones in the same format.    !
!                                                                        !
!     To compile:      g77 -o iso_interp_feh iso_interp_feh.f            !
!                                                                        !
!     The program must be able to access the isochrone grids via the     !
!     following path: ./isochrones/[bvri/jc2mass/hst_acs/etc.]           !
!     Interpolation is performed on [Fe/H] at constant Y and [a/Fe].     !
!     The order or degree of interpolation is set by the integer NISO.   !
!     For example, NISO=2 gives linear interpolation while NISO=4 gives  !
!     cubic interpolation. Default is cubic but reverts to linear near   !
!     the edges of the [Fe/H] grid.                                      !
!                                                                        !
!     Output from this program can be split into individual ages with    !
!     the isolf_split program.                                           !
!                                                                        !
!     This program may be run quietly from the command line              !
!     (if 'IArgC' and 'getarg' are supported) by uncommenting the lines  !
!     under 'run in quiet mode' below. To run from the command line      !
!     simply supply the arguments in the same order as requested by      !
!     the interactive mode, for example:                                 !
!                                                                        !
!             ./iso_interp_feh 1 1 2 -1.3 test.iso                       !
!                                                                        
!     will generate a new file called 'test.iso' using UBVRIJHKsKp       !
!     colors with Y=0.25+1.5Z, [alpha/Fe]=0, [Fe/H]=-1.3.                !
!------------------------------------------------------------------------!

      implicit none

      integer, parameter :: niso=4, max=50, nl=350, nz1=8, nz2=6
      integer, parameter :: isonum=1, nafe=6, ncol=15, maxcol=77
      integer :: i,j,k,l,ifeh,ieep,neweep,fn
      integer :: newimin,newimax,nz,iafe,nint,iclr,iy
      integer :: neeps(niso,max),eeps(niso,max,nl),nages(niso),
     >     ioff(niso),cols,ierr
      character(len= 5) :: astr
      character(len= 6) :: estr, magstr
      character(len=10) :: arg
      character(len=19) :: nagestr
      character(len=50) :: phot(niso)
      character(len=53) :: hdr1
      character(len=60) :: sep
      character(len=128) :: iso_dir, outfile
      character(len=1024) :: hdr4
      double precision, parameter :: ln10 = 2.3025850929940459d0
      double precision :: age(niso,max),mass(niso,max,nl),afe0(nafe),
     >     teff(niso,max,nl),gl(niso,max,nl),lum(niso,max,nl),mixl,
     >     color(niso,max,nl,maxcol),newfeh,qf(niso),ff(niso),newm,newt,
     >     newg,newl,newclr(maxcol),y(niso),z(niso),zeff(niso),newxeff,
     >     afe(niso),newzeff,newz,newy,newfeh0,feh0(nz1),feh(niso)
      logical :: lquiet = .false.

      iso_dir = 'isochrones/'
      !iso_dir = 'isochrones/'
      feh0=(/-2.5d0,-2.0d0,-1.5d0,-1.0d0,-0.5d0,0.0d0,0.3d0,0.5d0/)
      afe0=(/-0.2d0,0.0d0,0.2d0,0.4d0,0.6d0,0.8d0/)
      sep = '#----------------------------------------------------'
      
! run in quiet mode
      if(IArgC()>0) then
         lquiet=.true.
         call getarg(1,arg)
         read(arg,'(i2)') iclr
         call getarg(2,arg)
         read(arg,'(i1)') iy
         call getarg(3,arg)
         read(arg,'(i1)') iafe
         
         if(iy==1) then
            nz=nz1
         else
            nz=nz2
            ! make sure that only two options for [alpha/Fe]
            ! are allowed if Y=0.33 or 0.40
            if(iafe/=4.and.iafe/=2) then
               write(*,*) 'INVALID CHOICE OF [a/Fe]'
               write(*,*) 'For Y=0.33, 0.40, options are:'
               write(*,*) ' [a/Fe]=  [2]  0.0'
               write(*,*) '          [4] +0.4'
               stop
            endif
         endif

         call getarg(4,arg)
         read(arg,*) newfeh
         newfeh=newfeh*1d0
         call getarg(5,outfile)

      else
! run interactively
! obtain user input:
      write(*,*) '-----------------------------------------------------'
      write(*,*) '    Isochrone [Fe/H] Interpolation Program'
      write(*,*) '-----------------------------------------------------'
      write(*,*) ' Enter choice of color-Teff transformation: [1-14]'
      write(*,*) '     [1] UBVRI+2MASS+Kepler: Synthetic (Vega)    '
      write(*,*) '     [2] Washingon+DDO51   : Synthetic (Vega)    '
      write(*,*) '     [3] HST/WFCP2         : Synthetic (Vega)    '
      write(*,*) '     [4] HST/ACS-WFC       : Synthetic (Vega)    '
      write(*,*) '     [5] HST/ACS-HRC       : Synthetic (Vega)    '
      write(*,*) '     [6] HST/WFC3-UVIS+IR  : Synthetic (Vega)    '
      write(*,*) '     [7] Spitzer IRAC      : Synthetic (Vega)    '
      write(*,*) '     [8] UKIDSS            : Synthetic (Vega)    '
      write(*,*) '     [9] WISE              : Synthetic (Vega)    '
      write(*,*) '    [10] CFHT/MegaCam      : Synthetic (AB)      '
      write(*,*) '    [11] SDSS/ugriz        : Synthetic (AB)      '
      write(*,*) '    [12] PanSTARRS         : Synthetic (AB)      '
      write(*,*) '    [13] SkyMapper         : Synthetic (AB)      '
      write(*,*) '    [14] DECam             : Synthetic (AB)      '
      write(*,*) '    [15] BV(RI)c+uvby      : VandenBerg+Clem et al.'

      read  *, iclr
      if(iclr<1.or.iclr>ncol) stop 'INVALID CHOICE OF COLOR'

      write(*,*) '-----------------------------------------------'
      write(*,*) ' Enter choice of Y: [1-3]'
      write(*,*) ' Y=  [1] 0.245 + 1.5 Z'
      write(*,*) '     [2] 0.33           ([Fe/H]<=0 only)'
      write(*,*) '     [3] 0.40           ([Fe/H]<=0 only)'
      read  *, iy
      if(iy<1.or.iy>3) stop 'INVALID CHOICE OF Y'

      if(iy==1) then
         write(*,*) '-----------------------------------------------'
         write(*,*) ' Enter choice of [a/Fe]: [1-6]'
         write(*,*) ' [a/Fe]=  [1] -0.2'
         write(*,*) '          [2]  0.0'
         write(*,*) '          [3] +0.2'
         write(*,*) '          [4] +0.4   ([Fe/H]<=0 only)'
         write(*,*) '          [5] +0.6   ([Fe/H]<=0 only)'
         write(*,*) '          [6] +0.8   ([Fe/H]<=0 only)'
         read  *, iafe
         if(iafe<1.or.iafe>6) stop 'INVALID CHOICE OF [a/Fe]'
      else
         write(*,*) '-----------------------------------------------'
         write(*,*) ' Enter choice of [a/Fe]: [2,4]'
         write(*,*) ' [a/Fe]=  [2]  0.0'
         write(*,*) '          [4] +0.4'
         read  *, iafe
         if(iafe/=2 .and. iafe/=4) stop 'INVALID CHOICE OF [a/Fe]'
      endif
      
      if(iy==1) then
         nz=nz1
      else
         nz=nz2
      endif

      write(*,*) '-----------------------------------------------'
      write(*,*) ' Available [Fe/H] range is ',feh0(1),'->',feh0(nz)
      write(*,*) ' Enter new [Fe/H] to interpolate:'
      read  *, newfeh

      write(*,*) '-----------------------------------------------'
      write(*,*) ' Enter new filename:'
      read(*,'(a)') outfile
      write(*,*) '-----------------------------------------------'

      endif !if(IargC()...

      !check input values
      if(newfeh<feh0(1).or.newfeh>feh0(nz)) then
         write(*,*) 'EXTRAPOLATION WARNING - [Fe/H] out of bounds'
      endif
      if(newfeh>0.0d0.and.iafe>3) then
         write(*,*) 'INVALID COMBINATION OF [Fe/H] and [alpha/Fe]'
         stop
      endif

!   locate new [Fe/H] in grid
      call locate(feh0,nz,newfeh,ifeh)
!   nint is the number of points to be used for interpolation
      nint=niso                 
!   revert to linear interpolation if near the edge of the grid
      if(ifeh<nint/2.or.ifeh>nz-(nint/2)) nint=2
!   [a/Fe] > +0.2 only for [Fe/H]<=0
      if(iafe>3.and.newfeh>-0.5d0) nint=2

! this section of code reads in the isochrone files
      do j=1,nint                         
         fn=j+ifeh-nint/2
! read header
         open(isonum,file=trim(iso_file(fn,iy,iafe,iclr)),status='old')
         read(isonum,'(a16,i2,a6,i2)') nagestr,nages(j),magstr,cols
         read(isonum,*)
         read(isonum,'(a)') hdr1
         read(isonum,'(1x,f7.4,f8.4,1p2e11.4,0p2f7.2)') 
     *        mixl,y(j),z(j),zeff(j),feh(j),afe(j)
         read(isonum,*)
         read(isonum,'(25x,a65)') phot(j) 
         read(isonum,*)
         if(.not.lquiet) write(*,*) "READ IN [Fe/H]=",feh(j)
         
!   read lines of data
         do i=1,nages(j)
            read(isonum,'(a5,f6.3,a6,i3)') astr,age(j,i),estr,neeps(j,i)
            read(isonum,'(a)') hdr4
            do k=1,neeps(j,i)   !Mass,Teff,log g,Lum,colors of each EEP
               read(isonum,*) eeps(j,i,k),mass(j,i,k),teff(j,i,k),
     *              gl(j,i,k),lum(j,i,k),color(j,i,k,1:cols)
            enddo
            if(i<nages(j)) read(isonum,*)
            if(i<nages(j)) read(isonum,*)
         enddo
         close(isonum)
      enddo                     
! finished reading in isochrones

! make sure all files contain the same number and range ages:
      do j=2,nint
         if(nages(j)/=nages(j-1)) stop 'INCORRECT NUMBER OF AGES!'
         if(phot(j)/=phot(j-1))stop 'COLOR-TEFF NOT SAME IN ALL FILES'
         do i=1,nages(j)
            if(age(j,i)/=age(j-1,i)) stop 'AGES DO NOT MATCH!'
         enddo
      enddo

! interpolation in [Fe/H] gives smoother results
      do j=1,nint
         qf(j)=feh(j)
      enddo
      call interp(qf,ff,newfeh,nint)
      
! set new X, Y, Z, and calculate [Fe/H] from them
! if result differs from desired result by >5% print warning message
      newy=0.0d0
      newz=0.0d0
      do j=1,nint
         newy=newy+ff(j)*y(j)
         newz=newz+ff(j)*z(j)
      enddo

! adjust X and Z if [a/Fe] non-zero by introducing Xeff and Zeff
! calculate Zeff using relation of Salaris, Chieffi & Straniero
! (1993, ApJ, 414, 580)
      newzeff=newz/(6.38d-1*exp(ln10*afe0(iafe)) + 3.62d-1)
      newxeff=1.0d0-newy-newzeff
      newfeh0=log10(newzeff/(newxeff*2.29d-2))
      if(.not.lquiet .and.
     *     abs(newfeh-newfeh0)/abs(newfeh)>5.0d-2) then
         write(*,*) '***WARNING - POOR X,Z INTERPOLATION***'
         write(*,*) '[Fe/H](Input) - [Fe/H](Interp)=',newfeh-newfeh0
      endif

! open new file and write header information
      open(isonum,file=outfile,iostat=ierr)
      if(ierr/=0) write(*,*) outfile
      write(isonum,'(a16,i2,a6,i2)') nagestr,nages(1),magstr,cols
      write(isonum,'(a)') sep
      write(isonum,'(a)') hdr1
      write(isonum,'(a1,f7.4,f8.4,1p2e11.4,0p2f7.2)') 
     *     '#',mixl,newy,newz,newzeff,newfeh,afe0(iafe)
      write(isonum,'(a)') sep
      write(isonum,'(a25,a)') "#**PHOTOMETRIC SYSTEM**: ",trim(phot(1))
      write(isonum,'(a)') sep
      
! loop through ages 
      do j=1,nages(1)
         newimax=eeps(1,j,neeps(1,j))
         newimin=eeps(1,j,1)
         do k=2,nint
            if(newimax>eeps(k,j,neeps(k,j)))
     *          newimax=eeps(k,j,neeps(k,j))
            if(newimin<eeps(k,j,1)) newimin=eeps(k,j,1)
         enddo
         neweep=newimax-newimin+1
         do k=1,nint
            ioff(k)=newimin-eeps(k,j,1) 
         enddo

         write(isonum,'(a5,f6.3,a6,i3)') astr,age(1,j),estr,neweep
         write(isonum,'(a)') trim(hdr4)

! interpolate and write results for each EEP
         do i=1,neweep
            do k=1,nint
               do l=k,nint
!   make sure EEPs are constant across different files
                  if(eeps(k,j,i+ioff(k))/=eeps(l,j,i+ioff(l))) then
                     write(*,*) "EEPS NOT EQUAL"
                  endif
               enddo
            enddo
!   zero new quantities
            newm=0.0d0
            newt=0.0d0
            newg=0.0d0
            newl=0.0d0
            newclr=0.0d0
!   perform interpolation            
            do k=1,nint
               newm=newm+ff(k)*mass(k,j,i+ioff(k))
               newt=newt+ff(k)*teff(k,j,i+ioff(k))
               newg=newg+ff(k)*  gl(k,j,i+ioff(k))
               newl=newl+ff(k)* lum(k,j,i+ioff(k))
               do l=1,cols
                  newclr(l)=newclr(l)+ff(k)*color(k,j,i+ioff(k),l)
               enddo
            enddo
! write new isochrone(s) to file
            ieep=i-1+newimin
            write(isonum,'(i4,f10.6,99f8.4)')ieep,newm,newt,newg,newl,
     *           newclr(1:cols)
         enddo
         if(j<nages(1)) write(isonum,*)
         if(j<nages(1)) write(isonum,*)
      enddo
      close(isonum)
      
! all done!
      contains
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! finds the index, j, of x in array xx(n). assuming xx is sorted, x lies
! between xx(j) and xx(j+1) -- from NUMERICAL RECIPES IN FORTRAN 77  
      subroutine locate(xx,n,x,j)
      integer, intent(in) :: n
      double precision, intent(in) :: xx(n), x
      integer, intent(out) :: j
      integer :: jl, ju, jm
      jl=0
      ju=n+1
      do while(ju-jl>1)
         jm=(ju+jl)/2
         if((xx(n)>xx(1)).eqv.(x>xx(jm)))then
            jl=jm
         else
            ju=jm
         endif
      enddo
      j=jl
      end subroutine locate
!------------------------------------------------------------------------
      subroutine interp(a,b,x,n)
! {a} are the tabulated values for use in interpolation
! {b} are coefficients of the interpolating polynomial
!  x  is the abscissa to be interpolated
!  n  is the number of points to be used, interpolating polynomial
!     has order n-1 
      integer, intent(in) :: n
      double precision, intent(in) :: a(n), x
      double precision, intent(out) :: b(n)
      integer :: i,j
      do i=1,n
         b(i)=1.0d0
         do j=1,n
            if(j/=i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
         enddo
      enddo
      end subroutine interp
!------------------------------------------------------------------------
      character(len=256) function iso_file(fn,iy,iafe,iclr)
! there are currently 30 isochrone files, this function returns
! the name of the proper file based on [Fe/H] and [a/Fe]
      integer, intent(in) :: fn, iafe, iclr, iy
      character(len=6) :: FeH(nz1)
      character(len=5) :: afe(nafe)
      character(len=11) :: suf(ncol)
      feh = (/'fehm25','fehm20','fehm15','fehm10','fehm05','fehp00',
     >     'fehp03','fehp05'/)
      afe = (/'afem2','afep0','afep2','afep4','afep6','afep8'/)
      
      suf = (/'UBVRIJHKsKp','WashDDOuvby','HST_WFPC2  ','HST_ACSWF  ',
     >        'HST_ACSHR  ','HST_WFC3   ','SPITZER    ','UKIDSS     ',
     >        'WISE       ','CFHTugriz  ','SDSSugriz  ','PanSTARRS  ',
     >        'SkyMapper  ','DECam      ','BVRIuvby   '             /)

      iso_file = ''
      iso_file = trim(iso_dir) // trim(suf(iclr)) // '/' // FeH(fn) // 
     >     afe(iafe)
      if(iy==2) iso_file = trim(iso_file) // 'y33'
      if(iy==3) iso_file = trim(iso_file) // 'y40'
      iso_file = trim(iso_file) // '.' // trim(suf(iclr))
      end function iso_file
!------------------------------------------------------------------------
      end program iso_interp_feh
