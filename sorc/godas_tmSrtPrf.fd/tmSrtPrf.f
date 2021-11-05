!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  TMSRTPRF
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2003-07-24
!
! ABSTRACT:  Sort by date and time the profiles in a file of
!            subsurface temperature data
!
! PROGRAM HISTORY LOG:
! 2003-07-24  David W. Behringer
! 2004-01-30  David W. Behringer - cosmetic change in declaring and setting
!                                  the array ndxd
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - TEMPERATURE PROFILE DATA IN IEEE
!
!   OUTPUT FILES:
!     UNIT 51  - SORTED TEMPERATURE PROFILE DATA IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     UNIT 80  - SCRATCH FILE FOR PROFILE DATA
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - qSort, swap
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - qSort, swap
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     None
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program tmSrtPrf
!
!  tmSrtPrf sorts profiles by time
!
      character csign*8,sid*2,dtyp*2,qkey*1
      character str*80
      real, allocatable, dimension(:) :: pt, pz
      integer(kind=8), allocatable, dimension(:) :: ndxd
      integer(kind=8), parameter :: iyfct = 100000000
      integer, allocatable, dimension(:) :: nprfl
!
      call w3tagb('GODAS_TMSRTPRF',2003,0164,0164,'NP23')
!
! open profile file 
!
      open (11, form='unformatted', status='old', &
                            & access='sequential', err=110)
!
      nprf = 0
      npmx = 0
      do while (.true.)
        read (11, end=100, err=120) iyear,idate,csign,sid, &
                                       &  dtyp,qkey,yp,xp,np
        nprf = nprf + 1
        if (np .gt. npmx) npmx = np
      end do
  100 continue
!
      rewind 11
!
      allocate(pt(npmx))
      allocate(pz(npmx))
      allocate(ndxd(nprf))
      allocate(nprfl(nprf))
!
! open direct access scratch file 
!
      nb = 4*(2*npmx + 5) + 13
      open (80, status='scratch', form='unformatted', &
                                     & access='direct',recl=nb)
!
! begin loop on profile file
!
      do n=1,nprf
        read (11, err=120) iyear,idate,csign,sid, &
                     &  dtyp,qkey,yp,xp,np,(pz(k),pt(k),k=1,np)
        ndxd(n) = iyear*iyfct + idate
        nprfl(n) = n
        write (80,rec=n) iyear,idate,csign,sid,dtyp,qkey,yp,xp, &
                           &  np,pz,pt
      end do
!
      close (11)
!
! sort on ndxd
!
      imn = 1
      imx = nprf
      call qSort(ndxd,nprfl,imn,imx)
!
! write a new profile file in time sequence
!
      open (51, form='unformatted', access='sequential', err=510)
!
! begin loop on list file
!
      do n=1,nprf
        read (80,rec=nprfl(n)) iyear,idate,csign,sid,dtyp,qkey, &
                           &  yp,xp,np,pz,pt
        write (51, err=520) iyear,idate,csign,sid,dtyp,qkey,yp,xp, &
                           &  np,(pz(k),pt(k),k=1,np)
      end do
!
      close (51)
      close (80)
!
      call w3tage('GODAS_TMSRTPRF')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a)') 'Error opening profile file on unit 11'
      call w3tage('GODAS_TMSRTPRF')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading profile file on unit 11'
      call w3tage('GODAS_TMSRTPRF')
      call errexit(12)
!     call exit(12)
!
  510 write(6,'(a)') 'Error opening profile file on unit 51'
      call w3tage('GODAS_TMSRTPRF')
      call errexit(51)
!     call exit(51)
!
  520 write(6,'(a)') 'Error writing profile file on unit 51'
      call w3tage('GODAS_TMSRTPRF')
      call errexit(52)
!     call exit(52)
!
      contains
!
! -------------------------------------------------------------- 
!
      recursive subroutine qSort(a, b, lo0, hi0)
!
      integer*8 a(*), mid
      integer b(*), lo0, hi0, lo, hi
!
      lo = lo0
      hi = hi0
!
      if (hi0 .gt. lo0) then
!
        mid = a( ( lo0 + hi0 ) / 2 )
!
        do while( lo .le. hi )
!
          do while( (lo .lt. hi0) .and. (a(lo) .lt. mid) )
            lo = lo + 1
          end do
!
          do while( (hi .gt. lo0) .and. (a(hi) .gt. mid) )
            hi = hi - 1
          end do
!
          if ( lo .le. hi ) then
            call swap(a, b, lo, hi)
            lo = lo + 1
            hi = hi - 1
          end if
!
        end do
!
        if (lo0 .lt. hi) then
          call qSort(a, b, lo0, hi)
        end if
!
        if (lo .lt. hi0) then
          call qSort(a, b, lo, hi0)
        end if
!
      end if
!
      end subroutine qSort
!
! -------------------------------------------------------------------
!
      subroutine swap(a, b, i, j)
!
      integer*8 a(*), Ta
      integer b(*), i, j, Tb
!
      Ta = a(i)
      a(i) = a(j)
      a(j) = Ta
!
      Tb = b(i)
      b(i) = b(j)
      b(j) = Tb
!
      end subroutine swap
!
      end program tmSrtPrf
