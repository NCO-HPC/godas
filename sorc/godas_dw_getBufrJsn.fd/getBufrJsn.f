!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  GETBUFRJSN
!   PRGMMR: Diane Stokes   ORG: NP23        DATE: 2007-02-15
!
! ABSTRACT:  Extract Jason data from BUFR files.
!
! PROGRAM HISTORY LOG:
! unknown     Vera Gerald (original alps.f code)
! 2006-08-08  David W. Behringer - code clean-up
! 2006-11-30  Diane C. Stokes    - sort data by time (needed for later steps), 
!                                  add sanity checks on data
! 2008-02-15  Diane C. Stokes    - upgrade code to process multiple input files
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - LIST OF DAILY ALTIMETER BUFR FILES
!     UNIT 21  - ALTIMETER DATA IN BUFR. SEE REMARKS.
!
!   OUTPUT FILES:
!     UNIT 51  - DAILY ALTIMETRY FILE IN ASCII.  SEE REMARKS.
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  none
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - qsort
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!       BUFRLIB  - ufbtab, closbf
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - qsort
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!       BUFRLIB  - openbf, ufbget, ufbtab, closbf
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
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
      program getBufrJsn
!
      implicit none

      character(len=80) :: str, dayFile
      character(len=80),parameter :: tstr='YEAR MNTH DAYS HOUR MINU &
           & CLATH CLONH SLHD1 ORBN OBQL SACYLN'
      character(len=8) :: cdate

      real(8),allocatable :: TAB_8(:,:)
      real(8) :: bmiss   ! will define bmiss with first call to ufbtab

      integer(8),allocatable :: iwork(:)
      integer,allocatable :: iord(:)
      integer, parameter :: lunin=12, mxts=11, imn=1, max_prt_ibad_orb=50
      integer :: mxtb, ntab, istat, ierr
      integer :: kdate,iymd
! various counters
      integer :: ird=0, itx=0
      integer :: imiss=0, ibad_slhd=0, ibad_date=0, ibad_loc=0, ibad_orb=0, inot_fnl=0
      integer :: iob,irec

      real :: xtab, xlog
      integer :: keysize

      real :: slhd, xlat, xlon 
      integer :: iobql,iyr,mon,iday,ihr,min,iorb,isacyc


      call w3tagb('GODAS2_GETBUFRJSN',2008,0045,0045,'NP23')

!
!  read a daily file name from unit 11 and open the daily file. process
!  each file until the list is exhausted
!
      do while (.TRUE.)
        read(11,'(i8,1x,a)',end=200,err=900) kdate,dayFile
        open(lunin,file=trim(dayFile),form='unformatted',status='old',err=910)
        call openbf(lunin,'IN',lunin)
        write(cdate,'(i8.8)')kdate
        open(51,file=cdate//'.sshd.qc2048',form='formatted',err=915)

! Get data count in bufr file and allocate arrays 
!  (negative unit number in call to ufbtab returns ob count only)
!  At same time, getting BMISS value returned from ufbtab.
!  The actual missing value returned from bufrlib routines may depend 
!  on the machine the lib was compiled on (roundoff issues).  
!  So, get it on the fly.  Calling ufbtab with negative unit number returns 
!  missing in array passed as 2nd argument.
       call ufbtab (-lunin,bmiss,1,1,ntab,' ')
       mxtb=ntab
       print*,'BMISS is:',bmiss

       allocate(tab_8(mxts,mxtb),stat=istat);if(istat.ne.0) goto 930
       allocate(iwork(mxtb)     ,stat=istat);if(istat.ne.0) goto 930
       allocate(iord(mxtb)      ,stat=istat);if(istat.ne.0) goto 930

! store all bufr obs in table tab_8
       call ufbtab (lunin,tab_8,mxts,mxtb,ntab,tstr)
       call closbf(lunin)
       ird=ntab
 
! sort obs by date/time.  order is kept in array iord
 
!   include current record number in key for uniqueness
!   if ntab gt 9999999, we need to modify code because iwork will not fit in int(8)
       if(ntab.gt.9999999)then
         write(6,*) 'Need to work on sort method.  Current key won''t fit in int(8)'
         write(6,*) 'See if "orders" available.'
!        call errexit(12)
         call exit(12)
       endif
       xtab=ntab
       xlog=log10(xtab)
       keysize=xlog+1
      
      itx=0
      imiss=0
      ibad_slhd=0
      ibad_date=0
      ibad_loc=0
      ibad_orb=0
      inot_fnl=0

      do iob=1,ntab
        iord(iob)=iob
        iwork(iob)=tab_8(1,iob)*10E7+tab_8(2,iob)*10E5+tab_8(3,iob)*10E3 &
      &   + tab_8(4,iob)*100+tab_8(5,iob)
        iwork(iob)=iwork(iob)*(10**keysize)+iob
      enddo
      call qSort(iwork,iord,imn,ntab)
! don't need iwork anymore
      deallocate(iwork,stat=istat);if(istat.ne.0) goto 940

      do iob=1,ntab
        irec=iord(iob)

! skip records that are missing any of the required fields
        if(any(tab_8(:,irec).ge.bmiss))then
          imiss=imiss+1
          cycle
        endif

! check that slhd within reasonable range.
        slhd=tab_8(8,irec)
        if (slhd .le. -10.00 .or. slhd .ge. 10.00) then  
          ibad_slhd = ibad_slhd + 1
          cycle
        endif
 
! obql = 2048 : final product.  keep these only.
!
        iobql = nint(tab_8(10,irec))
        if (iobql .eq. 2048) then

          iyr = nint(tab_8(1,irec))
          mon = nint(tab_8(2,irec))
          iday = nint(tab_8(3,irec))
          ihr = nint(tab_8(4,irec))
          min = nint(tab_8(5,irec))
          xlat = tab_8(6,irec)
          xlon = tab_8(7,irec)
          if(xlat.lt.-90.or.xlat.gt.90.or.xlon.lt.-180.or.xlon.gt.180)then
            write(6,'(a,4x,2f15.2)') 'Skip ob with bad lat/lon:',xlat,xlon
            ibad_loc=ibad_loc+1
            cycle
          endif
          iorb = nint(tab_8(9,irec))
          isacyc = nint(tab_8(11,irec))
          if(iorb.gt.99999.or.isacyc.gt.99999)then
            if(ibad_orb.lt.max_prt_ibad_orb)then
              write(6,'(a,2x,2i12)') 'Skip ob with bad orbit or cycle number:',iorb, isacyc
            else
              if(ibad_orb.eq.max_prt_ibad_orb) then
                write(6,*) ''
                write(6,*) '** NO ADDITIONAL BAD ORBIT WARNINGS PRINTED **'
                write(6,*) '**    CHECK TOTAL COUNT IN SUMMARY TABLE    **'
              endif
            endif
            ibad_orb=ibad_orb+1
            cycle
          endif
!
          iymd = iyr*10000 + mon*100 + iday
          if(iymd.ne.kdate)then
            print*,'Skip report with wrong date:',iymd
            ibad_date=ibad_date+1
            cycle
          endif
          write (51,71)iymd,ihr,min,isacyc,iorb,xlat,xlon,slhd
          itx = itx + 1
        else
         print*,'Skip ob with iobql=',iobql
         inot_fnl=inot_fnl+1
        endif
!
      enddo
71    format(i9,2i3,i5,1x,i4,2f12.5,1x,f7.3)

! deallocate arrays
      deallocate(tab_8,stat=istat);if(istat.ne.0) goto 940
      deallocate(iord ,stat=istat);if(istat.ne.0) goto 940

      close(51)
!
      write(6,'(/)') 
      write(6,'(a,i12)') 'For date: ',kdate
      write(6,'(a,i12)') 'Records read:                    ', ird
      write(6,'(a,i12)') 'Records missing required field:  ', imiss
      write(6,'(a,i12)') 'Records with BAD SLHD:           ', ibad_slhd
      write(6,'(a,i12)') 'Records with BAD date:           ', ibad_date
      write(6,'(a,i12)') 'Records with BAD lat or lon:     ', ibad_loc
      write(6,'(a,i12)') 'Records with BAD orbit or cycle: ', ibad_orb
      write(6,'(a,i12)') 'Records extracted:               ', itx
      write(6,'(//)') 

      if(ibad_date.ne.0)then
        write(6,'(a)') ' ************************************************' 
        write(6,*)     ' WARNING:  ',ibad_date, 'WRONG DATES IN BUFR FILE ',dayFile 
        write(6,'(a)') ' ************************************************' 
      endif

      enddo

  200 continue
!
      call w3tage('GODAS2_GETBUFRJSN')
      call errexit(0)
!     call exit(0)

  900 continue
      ierr=90
      write(6,*)  'ERROR OPENING FILE LIST.  EXIT ',ierr
      go to 990
  910 continue
      ierr=91
      write(6,*)  'ERROR OPENING BUFR FILE.  EXIT ',ierr
      go to 990
  915 continue
      ierr=91
      write(6,*)  'ERROR OPENING ASCII OUTPUT FILE.  EXIT ',ierr
      go to 990
  920 continue
      ierr=92
      write(6,*) 'ERROR FINDING BUFR VALUE FOR MISSING DATA.  EXIT '
      go to 990
  930 continue
      ierr=93
      write(6,*) 'ERROR ',istat,' ALLOCATING ARRAYS.  EXIT ',ierr
      go to 990
  940 continue
      ierr=94
      write(6,*) 'ERROR ',istat,' DEALLOCATING ARRAYS.  EXIT ',ierr
      go to 990
  990 continue
      call w3tage('GODAS2_GETBUFRJSN')
      call errexit(ierr)
!     call exit(ierr)

!

      contains
!
! -------------------------------------------------------------- 
!
      recursive subroutine qSort(a, b, lo0, hi0)

!
      integer(8), intent(inout) ::  a(:)
      integer(8) :: mid
      integer, intent(inout) :: b(:)
      integer, intent(in) :: lo0, hi0
      integer :: lo, hi
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
      integer, intent(in) :: i, j
      integer(8), intent(inout) :: a(:)
      integer(8) :: Ta
      integer, intent(inout) :: b(:)
      integer :: Tb
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
      end program getBufrJsn
