!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  MKASMALT
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2006-08-08
!
! ABSTRACT:  Write altimetry observations in form used by GODAS
!
! PROGRAM HISTORY LOG:
! 2006-08-06  David W. Behringer
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - LIST OF ALTIMETRY CYCLE FILES IN ASCII
!     UNIT 21  - GRID/MASK FOR OCEAN MODEL IN IEEE
!     UNIT 31  - ALTIMETRY CYCLE DATA IN ASCII (SEE REMARKS)
!
!   OUTPUT FILES:
!     UNIT 51  - ALTIMETRY DATA FOR ASSIMILATION IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     UNIT 81  - SCRATCH FILE FOR ALTIMETRY DATA
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - real_indx, inbnds, qSort, swap
!                  jlnDayR, jlnDayC, calDayR, calDayC, weekDay
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - real_indx, inbnds, qSort,
!                  jlnDayR, calDayR, weekDay
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  21 - ERROR OPENING UNIT 21
!     COND =  22 - ERROR READING UNIT 21
!     COND =  31 - ERROR OPENING UNIT 31
!     COND =  32 - ERROR READING UNIT 31
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!     COND =  61 - ERROR READING FROM COMMAND LINE
!     COND =  81 - ERROR OPENING UNIT 81
!     COND =  82 - ERROR WRITING UNIT 81
!     COND =  83 - ERROR READING UNIT 81
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     Reads a date YYYYMMDD from command line
!     The names of the cycle files (UNIT 31) are read from a list on UNIT 11,
!      they are not specified externally.
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program mkAsmAlt
!
!  mkAsmAlt prepares SSH for use in assimilation
!
      use tmutil_mod
!
!  if std. error = 2cm, then EVR = 0.25
!                  3cm             0.1111
!                  4cm             0.0625
!                  5cm             0.04
!                 10cm             0.01
!
      real, save :: evr0 = 0.0625
!
!  latitude bounds are set at 30S and 40N. the std error is increased as
!  these limits are approached.  the effect is to smoothly reduce the 
!  influence of the ssh data approaching these bounds.
!
      real, save :: bndw = 5.0
!
      integer, parameter :: max_cyc_files=7
      integer :: imx, jmx, kmx, nwks, npts, ngpts, np
      real, allocatable, dimension(:) :: ar, xi, yj, xt, yt, zt, zw, tr, evr
      integer, allocatable, dimension(:,:) :: mskt
      integer(8), allocatable, dimension(:) :: ndxp
      integer, allocatable, dimension(:) :: ndxr
      integer :: year, month, day, yday
      integer :: year1, month1, day1, yday1
      integer :: year2, month2, day2, yday2
      integer :: yrc1, mnc1, dyc1, ydc1
      integer :: yrc2, mnc2, dyc2, ydc2
      real :: t1, tw, t2, tc, hgt
      real :: ax, ay, xin, yjn
      integer ::  week0, week, dow, year0, yr85 = 1985
      integer ::  nfl, kmin = -1
      character(len=50) inFile(max_cyc_files), fl
      character(len=32) stamp
      integer :: nob, wid
      real :: Lts, Ltn, Ltsb, Ltnb
      logical :: ltbnd = .TRUE.
      character(len=8) :: jsnHt = 'JSNHT   '
      character(len=8) :: chDate, wkDate(7)
      character(len=81) :: str
      real :: rwtr, rwa, rwxi, rwyj, rwevr
      character(len=8) :: rwtxt
      integer :: npts_total=0
!
      call w3tagb('GODAS_MKASMALT',2006,0164,0164,'NP23')
!
!     kmin = minimum number of model levels for assimilation of altimetry
!     nwks = time window (weeks) for altimetry data
!     Lts, Ltn - southern and northern limits for assimilation of altimetry
!
      kmin = 30
      nwks = 5
      Lts = -30.0
      Ltn = 40.0
      Ltsb = Lts + bndw
      Ltnb = Ltn - bndw
!
! get date of the desired analysis from command line and store it in chDate
!
      narg = iargc()
      if (narg .eq. 0) go to 610
      call getarg(1,str)
      chDate = str(1:8)
!
! read model coordinates and mask
!
      open(unit=21,form='UNFORMATTED', &
     &                      access='sequential', err=210)
      read(21, err=220) imx,jmx,kmx
      allocate (xt(imx))
      allocate (yt(jmx))
      allocate (zt(kmx))
      allocate (zw(kmx))
      allocate (mskt(imx,jmx))
      read(21, err=220) xt, yt, zt, zw
      read(21, err=220) mskt
      close(21)
!
! compute dates for nwks weeks of data (dates are for Wednesday of each week)
! the date in chDate (yyyymmdd) is the date of the desired analysis
!
      read(chDate,'(i4,i2,i2)') year, month, day
      year0 = year - 1
      yday = jlnDayR (year, month, day, year0)
!
      dow = weekDay(year, month, day)
      if (dow .ne. 4) then
        yday = yday + 4 - dow
        call calDayR(year, month, day, year0, yday)
      endif
!
      nw2 = nwks / 2
      yday = yday - 7*nw2
!
      do n=1,nwks
        call calDayR(year, month, day, year0, yday)
        write(wkDate(n),'(i4,i2.2,i2.2)') year,  month, day
        yday = yday + 7
      enddo
!
!  open file for altimetry data in GODAS format
!
      open(unit=51,access='SEQUENTIAL',form='UNFORMATTED',err=510)
!
!  begin main loop for assmbling weekly data files
!
      do nw=1,nwks
        read(wkDate(nw),'(i4)') year
        month = 1
        day = 1
        dow = weekDay(year, month, day)
        read(wkDate(nw),'(i4,i2,i2)') year, month, day
        year0 = year
        yday = jlnDayR (year, month, day, year0)
        week0 = 1 + (yday + dow - 2) / 7
        dow = weekDay(year, month, day)
        if (dow .ne. 4) then
          yday = yday + 4 - dow
          if (yday .lt. 1) then
            year0 = year0 - 1
            yday = jlnDayR (year, month, day, year0) + 4 - dow
          endif
          call calDayR(year, month, day, year0, yday)
        endif
        year0 = yr85
        yday = jlnDayR (year, month, day, year0)
        yday1 = yday - 3
        call calDayR(year1, month1, day1, year0, yday1)
        yday2 = yday + 3
        call calDayR(year2, month2, day2, year0, yday2)
!
        t1 = float(yday1)
        t2 = float(yday2) + 1.0
        tw = 0.5*(t1 + t2)
        wid = day + 100*(month + 100*mod(year,100))
!
! print date information to stdout
!
        write(6,'(/,a,i2,3(a,i4,a,i2.2,a,i2.2))') 'Heights for week ', &
     &     week0,' (',year,'/',month,'/',day,') between dates: ',  &
     &     year1,'/',month1,'/',day1,' -> ',year2,'/',month2,'/',day2
!
! open file with table listing information on altimeter cycle files
! the table file is TEXT and each line has the start date, end date, 
!   number of records, and name of an altimeter cycle file
!
! the subsequent loop determines which cycle files have data within the
!  current week
!
        open(unit=11,access='SEQUENTIAL',form='FORMATTED',err=110)
        nfl = 0
        do while (.TRUE.)
          read(11,'(i4,1x,i2,1x,i2,2x,i4,1x,i2,1x,i2,i7,a)',end=100,err=120) &
     &    yrc1, mnc1, dyc1, yrc2, mnc2, dyc2, n, fl
          ydc1 = jlnDayR (yrc1, mnc1, dyc1, year0)
          ydc2 = jlnDayR (yrc2, mnc2, dyc2, year0)
          if ((ydc1 .ge. yday1 .and. ydc1 .le. yday2) .or.   &
     &              (ydc2 .ge. yday1 .and. ydc2 .le. yday2)) then
            if(nfl.lt.max_cyc_files)then
              nfl = nfl + 1
              inFile(nfl) = adjustl(fl)
            else
              write(6,'(a)') 'WARNING:  More cycle files than expected.  Increase max_cyc_files.'
              write(6,'(a,a,a)') 'WARNING:  File ',trim(adjustl(fl)),' skipped.'
            endif
!         else if (ydc1 > yday2) then   ! commented in case files out of order, as in cajb_c000.txt prob 
!           exit
          endif
        enddo
100     continue
        close(11)
!
! open a scratch file and use it to accumulate data for the current week
!  from the cycle files
! the cycle files have 1-degree average, along-track data in TEXT format
!
        if (nfl .gt. 0) then
          ngpts = 0
          open(unit=81,status='SCRATCH',form='UNFORMATTED',  &
     &                               access='DIRECT',recl=16,err=810)
          do n=1,nfl
            open(unit=31,file=trim(inFile(n)),access='SEQUENTIAL', &
     &                                           form='FORMATTED',err=310)
            do while (.TRUE.)
              read(31,'(f9.2,f8.2,f9.2,2f8.2)',end=200,err=320)  &
     &                                            ay, ax, dum, tc, hgt
              yday = tc
              if (yday .ge. yday1 .and. yday .le. yday2) then
                ay = ay + 90.0
                ngpts = ngpts + 1
                write(81,rec=ngpts,err=820) ay, ax, tc, hgt
              endif
            enddo
200         continue
            close(31)
          enddo
!
! create a sorting index based on latitude and longitude
!
          allocate (ndxp(ngpts))
          allocate (ndxr(ngpts))
          do n=1,ngpts
            read(81,rec=n,err=830) ay, ax, tc, hgt
            ndxp(n) = int(ay*100.0)
            ndxp(n) = ndxp(n)*100000 + int(ax*100.0)
            ndxr(n) = n
          enddo
!
! sort on ndxp
!
          ns = 1
          ne = ngpts
          call qSort(ndxp,ndxr,ns,ne)
!
! allocate space to decimal indexes, ssh value, relative time, and error
!
          allocate (xi(ngpts))
          allocate (yj(ngpts))
          allocate (ar(ngpts))
          allocate (tr(ngpts))
          allocate (evr(ngpts))
          np = 0
          nob = 0
!
! read the sorted data and compute the decimal indexes and the assigned
!  errors of the observations.  The observation errors are made to increase
!  approaching the boundaries Lts and Ltn in order to taper the influence
!  of the observations toward those limits.
!
          do ng=1,ngpts
            read(81,rec=ndxr(ng),err=830) ay, ax, tc, hgt
            ay = ay - 90.0
            xin = real_indx(ax,xt,imx)
            yjn = real_indx(ay,yt,jmx)
            if (inbnds(xin,yjn,ax,ay)) then
              np = np + 1
              xi(np) = xin
              yj(np) = yjn
              tr(np) = tc - tw
              ar(np) = hgt
              evr(np) = evr0
              if (ltbnd) then
                if ((Ltsb-ay) .gt. 0.0) then
                  evr(np) = evr0 * exp(2.0*(ay-Ltsb)/bndw)
                 else if ((ay-Ltnb) .gt. 0.0) then
                  evr(np) = evr0 * exp(2.0*(Ltnb-ay)/bndw)
                endif
              endif
            else
              nob = nob + 1
            endif
          enddo
300       continue
          npts = np
          npts_total = npts_total+npts
!
          close(81)
!
! write data for one week in the GODAS IEEE format: time stamp and number
! of data records followed by the data records
!
          write(stamp,'(a,i2,a,i2,a,i4,a)') 'm/d/y=', &
     &             month, '/', day, '/', year, ', h:m:s=12: 0: 0'
          write(51,err=520) stamp, npts
!
          do n=1,npts
            write(51,err=520) jsnHt,tr(n),ar(n),xi(n),yj(n),evr(n)
          enddo
!
          write(6,'(a,i7)') 'Points Kept:   ', npts
          write(6,'(a,i7)') 'Out of Bounds: ', nob
!
          deallocate (ndxp)
          deallocate (ndxr)
          deallocate (xi)
          deallocate (yj)
          deallocate (ar)
          deallocate (tr)
          deallocate (evr)
!
        else
!
          npts = 0
!
! if there is no data for a particular week, then write a record with 
! a time stamp and a zero data count
!
          write(stamp,'(a,i2,a,i2,a,i4,a)') 'm/d/y=', &
     &             month, '/', day, '/', year, ', h:m:s=12: 0: 0'
          write(51,err=520) stamp, npts
!
          write(6,'(a,i4,a,i2.2,a,i2.2,a,i4,a,i2.2,a,i2.2)') &
     &     'There is no data for ',  &
     &     year1,'/',month1,'/',day1,' -> ',year2,'/',month2,'/',day2
        endif
      enddo
!
      close(51)
      if(npts_total.eq.0) &
          write(6,'(/,a)') '*** WARNING: NO VALID ALTIMETER DATA WAS KEPT ***'
!
      call w3tage('GODAS_MKASMALT')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a)') 'Error opening ave altimetry list file on unit 11'
      call w3tage('GODAS_MKASMALT')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading ave altimetry list file on unit 11'
      call w3tage('GODAS_MKASMALT')
      call errexit(12)
!     call exit(12)
!
  210 write(6,'(a)') 'Error opening grid/mask file on unit 21'
      call w3tage('GODAS_MKASMALT')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a)') 'Error reading grid/mask file on unit 21'
      call w3tage('GODAS_MKASMALT')
      call errexit(22)
!     call exit(22)
!
  310 write(6,'(a)') 'Error opening altimetry cycle file on unit 31'
      call w3tage('GODAS_MKASMALT')
      call errexit(31)
!     call exit(31)
!
  320 write(6,'(a)') 'Error reading altimetry cycle file on unit 31'
      call w3tage('GODAS_MKASMALT')
      call errexit(32)
!     call exit(32)
!
  510 write(6,'(a)') 'Error opening GODAS altimetry file on unit 51'
      call w3tage('GODAS_MKASMALT')
      call errexit(51)
!     call exit(51)
!
  520 write(6,'(a)') 'Error writing GODAS altimetry file on unit 51'
      call w3tage('GODAS_MKASMALT')
      call errexit(52)
!     call exit(52)
!
  610 write(6,'(a)') 'Error reading date from command line'
      call w3tage('GODAS_MKASMALT')
      call errexit(61)
!     call exit(61)
!
  810 write(6,'(a)') 'Error opening scratch file on unit 81'
      call w3tage('GODAS_MKASMALT')
      call errexit(81)
!     call exit(81)
!
  820 write(6,'(a)') 'Error writing scratch file on unit 81'
      call w3tage('GODAS_MKASMALT')
      call errexit(82)
!     call exit(82)
!
  830 write(6,'(a)') 'Error reading scratch file on unit 81'
      call w3tage('GODAS_MKASMALT')
      call errexit(83)
!     call exit(83)
!
      contains
!
! --------------------------------------------------------------
!
      function real_indx (p, pt, nmx)
!
      real p, pt(*)
      integer nmx
!
      if (p .ge. pt(1) .and. p .le. pt(nmx)) then
        n = 1
        do while (pt(n) .lt. p)
          n = n + 1
        end do
        real_indx = float(n-1) + (p-pt(n-1)) / (pt(n)-pt(n-1))
      else
        real_indx = -1.0
      end if
!
      end function real_indx
!
! --------------------------------------------------------------
!
      function inbnds (x, y, xl, yl)
!
      real x, y, xl, yl
      logical inbnds
!
      if (x .lt. 0.0 .or. y .lt. 0.0) then
        inbnds = .false.
      else
        i = int(x)
        j = int(y)
        if (mskt(i,j) .lt. kmin .or. mskt(i+1,j) .lt. kmin .or. &
     &         mskt(i+1,j+1) .lt. kmin .or. mskt(i,j+1) .lt. kmin) then
          inbnds = .false.
        else
          if (yl .ge. Lts .and. yl .le. Ltn) then
            inbnds = .true.
          else
            inbnds = .false.
          endif
        endif
      endif
!
      end function inbnds
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
      end program mkAsmAlt
