!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  SST2G3I
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2003-07-28
!
! ABSTRACT:  Interpolate SST analysis to ocean grid for use in 
!            Global Ocean Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2003-07-28  David W. Behringer
! 2004-04-21  David W. Behringer - correct time stamp in output file
! 2004-05-26  Diane C. Stokes    - nrec changed from 3 (1-day run) to
!                                  vary according to input parm
! 2004-08-18  Diane C. Stokes    - correct starting lat of input sst 
!                                  (yg0 changed from -89.0 to -89.5)
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - SST ANALYSIS
!     UNIT 12  - SST INFO FILE FROM GRIB
!     UNIT 14  - MASK/GRID FOR OCEAN MODEL
!     UNIT 31  - NAMELIST cntrlatm
!
!   OUTPUT FILES:
!     UNIT 51  - SST FIELD INTERPOLATED TO OCEAN GRID
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     None
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:  - rvrs_ns, shft_180, jlnDayR, jlnDayC, calDayR, intrp2d
!     LIBRARY:
!       W3LIB  - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:  - rvrs_ns, shft_180, jlnDayR, calDayR, intrp2d
!     LIBRARY:
!       W3LIB  - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  21 - ERROR OPENING UNIT 12
!     COND =  22 - ERROR READING UNIT 12
!     COND =  41 - ERROR OPENING UNIT 14
!     COND =  42 - ERROR READING UNIT 14
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!     COND =  61 - ERROR READING COMMAND LINE
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     Requires command line argument of date (YYYYMMDD)
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program sst2g3i
!
      use gnutil_mod
      use tmutil_mod
!
      integer, parameter :: imx=362, jmx=202, kmx=40, nOvr=2
      integer, parameter :: img=360, jmg=180
      character(len=80) :: str
      real, dimension(imx,jmx) :: at, a2
      real, dimension(imx) :: xt
      real, dimension(jmx) :: yt
      real, dimension(kmx) :: zt, zw
      integer, dimension(imx,jmx) :: mskt
      real, dimension(img,jmg) :: ag, wrk
      real, parameter :: xg0 = 0.5, yg0 = -89.5, tk0 = 273.15
      real, dimension(img) :: xg, wg
      real, dimension(jmg) :: yg
      integer, dimension(img,jmg) :: mask, msk0, msk1
      integer :: year, month, day
      integer :: yr, mo, dy, ndys, indx, yr0
      integer, dimension(:), allocatable :: jdy
      integer :: oset, jdyM, jdys
      logical :: data_found, lst_src_dtes
      character(len=80) :: iotext
      character(len=14) :: iot0
      data iot0 /'dim x(i),y(j),'/
      character(len=40) :: iot1
      data iot1 /'(i,j);read(nu)stamp,avg,i,j,k,z,mon,x,y,'/
      character(len=22) :: iot2
      data iot2 /'                      '/
      character(len=2) :: iotx
      data iotx /'st'/
      character(len=32) :: stamp
      character(len=6) :: st0
      data st0 /'m/d/y='/
      character(len=16) :: st1
      data st1 /', h:m:s= 0: 0: 0'/
      integer :: km = 1
      real :: z0 = 0.0
      real :: prd1 = 1.0
!
      integer :: i, j, k, n, len, ln1, narg
!
      namelist /cntrlatm/ ndysbc
!
      call w3tagb('GODAS_SST2G3I',2003,0164,0164,'NP23')
!
      narg = iargc()
      if (narg .lt. 1) go to 610
!
      call getarg(1,str)
      read(str,'(i4,2i2)',err=610) year, month, day
!
!  Open SST info file
!
      open(unit=12,status='old',access='SEQUENTIAL', &
                        & form='FORMATTED', err=210)
!
      ndysbc = 1
      read  (31, cntrlatm)
      nrec = ndysbc + 2
      allocate (jdy(nrec))
!
      yr0 = 9999
      do n = 1,nrec
        read(12,'(a)', err=220) str
        do i=1,79
          if (str(i:i+1) .eq. 'd=') then
            read(str(i+2:i+5), '(i4)') yr
            exit
          endif
        enddo
        if (yr .lt. yr0) yr0 = yr
      enddo
      rewind(12)
      jdyM = jlnDayR(year,month,day,yr0)
      do n = 1,nrec
        read(12,'(a)', err=220) str
        do i=1,79
          if (str(i:i+1) .eq. 'd=') then
            read(str(i+2:i+9), '(i4,i2,i2)') yr, mo, dy
            exit
          endif
        enddo
        jdy(n) = jlnDayR(yr,mo,dy,yr0)
      enddo
      close(12)
      if (jdy(1) .le. jdyM .and. jdy(nrec) .ge. jdyM+nrec-3) then
        data_found = .true.
      else
        data_found = .false.
      endif
!
      if (.not. data_found) then
        write(6,'(a,i2,a,i4,a,i2.2,a,i2.2)') 'Forcing data unavailable&
                  & for ', nrec-2, ' days starting ', year, '-', &
                  & month, '-', day
        go to 120
      endif
!
!  Compute SST analysis grid
!
      dx = 1.0
      do i=1,img
        xg(i) = xg0 + float(i-1)*dx
      enddo
      dy = 1.0
      do j=1,jmg
        yg(j) = yg0 + float(j-1)*dy
      enddo
!
!  Open ocean grid/mask file
!
      open(unit=14,status='old',access='SEQUENTIAL', &
                    & form='UNFORMATTED', err=410)
      read(14, err=420) i,j,k
      read(14, err=420) xt, yt, zt, zw
      read(14, err=420) mskt
      close(14)
!
!  Open SST analysis file
!
      open(unit=11,status='old',access='SEQUENTIAL', &
                    & form='UNFORMATTED', err=110)
!
!  Open SST ocean file
!
      open(unit=51,access='SEQUENTIAL',form='UNFORMATTED',err=510)
      write(iotext,'(a,a,a,a,a)') iot0,iotx,iot1,iotx,iot2
!
      do n = 1,nrec
        read(11, err=120) ag
        call rvrs_ns(ag,img,jmg,wrk)
!
!       call shft_180(ag,img,jmg,wg)
        call intrp2d(ag,xg,yg,img,jmg,at,xt,yt,imx,jmx,nOvr)
        jdys = jdy(n) + 1
        call calDayR(yr,mo,dy,yr0,jdys)
        write(6,'(i4,a,i2.2,a,i2.2)') yr,'-',mo,'-',dy
        write(stamp,'(a,i2,a,i2,a,i4.4,a)') st0,mo,'/',dy,'/',yr,st1
        write(51, err=520) iotext
        write(51, err=520) stamp,prd1,imx,jmx,km,z0,month,xt,yt,at
      enddo
!
      close(11)
      close(51)
!
      call w3tage('GODAS_SST2G3I')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a)') 'Error opening SST analysis file'
      call w3tage('GODAS_SST2G3I')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading SST analysis file'
      call w3tage('GODAS_SST2G3I')
      call errexit(12)
!     call exit(12)
!
  210 write(6,'(a)') 'Error opening SST info file'
      call w3tage('GODAS_SST2G3I')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a)') 'Error reading SST info file'
      call w3tage('GODAS_SST2G3I')
      call errexit(22)
!     call exit(22)
!
  410 write(6,'(a)') 'Error opening ocean mask file'
      call w3tage('GODAS_SST2G3I')
      call errexit(41)
!     call exit(41)
!
  420 write(6,'(a)') 'Error reading ocean mask file'
      call w3tage('GODAS_SST2G3I')
      call errexit(42)
!     call exit(42)
!
  510 write(6,'(a)') 'Error opening ocean SST file'
      call w3tage('GODAS_SST2G3I')
      call errexit(51)
!     call exit(51)
!
  520 write(6,'(a)') 'Error writing ocean SST file'
      call w3tage('GODAS_SST2G3I')
      call errexit(52)
!     call exit(52)
!
  610 write(6,'(a)') 'Error reading command line'
      call w3tage('GODAS_SST2G3I')
      call errexit(61)
!     call exit(61)
!
    contains
!
      subroutine rvrs_ns(a,im,jm,w)
!
      integer :: im, jm
      real, dimension(im,jm) :: a
      real, dimension(im,jm) :: w
      integer :: i, j, jj
!
      do j=1,jm
        jj = jm - j + 1
        do i=1,im
          w(i,jj) = a(i,j)
        enddo
      enddo
!
      do j=1,jm
        do i=1,im
          a(i,j) = w(i,j) - tk0
        enddo
      enddo
!
      end subroutine rvrs_ns
!

      subroutine shft_180(a,im,jm,w)
!
      integer :: im, jm
      real, dimension(im,jm) :: a
      real, dimension(im) :: w
      integer :: i, ii, j
!
      do j=1,jm
        do i=1,im
          ii = i+180
          if (ii > im) ii = ii - im
          w(ii) = a(i,j)
        enddo
        do i=1,im
          a(i,j) = w(i)
        enddo
      enddo
!
      end subroutine shft_180
!
    end program sst2g3i
