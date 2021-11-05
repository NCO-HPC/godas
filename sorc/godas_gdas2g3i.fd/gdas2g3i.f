!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  GDAS2G3I
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2003-07-28
!
! ABSTRACT:  Interpolate GDAS forcing field to ocean grid for use in 
!            Global Ocean Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2003-07-28  David W. Behringer
! 2004-04-21  David W. Behringer - correct time stamp in output file
! 2004-05-26  Diane C. Stokes    - nrec changed from 3 (1-day run) to
!                                  vary according to input parm
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - FORCING FIELD FROM GDAS
!     UNIT 12  - 2ND FORCING FIELD FROM GDAS IF NEEDED
!     UNIT 13  - MASK FOR ATMOSPHERE MODEL
!     UNIT 14  - MASK/GRID FOR OCEAN MODEL
!     UNIT 31  - NAMELIST cntrlatm
!
!   OUTPUT FILES:
!     UNIT 51  - FORCING FIELD INTERPOLATED TO OCEAN GRID
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     None
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:  - jlnDayR, jlnDayC, calDayR, expndLndMsk, expndOcnFld, intrp2d
!     LIBRARY:
!       W3LIB  - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:  - jlnDayR, calDayR, expndLndMsk, expndOcnFld, intrp2d
!     LIBRARY:
!       W3LIB  - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  21 - ERROR OPENING UNIT 12
!     COND =  22 - ERROR READING UNIT 12
!     COND =  31 - ERROR OPENING UNIT 13
!     COND =  32 - ERROR READING UNIT 13
!     COND =  41 - ERROR OPENING UNIT 14
!     COND =  42 - ERROR READING UNIT 14
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!     COND =  61 - ERROR READING COMMAND LINE
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     Requires command line arguments of date (YYYYMMDD) and 
!     two character field code that is one of the following:
!      sh  -  sensible heat
!      lh  -  latent heat
!      sw  -  short wave
!      lw  -  long wave
!      sf  -  salt flux (requires 2 fields: precip & latent heat)
!      tx  -  zonal wind stress
!      ty  -  meridional wind stress
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program gdas2g3i
!
      use gnutil_mod
      use tmutil_mod
!
      integer, parameter :: imx=362, jmx=202, kmx=40, nOvr=2
      integer, parameter :: img=360, jmg=181
      character(len=80) :: str
      real, dimension(imx,jmx) :: at, a2
      real, dimension(imx) :: xt
      real, dimension(jmx) :: yt
      real, dimension(kmx) :: zt, zw
      integer, dimension(imx,jmx) :: mskt
      real, dimension(img,jmg) :: ag, wrk
      real, dimension(img) :: xg
      real, dimension(jmg) :: yg
      integer, dimension(img,jmg) :: mask, msk0, msk1
      integer :: year, month, day, wk, yr, mo, dy, yr0
      integer, dimension(:), allocatable :: jdy
      integer :: jdyM
      integer :: yrs, mos, dys, jdys
      logical :: data_found, two, lst_src_dtes
      character(len=80) :: iotext
      character(len=14) :: iot0
      data iot0 /'dim x(i),y(j),'/
      character(len=36) :: iot1
      data iot1 /'(i,j);read(nu)stamp,avg,i,j,mon,x,y,'/
      character(len=26) :: iot2
      data iot2 /'                          '/
      character(len=2) :: iotx
      character(len=32) :: stamp
      character(len=6) :: st0
      data st0 /'m/d/y='/
      character(len=16) :: st1
      data st1 /', h:m:s= 0: 0: 0'/
!
      integer :: i, j, k, n, len, ln1, narg
      real :: scale, scal2, depth
      real :: prd1 = 1.0
!
      namelist /cntrlatm/ ndysbc
!
      call w3tagb('GODAS_GDAS2G3I',2003,0164,0164,'NP23')
!
      narg = iargc()
      if (narg .lt. 2) go to 610
!
      two = .false.
      call getarg(1,str)
      read(str,'(i4,2i2)',err=610) year, month, day
      call getarg(2,str)
      if (str(1:2) .eq. 'sh' .or. str(1:2) .eq. 'SH') then
        iotx = 'sh'
        scale = 2.39e-5
      else if (str(1:2) .eq. 'lh' .or. str(1:2) .eq. 'LH') then
        iotx = 'lh'
        scale = 2.39e-5
      else if (str(1:2) .eq. 'lw' .or. str(1:2) .eq. 'LW') then
        iotx = 'lw'
        scale = 2.39e-5
      else if (str(1:2) .eq. 'sw' .or. str(1:2) .eq. 'SW') then
        iotx = 'sw'
        scale = 2.39e-5
      else if (str(1:2) .eq. 'sf' .or. str(1:2) .eq. 'SF') then
        iotx = 'sw'
        scale = -1.4e-9
        scal2 = -3.5e-3
        two = .true.
      else if (str(1:2) .eq. 'tx' .or. str(1:2) .eq. 'TX') then
        iotx = 'tx'
        scale = 10.0
      else if (str(1:2) .eq. 'ty' .or. str(1:2) .eq. 'TY') then
        iotx = 'ty'
        scale = 10.0
      else
        write (6,'(a)')  &
                &  'Field code must be sh, lh, sw, lw, sf, tx or ty'
        go to 610
      end if
!
!  Read atmosphere mask file
!
      open(unit=13,status='old',access='SEQUENTIAL', &
                    & form='UNFORMATTED', err=310)
      read(13, err=320) mask
      close(13)
      call expndLndMsk(mask,img,jmg,2,msk0)
!
      dx = 360./float(img)
      do i=1,img
        xg(i)=float(i-1)*dx
      enddo
      dy = 1.0
      do j=1,jmg
        yg(j)=float(j-1)*dy-90.
      enddo
!
!  Read ocean mask file
!
      open(unit=14,status='old',access='SEQUENTIAL', &
                    & form='UNFORMATTED', err=410)
      read(14, err=420) i,j,k
      read(14, err=420) xt, yt, zt, zw
      read(14, err=420) mskt
      close(14)
!
!  Open atmospheric forcing file
!
      open(unit=11,status='old',access='SEQUENTIAL', &
                    & form='UNFORMATTED', err=110)
!
      ndysbc = 1
      read  (31, cntrlatm)
      nrec = ndysbc + 2
      allocate (jdy(nrec))
!
      yr0 = 9999
      do n = 1,nrec
        read(11, err=120) yr
        if (yr .lt. yr0) yr0 = yr
      enddo
      rewind(11)
      jdyM = jlnDayR(year,month,day,yr0)
      do n = 1,nrec
        read(11, err=120) yr,mo,dy
        jdy(n) = jlnDayR(yr,mo,dy,yr0)
      enddo
      rewind(11)
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
!  Open 2nd atmospheric forcing file if necessary
!
      if (two) then
        open(unit=12,status='old',access='SEQUENTIAL', &
                    & form='UNFORMATTED', err=210)
      endif
!
      open(unit=51,access='SEQUENTIAL',form='UNFORMATTED',err=510)
      write(iotext,'(a,a,a,a,a)') iot0,iotx,iot1,iotx,iot2
!
      do n = 1,nrec
        read(11, err=120) yr,mo,dy,wk,ag
        call expndOcnFld(ag,img,jmg,xg,yg,mask,-5,wrk,msk0,msk1)
        call intrp2d(ag,xg,yg,img,jmg,at,xt,yt,imx,jmx,nOvr)
        if (two) then
          read(12, err=220) yr,mo,dy,wk,ag
          call expndOcnFld(ag,img,jmg,xg,yg,mask,-5,wrk,msk0,msk1)
          call intrp2d(ag,xg,yg,img,jmg,a2,xt,yt,imx,jmx,nOvr)
          at = scale*at + scal2*a2
        else
          at = scale*at
        endif
        jdys = jlnDayR(yr,mo,dy,yr0) + 1
        call calDayR(yrs,mos,dys,yr0,jdys)
        write(6,'(i4,a,i2.2,a,i2.2)') yrs,'-',mos,'-',dys
        write(stamp,'(a,i2,a,i2,a,i4.4,a)') st0,mos,'/',dys,'/',yrs,st1
        write(51, err=520) iotext
        write(51, err=520) stamp,prd1,imx,jmx,month,xt,yt,at
      enddo
!
      close(11)
      if (two) close(12)
      close(51)
!
      call w3tage('GODAS_GDAS2G3I')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a)') 'Error opening GDAS forcing file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading GDAS forcing file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(12)
!     call exit(12)
!
  210 write(6,'(a)') 'Error opening 2nd GDAS forcing file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a)') 'Error reading 2nd GDAS forcing file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(22)
!     call exit(22)
!
  310 write(6,'(a)') 'Error opening atmosphere mask file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(31)
!     call exit(31)
!
  320 write(6,'(a)') 'Error reading atmosphere mask file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(32)
!     call exit(32)
!
  410 write(6,'(a)') 'Error opening ocean mask file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(41)
!     call exit(41)
!
  420 write(6,'(a)') 'Error reading ocean mask file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(42)
!     call exit(42)
!
  510 write(6,'(a)') 'Error opening ocean forcing file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(51)
!     call exit(51)
!
  520 write(6,'(a)') 'Error writing ocean forcing file'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(52)
!     call exit(52)
!
  610 write(6,'(a)') 'Error reading command line'
      call w3tage('GODAS_GDAS2G3I')
      call errexit(61)
!     call exit(61)
!
      end program gdas2g3i
