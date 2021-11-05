!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  GODAS_MKLSACHV
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2003-07-28
!
! ABSTRACT:  Perform quality control on subsurface temperature data
!   for use in Global Ocean Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2003-07-28  David W. Behringer
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - TEMPERATURE PROFILE DATA IN IEEE
!
!   OUTPUT FILES:
!     UNIT 51  - SYNTHETIC SALINITY PROFILE DATA IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     UNIT 80  - SCRATCH FILE FOR GODAS PROFILE DATA
!     UNIT 81  - SCRATCH FILE FOR GODAS - OBS PROFILE DATA
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - getLevData, getS, indx
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - getLevData, getS
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  21 - ERROR OPENING UNIT 12
!     COND =  22 - ERROR READING UNIT 12
!     COND =  31 - ERROR OPENING UNIT 13
!     COND =  32 - ERROR READING UNIT 13
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     Requires a GODAS time_mean file in netCDF format
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
      program mkLSAchv
!
!  mkLSAchv reads a T(z) archive and creates a comparable S(z) archive
!  via a geographically varying TS relationship computed from Levitus.
!  It uses the annual Levitus files distributed with MOM. It should be
!  run with a "clean" version of a standard T(z) archive.
!
      integer imx,jmx,kmx
      real, allocatable, dimension(:) :: xt, yt, zt
      real, allocatable, dimension(:,:,:) :: t, s
      real, allocatable, dimension(:) :: tt, st
      real, allocatable, dimension(:,:) :: b
      real dpm
      integer year, month, day, hour, kd
      real, allocatable, dimension(:) :: tr, zr, sr
      real x, y
      real dtml, s0mn
      character str*30,csign*8,sid*2,dtyp*2,qkey*1
      character*50 TFile, SFile, lvTFile, lvSFile
!
      real pacw, pace, yNP, sNP, atlw, atle, yNA, sNA, ySG, sSG
      real, allocatable, dimension(:) :: bndN, bndS
!
      call w3tagb('GODAS_MKLSACHV',2003,0164,0164,'NP23')
!
      dtml = 0.25
      s0mn = 0.0
!
      pacw=110.0
      pace=260.0
      yNP=25.0
      sNP=33.8
      atlw=280.0
      atle=355.0
      yNA=32.0
      sNA=36.0
      ySG=-30.0
      sSG=34.4
!
      call getLevData(istat)
!
      if (istat .eq. 1) then
        go to 210
      else if (istat .eq. 2) then
        go to 220
      else if (istat .eq. 3) then
        go to 310
      else if (istat .eq. 4) then
        go to 320
      end if
!
! open temperature profile file
!
      open (11, form='unformatted', status='old', &
                 & access='sequential', err=110)
!
      nprf = 0
      kdmx = 0
      do while (.true.)
        read (11, end=10, err=120) iyear,idate,csign, &
                             &    sid,dtyp,qkey,y,x,kd
        nprf = nprf + 1
        if (kd .gt. kdmx) kdmx = kd
      end do
  10  continue
!
      rewind 11
!
      allocate(tr(kdmx))
      allocate(zr(kdmx))
      allocate(sr(kdmx))
!
! open synthetic salinity profile file
!
      open (51, form='unformatted', access='sequential', err=510)
!
      nsc = 0
      do n=1,nprf
!
        read (11, err=120) iyear,idate,csign,sid,dtyp,qkey,y,x, &
                           &  kd,(zr(k),tr(k),k=1,kd)
!
        call getS
!
        if (sr(1) > s0mn) then
          write (51,err=520) iyear,idate,csign,sid,dtyp,qkey,y,x, &
                           &  kd,(zr(k),sr(k),k=1,kd)
          nsc = nsc + 1
        end if
      end do
!
      close(11)
      close(51)
!
      write (6,'(i6,a)') nsc,' salinity profiles created.'
!
      call w3tage('GODAS_MKLSACHV')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a)') 'Error opening temperature profile file'
      call w3tage('GODAS_MKLSACHV')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading temperature profile file'
      call w3tage('GODAS_MKLSACHV')
      call errexit(12)
!     call exit(12)
!
  210 write(6,'(a)') 'Error opening Levitus temperature file'
      call w3tage('GODAS_MKLSACHV')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a)') 'Error reading Levitus temperature file'
      call w3tage('GODAS_MKLSACHV')
      call errexit(22)
!     call exit(22)
!
  310 write(6,'(a)') 'Error opening Levitus salinity file'
      call w3tage('GODAS_MKLSACHV')
      call errexit(31)
!     call exit(31)
!
  320 write(6,'(a)') 'Error reading Levitus salinity file'
      call w3tage('GODAS_MKLSACHV')
      call errexit(32)
!     call exit(32)
!
  510 write(6,'(a)') 'Error opening salinity profile file'
      call w3tage('GODAS_MKLSACHV')
      call errexit(51)
!     call exit(51)
!
  520 write(6,'(a)') 'Error writing salinity profile file'
      call w3tage('GODAS_MKLSACHV')
      call errexit(52)
!     call exit(52)
!
      contains
!
! ---------------------------------------------------------------- 
!
      subroutine getLevData(stat)
!
      integer stat
      character stamp*32, str*64
!
      open(unit=12,status='old',form='UNFORMATTED', &
                            & access='sequential', err=210)
!
      read(12,err=220) str
      read(12,err=220) stamp,dpm,imx,jmx,kmx
!
      allocate(xt(imx))
      allocate(yt(jmx))
      allocate(zt(kmx))
      allocate(t(imx,jmx,kmx))
      allocate(s(imx,jmx,kmx))
      allocate(b(imx,jmx))
      allocate(tt(kmx))
      allocate(st(kmx))
      allocate(bndN(imx))
      allocate(bndS(imx))
!
      rewind(12)
!
      do k=1,kmx
        read(12,err=220) str
        read(12,err=220) stamp,dpm,i,j,km,k0,month,xt,yt,zt,b
!
        do j=1,jmx
          do i=1,imx
            t(i,j,k) = b(i,j)
          end do
        end do
!
      end do
!
      day = 15
      hour = 120
!
      close(12)
!
      open(unit=13,status='old',form='UNFORMATTED', &
                            & access='sequential', err=310)
!
      do k=1,kmx
        read(13,err=320) str
        read(13,err=320) stamp,dpm,i,j,km,k0,month,xt,yt,zt,b
!
        do j=1,jmx
          do i=1,imx
            s(i,j,k) = b(i,j)
          end do
        end do
!
      end do
!
      close(13)
!
      do k=1,kmx
        zt(k) = 0.01*zt(k)
      end do
!
      do i=1,imx
        if (xt(i) .ge. pacw .and. xt(i) .le. pace) then
          jso = indx(yNP, yt, jmx)
          do j=jso,jmx
            if (s(i,j,1) .le. sNP) then
              bndN(i) = yt(j)
              exit
            end if
          end do
        else if (xt(i) .ge. atlw .and. xt(i) .le. atle) then
          jso = indx(yNA, yt, jmx)
          do j=jso,jmx
            if (s(i,j,1) .le. sNA) then
              bndN(i) = yt(j)
              exit
            end if
          end do
        else
          bndN(i) = 30.0
        end if
      end do
!
      jno = indx(ySG, yt, jmx)
      do i=1,imx
        do j=jno,1,-1
          if (s(i,j,1) .le. sSG) then
            bndS(i) = yt(j)
            exit
          end if
        end do
      end do
!
      stat = 0
      return
!
  210 stat = 1
      return
!
  220 stat = 2
      return
!
  310 stat = 3
      return
!
  320 stat = 4
      return
!
      end subroutine getLevData
!
! ---------------------------------------------------------------- 
!
      subroutine getS
!
      integer i, j, k, kr, krp, n
      real rm, rp, dt, tx, dxm, dxp, dym, dyp
!
      if (y .lt. yt(1) .or. y .gt. yt(jmx)) then
        sr(1) = -1.0
        return
      endif
!
      i = indx(x, xt, imx)
      j = indx(y, yt, jmx)
!
      dxm = (x - xt(i)) / (xt(i+1) - xt(i))
      dxp = (xt(i+1) - x) / (xt(i+1) - xt(i))
      dym = (y - yt(j)) / (yt(j+1) - yt(j))
      dyp = (yt(j+1) - y) / (yt(j+1) - yt(j))
!
      do k=1,kmx
        tt(k) = t(i,j,k)*dxp*dyp + t(i+1,j,k)*dxm*dyp + &
              & t(i,j+1,k)*dxp*dym + t(i+1,j+1,k)*dxm*dym
        st(k) = s(i,j,k)*dxp*dyp + s(i+1,j,k)*dxm*dyp + &
              & s(i,j+1,k)*dxp*dym + s(i+1,j+1,k)*dxm*dym
      end do
!
      if (y .lt. bndS(i) .or. y .gt. bndN(i)) then
        do k=1,kd
          if (zr(k) .le. zt(1)) then
            sr(k) = st(1)
          else if (zr(k) .ge. zt(kmx)) then
            sr(k) = st(kmx)
          else
            do n=2,kmx
              if (zr(k) .ge. zt(n-1) .and. zr(k) .le. zt(n)) then
                kr = n
                exit
              end if
            end do
            rm = (zt(kr) - zr(k)) / (zt(kr) - zt(kr-1))
            rp = (zr(k) - zt(kr-1)) / (zt(kr) - zt(kr-1))
            sr(k) = rm * st(kr-1) + rp * st(kr)
          end if
        end do
      else
        tx = tr(1)
        kr = -1
        do k=2,kd
          dt = tr(1) - tr(k)
          if (dt .le. dtml) then
            kr = k
            if (tr(k) .gt. tx) tx = tr(k)
          else
            exit
          end if
        end do
        do k=1,kr
          tr(k) = tx
        end do
!
        do k=1,kd
          kr = 0
          krp = 0
          do n=2,kmx
            if (tr(k) .le. tt(n-1) .and. tr(k) .ge. tt(n)) then
              kr = n
            end if
            if (tr(k) .ge. tt(n-1) .and. tr(k) .le. tt(n)) then
              krp = n
            end if
          end do
          if (kr .gt. 0 .and. krp .gt. 0) then
            if (abs(zr(k) - zt(krp)) .lt. abs(zr(k) - zt(kr))) then
              kr = krp
            end if
          else if (krp .gt. 0) then
            kr = krp
          end if
          if (kr .gt. 0) then
            rm = (tt(kr) - tr(k)) / (tt(kr) - tt(kr-1))
            rp = (tr(k) - tt(kr-1)) / (tt(kr) - tt(kr-1))
            sr(k) = rm * st(kr-1) + rp * st(kr)
          else
            if (abs(tr(k) - tt(1)) < abs(tr(k) - tt(kmx))) then
              sr(k) = st(1)
            else
              sr(k) = st(kmx)
            end if
          end if
        end do
      end if
!
      end subroutine getS
!
! ---------------------------------------------------------------- 
!
      function indx(p, pt, nmx)
!
      integer indx
      real p, pt(*)
      integer nmx, n
!
      if (p .ge. pt(1) .and. p .le. pt(nmx)) then
        n = 1
        do while (pt(n) .lt. p .and. n .lt. nmx)
          n = n + 1
        end do
        indx = n - 1
      else
        indx = -1
      end if
!
      end function indx
!
      end program mkLSAchv
