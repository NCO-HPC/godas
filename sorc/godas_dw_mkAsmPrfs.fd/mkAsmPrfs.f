!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  GODAS_MKASMPRFS
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2003-07-25
!
! ABSTRACT:  Write salinity observations in form used by GODAS
!
! PROGRAM HISTORY LOG:
! 2003-07-25  David W. Behringer
! 2004-09-13  David W. Behringer - Fix of function real_indx:
!              real_indx = float(n) + ... is changed to
!              real_indx = float(n-1) + ...
! 2007-08-28  David W. Behringer - Extend profiles from 30 levels max
!              to 35 levels, using climatological temperature data.
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - SYNTHETIC SALINITY PROFILE DATA IN IEEE
!                 BASED TEMPERATURE PROFILES ALREADY FILTERED TO
!                 BE INBOUNDS, Z-INTERPOLATED AND Z-EXTENDED
!     UNIT 12  - GRID/MASK FOR OCEAN MODEL
!
!   OUTPUT FILES:
!     UNIT 51  - SALINITY PROFILE DATA READY FOR ASSIMILATION IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     UNIT 80  - SCRATCH FILE FOR PROFILE DATA
!     UNIT 81  - SCRATCH FILE FOR PROFILE DATA
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - real_indx, qSort, swap
!                  jlnDayR, jlnDayC, calDayR, calDayC, weekDay
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - real_indx, qSort, swap
!                  jlnDayR, calDayR, weekDay
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  21 - ERROR OPENING UNIT 12
!     COND =  22 - ERROR READING UNIT 12
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!     COND =  61 - ERROR READING FROM COMMAND LINE
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     Reads a date YYYYMMDD from command line
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program mkAsmPrfs
!
!  mkAsmPrfs prepares profiles for use in assimilation
!
      use tmutil_mod
!
      integer imx, jmx, kmx, kax
      character csign*8,sid*2,dtyp*2,qkey*1,csignS*8
      character str*80, cdate*8
      real, allocatable, dimension(:) :: pt, pz, ptS, pzS
      real, allocatable, dimension(:) :: xt, yt, zt, zw, dz
      real stde, evr
      integer, allocatable, dimension(:,:) :: mskt
      integer ayear, amonth, aday, dowA, dow0, dow
      character*32 stamp
      logical okay
      integer(8), allocatable, dimension(:) :: ndxp
      integer, allocatable, dimension(:) :: nprfl
!
      call w3tagb('GODAS_MKASMPRFS',2003,0164,0164,'NP23')
!
!     kax = maximum number of model levels for assimilation
!     nwks = time window (weeks) for profile data
!
      kax = 35
      nwks = 5
!
!  set std errors, then 1/e.v. = 1/(se*se)
!  if std error = 0.5, then EVR = 4.0
!                 1.0             1.0
!                 1.5             0.4444
!  multiplying stde by 0.001 rescales 0/00 to psu
!
      stde = 0.1
      stde = 0.001 * stde
      evr = 1.0/(stde*stde)
!
! get date from command line
!
      narg = iargc()
      if (narg .eq. 0) go to 610
      call getarg(1,str)
      read(str,'(i4,2i2)',err=610) ayear, amonth, aday
!
! read model coordinates and mask
!
      open(unit=12,status='old',form='UNFORMATTED', &
                 & access='sequential', err=210)
      read(12, err=220) imx,jmx,kmx
      allocate(xt(imx))
      allocate(yt(jmx))
      allocate(zt(kmx))
      allocate(zw(kmx))
      allocate(mskt(imx,jmx))
      read(12, err=220) xt, yt, zt, zw
      read(12, err=220) mskt
      close(12)
!
      allocate(dz(kmx))
      dz(1) = 2.0*zt(1)
      do k=2,kmx
        dz(k) = 2.0*(zt(k) - zt(k-1)) - dz(k-1)
      end do
!
! set time limits for profile data
!
      iyear0 = ayear - 1
      dowA = weekDay(ayear, amonth, aday)
      jdayA = jlnDayR(ayear, amonth, aday, iyear0)
      jdayS = jdayA - dowA - 7*(nwks/2) + 1
      jdayE = jdayS + 7*nwks - 1
      call calDayR(iyrS,imoS,idyS,iyear0,jdayS)
      call calDayR(iyrE,imoE,idyE,iyear0,jdayE)
      write (6,'(a,i4,a,i2,a,i2,a,i4,a,i2,a,i2)') &
       & 'Profiles will be retained between ',iyrS,'/',imoS,'/',idyS, &
       & ' and ',iyrE,'/',imoE,'/',idyE
!
! open assimilation profile file
!
      open (51,form='unformatted',access='sequential',err=510)
!
! open the input profile file
!
      open (11,status='old',form='unformatted', &
                 & access='sequential', err=110)
!
      nprf = 0
      npmx = 0
      do while (.true.)
        read(11,end=5,err=120) iyear,idate,csign,sid, &
                                       &  dtyp,qkey,yp,xp,np
        nprf = nprf + 1
        if (np .gt. npmx) npmx = np
      end do
    5 continue
      rewind (11)
!
      allocate(pt(npmx))
      allocate(pz(npmx))
      allocate(ptS(npmx))
      allocate(pzS(npmx))
      allocate(ndxp(nprf))
      allocate(nprfl(nprf))
!
! open direct access scratch file and a sequential scratch file
!
      nb = 4*(2*npmx + 5) + 8
      open (80,status='scratch',form='unformatted', &
                               & access='direct',recl=nb)
!
      open (81,status='scratch',form='unformatted', &
                               & access='sequential')
!
! get reference dates for first week
!
      read (11, err=120) iyear,idate,csign,sid,dtyp, &
                           &  qkey,yp,xp,np,(pz(k),pt(k),k=1,np)
      write (cdate,'(i8.8)') idate
      read (cdate,'(4i2)') imon, iday, ihr, imin
      jday0 = jlnDayR(iyear, imon, iday, iyear0)
!
      write (cdate,'(i8.8)') idate
      read (cdate,'(4i2)') imon, iday, ihr, imin
      dow0 = weekDay(iyear, imon, iday)
!
      jday0 = jlnDayR(iyear, imon, iday, iyear0)
      jday1 = jday0 - dow0 + 1
      jday2 = jday1 + 6
      jdayR = jday1 + 3
!
      xi = real_indx(xp,xt,imx)
      yj = real_indx(yp,yt,jmx)
      nrc = 1
      ndxp(nrc) = int(yj*1000.0)
      ndxp(nrc) = ndxp(nrc)*1000000 + int(xi*1000.0)
      nprfl(nrc) = nrc
      write (80,rec=nrc) iyear,idate,csign,xi,yj,np,pz,pt
!
! loop on profile file
!
      do while (.true.)
        read (11, end=100, err=120) iyear,idate,csign,sid,dtyp, &
                           &  qkey,yp,xp,np,(pz(k),pt(k),k=1,np)
        write (cdate,'(i8.8)') idate
        read (cdate,'(4i2)') imon, iday, ihr, imin
        dow = weekDay(iyear, imon, iday)
        jday = jlnDayR(iyear, imon, iday, iyear0)
!
        if (jday .lt. jday0) then
          write(6,'(a)') 'Reading profile file (unit 11) :'
          write(6,'(a)') ' time sequence error.'
          go to 120
        end if
!
        if (jday .le. jday2) then
          xi = real_indx(xp,xt,imx)
          yj = real_indx(yp,yt,jmx)
          nrc = nrc + 1
          ndxp(nrc) = int(yj*1000.0)
          ndxp(nrc) = ndxp(nrc)*1000000 + int(xi*1000.0)
          nprfl(nrc) = nrc
          write (80,rec=nrc) iyear,idate,csign,xi,yj,np,pz,pt
        else
!
! save data
!
          iyearS = iyear
          idateS = idate
          csignS = csign
          nps = np
          do k=1,np
            pzS(k) = pz(k)
            ptS(k) = pt(k)
          end do
!
! sort on ndxp
!
          ns = 1
          ne = nrc
          call qSort(ndxp,nprfl,ns,ne)
!
          nobs = 0
          do n=1,nrc
            read (80,rec=nprfl(n)) iyear,idate,csign,xi,yj,np,pz,pt
            write (cdate,'(i8.8)') idate
            read (cdate,'(4i2)') imon, iday, ihr, imin
            dow = weekDay(iyear, imon, iday)
            fdow = float(dow-1) + (float(ihr) + float(imin)/60.0)/24.0
            rday = fdow - 3.5
            kt = np
!
            i = xi
            if (i .lt. 1) i = 1
            j = yj
            km = min(mskt(i,j),mskt(i,j+1),mskt(i+1,j+1),mskt(i+1,j))
            kt = min(kt,km,kax)
!
            do k=1,kt
              zk = float(k)
              ttk = 0.001 * (pt(k) - 35.0)
              write (81) csign,rday,ttk,xi,yj,zk,evr
              nobs = nobs + 1
            end do
          end do
!
          rewind (81)
!
          call calDayR(iyr,imo,idy,iyear0,jdayR)
          write (stamp,'(a,i2,a,i2,a,i4,a)') &
                 & 'm/d/y=',imo,'/',idy,'/',iyr,', h:m:s=12: 0: 0'
          write(6,'(a,i8)') stamp,nobs
          write(51,err=520) stamp,nobs
          do n=1,nobs
            read(81) csign,rday,ttk,xi,yj,zk,evr
            write(51,err=520) csign,rday,ttk,xi,yj,zk,evr
          end do
!
          rewind (81)
!
          jday1 = jday1 + 7
          jday2 = jday1 + 6
          jdayR = jday1 + 3
!
          nrc = 0
          xi = real_indx(xp,xt,imx)
          yj = real_indx(yp,yt,jmx)
          nrc = nrc + 1
          ndxp(nrc) = int(yj*1000.0)
          ndxp(nrc) = ndxp(nrc)*1000000 + int(xi*1000.0)
          nprfl(nrc) = nrc
          write (80,rec=nrc) iyearS,idateS,csignS,xi,yj,npS,pzS,ptS
!
        end if
!
      end do
  100 continue
!
      close (11)
!
! sort on ndxp
!
      ns = 1
      ne = nrc
      call qSort(ndxp,nprfl,ns,ne)
!
      nobs = 0
      do n=1,nrc
        read (80,rec=nprfl(n)) iyear,idate,csign,xi,yj,np,pz,pt
        write (cdate,'(i8.8)') idate
        read (cdate,'(4i2)') imon, iday, ihr, imin
        dow = weekDay(iyear, imon, iday)
        fdow = float(dow-1) + (float(ihr) + float(imin)/60.0)/24.0
        rday = fdow - 3.5
        kt = np
!
        i = xi
        if (i .lt. 1) i = 1
        j = yj
        km = min(mskt(i,j),mskt(i,j+1),mskt(i+1,j+1),mskt(i+1,j))
        kt = min(kt,km,kax)
!
        do k=1,kt
          zk = float(k)
          ttk = 0.001 * (pt(k) - 35.0)
          write (81) csign,rday,ttk,xi,yj,zk,evr
          nobs = nobs + 1
        end do
      end do
!
      close (80)
      rewind (81)
!
      call calDayR(iyr,imo,idy,iyear0,jdayR)
      write (stamp,'(a,i2,a,i2,a,i4,a)') &
             & 'm/d/y=',imo,'/',idy,'/',iyr,', h:m:s=12: 0: 0'
      write(6,'(a,i8)') stamp,nobs
      write(51,err=520) stamp,nobs
      do n=1,nobs
        read(81) csign,rday,ttk,xi,yj,zk,evr
        write(51,err=520) csign,rday,ttk,xi,yj,zk,evr
      end do
!
      close (81)
      close (51)
!
      call w3tage('GODAS_MKASMPRFS')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a)') 'Error opening profile file on unit 11'
      call w3tage('GODAS_MKASMPRFS')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading profile file on unit 11'
      call w3tage('GODAS_MKASMPRFS')
      call errexit(12)
!     call exit(12)
!
  210 write(6,'(a)') 'Error opening grid/mask file on unit 12'
      call w3tage('GODAS_MKASMPRFS')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a)') 'Error reading grid/mask file on unit 12'
      call w3tage('GODAS_MKASMPRFS')
      call errexit(22)
!     call exit(22)
!
  510 write(6,'(a)') 'Error opening profile file on unit 51'
      call w3tage('GODAS_MKASMPRFS')
      call errexit(51)
!     call exit(51)
!
  520 write(6,'(a)') 'Error writing profile file on unit 51'
      call w3tage('GODAS_MKASMPRFS')
      call errexit(52)
!     call exit(52)
!
  610 write(6,'(a)') 'Error reading date from command line'
      call w3tage('GODAS_MKASMPRFS')
      call errexit(61)
!     call exit(61)
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
      if (p .gt. pt(1) .and. p .le. pt(nmx)) then
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
      end program mkAsmPrfs
