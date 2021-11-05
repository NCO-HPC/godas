!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  GODAS_MKASMPRF
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2003-07-25
!
! ABSTRACT:  Write temperature observations in form used by GODAS
!
! PROGRAM HISTORY LOG:
! 2003-07-25  David W. Behringer
! 2004-09-13  David W. Behringer - Fix of function real_indx:
!              real_indx = float(n) + ... is changed to
!              real_indx = float(n-1) + ...
! 2015-09-29  David W. Behringer - If a platform reports mutltiple
!              profiles on a single day, they are averaged
! 2017-09-15  David W. Behringer - Allows daily averages of massive clusters
!              of profiles, e.g. those air-dropped into hurricanes
! 2018-01-15  David W. Behringer - Fix bug causing occasional zero-filledi
!              profiles
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - TEMPERATURE PROFILE DATA IN IEEE
!     UNIT 12  - GRID/MASK FOR OCEAN MODEL
!
!   OUTPUT FILES:
!     UNIT 51  - TEMPERATURE PROFILE DATA FOR ASSIMILATION IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     UNIT 80  - SCRATCH FILE FOR PROFILE DATA
!     UNIT 81  - SCRATCH FILE FOR PROFILE DATA
!     UNIT 82  - SCRATCH FILE FOR PROFILE DATA
!     UNIT 83  - SCRATCH FILE FOR FILTERED PROFILE DATA
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - real_indx, inbnds, intrpZ, qSort, swap, ssort
!                  jlnDayR, jlnDayC, calDayR, calDayC, weekDay
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - real_indx, inbnds, intrpZ, qSort, swap, ssort
!                  jlnDayR, calDayR, weekDay
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  13 - ERROR READING UNIT 83
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
      program mkAsmPrf
!
!  mkAsmPrf prepares profiles for use in assimilation
!
      use tmutil_mod
!
      integer imx, jmx, kmx
      character csign*8,sid*2,dtyp*2,qkey*1,csign0*8
      character str*80, cdate*8
      real, allocatable, dimension(:) :: pt, pz, w
      real, allocatable, dimension(:) :: xt, yt, zt, zw, dz, tt, ttz
      real, allocatable, dimension(:) :: sz, st, stz
      integer, allocatable, dimension(:) :: cnt
      integer, allocatable, dimension(:,:) :: mskt
      integer ayear, amonth, aday, dowA, dow0, dow
      character*32 stamp
      logical okay, sclDz, srTz
      real z0g, dZg, tzeps, mVal
      integer kav
      integer(8), allocatable, dimension(:) :: ndxp, ndxd
      integer, allocatable, dimension(:) :: nprfl
!
      call w3tagb('GODAS_MKASMPRF',2003,0164,0164,'NP23')
!
!     kax = maximum number of model levels for assimilation
!     nwks = time window (weeks) for profile data
!
      kax = 30
      nwks = 5
!
!  to protect against linear interpolation across large data gaps
!  as might occur with TAO dropouts, set a threshold dZ beyond which
!  a profile is truncated if it is encountered shallower than Z0
!  (applied in intrpZ)
!
      z0g = 300.0
      dZg = 300.0
!
      sclDz = .true.
      srTz = .true.
      kav = 5
      tzeps =  0.00005
!
!  to set std errors between 1.0 and 2.5
!    SE0 = 1.0 and SEF = 1.5
!  then 1/e.v. = 1/(se*se)
!  if std error = 0.5, then EVR = 4.0
!                 1.0             1.0
!                 1.5             0.4444
!
      se0 = 1.0
      seF = 1.5
! set missing value
      mVal = 999.999
!
! get date from command line
!
      narg = iargc()
      if (narg .eq. 0) go to 610
      call getarg(1,str)
      read(str,'(i4,2i2)',err=610) ayear, amonth, aday
!
! read model coordinates and mask
!  imx = 362, jmx = 202, kmx = 40
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
      allocate(tt(kmx))
      allocate(ttz(kmx))
      allocate(st(kmx))
      allocate(stz(kmx))
      allocate(sz(kmx))
!
! set time limits for profile data
!
      iyearR = ayear - 1
      dowA = weekDay(ayear, amonth, aday)
      jdayA = jlnDayR(ayear, amonth, aday, iyearR)
      jdayS = jdayA - dowA - 7*(nwks/2) + 1
      jdayE = jdayS + 7*nwks - 1
      call calDayR(iyrS,imoS,idyS,iyearR,jdayS)
      call calDayR(iyrE,imoE,idyE,iyearR,jdayE)
      write (6,'(a,i4,a,i2,a,i2,a,i4,a,i2,a,i2)') &
       & 'Profiles will be retained between ',iyrS,'/',imoS,'/',idyS, &
       & ' and ',iyrE,'/',imoE,'/',idyE
!
! begin by interpolating to model depths and sorting by call ID
! a search will be made for IDs associated with unusaully large
!  numbers of profiles
!
! open the input profile file
!
      open (11,status='old',form='unformatted', &
                 & access='sequential', err=110)
!
      nprf = 0
      npmx = 0
      do while (.true.)
        read(11,end=55,err=120) iyear,idate,csign,sid, &
                                       &  dtyp,qkey,yp,xp,np
        nprf = nprf + 1
        if (np .gt. npmx) npmx = np
      end do
   55 continue
      rewind (11)
!
      allocate(pt(npmx))
      allocate(pz(npmx))
      allocate(w(npmx))
      allocate(cnt(npmx))
!
! open direct access scratch file and a sequential scratch file
!
      nb = 4*(3*kax + 5) + 8
!     open (80,file='scrtch0',form='unformatted', &
!                              & access='direct',recl=nb)
      open (80,status='scratch',form='unformatted', &
                               & access='direct',recl=nb)
      open (82,file='scrtch2',form='formatted',blank='zero')
!
      nrc = 0
      do while (.true.)
        read (11, end=56, err=120) iyear,idate,csign,sid,dtyp, &
                           &  qkey,yp,xp,np,(pz(k),pt(k),k=1,np)
        write (cdate,'(i8.8)') idate
        read (cdate,'(2i2)') imon, iday
        kt = 0
        if (np .gt. 3) then
          call intrpZ(kt,pt,w,pz,np,tt,ttz,zt,kmx)
        end if                                       ! bug fix
        if (kt .gt. 3) then                          ! 3/15/2018
          if (kt .gt. kax) kt = kax
          do k=kt+1,kax
            tt(k) = mVal
            ttz(k) = mVal
          enddo
          nrc = nrc+1
          write(80,rec=nrc) iyear,idate,csign,yp,xp,kt, &
                                 & (tt(k),ttz(k),zt(k),k=1,kax)
          write(82,'(a8,i4,2i2.2,i10)') csign,iyear,imon,iday,nrc
        endif
      enddo
   56 continue
      close(11)
      close(82)
!
! sort on call ID and date as listed in scrtch2 (fort.82)
!
      ko = ssort_('scrtch2')
!
      open(82,file='scrtch2',form='formatted')
      rewind(82)
      rewind(80)
      open (83,status='scratch',form='unformatted', access='sequential')
!
      read(82,'(a8,i4,2i2,i10)') csign,iyear,imon,iday,n
      read(80,rec=n) iyear,idate,csign,yp,xp,kt, &
                               & (tt(k),ttz(k),zt(k),k=1,kt)
      csign0 = csign
      idate0 = idate
      iyear0 = iyear
      imon0 = imon
      iday0 = iday
      xp0 = xp
      yp0 = yp
      do k=1,kt
        cnt(k) = 1
        st(k) = tt(k)
        stz(k) = ttz(k)
        sz(k) = zt(k)
      enddo
      do k=kt+1,kax
        cnt(k) = 0
        st(k) = 0.0
        stz(k) = 0.0
        sz(k) = 0.0
      enddo
      do i=2,nrc
        read(82,'(a8,i4,2i2,i10)') csign,iyear,imon,iday,n
        read(80,rec=n) iyear,idate,csign,yp,xp,kt, & 
                                  & (tt(k),ttz(k),zt(k),k=1,kt)
        if (csign .ne. csign0 .or. iyear .ne. iyear0  .or. imon &
                                & .ne. imon0  .or. iday .ne. iday0) then
          do k=1,kax
            if (cnt(k) .gt. 0) then
              st(k) = st(k) / real(cnt(k))
              stz(k) = stz(k) / real(cnt(k))
              kts = k
            endif
          enddo
          write(83) iyear0,idate0,csign0,yp0,xp0,kts, &
                                      & (st(k),stz(k),sz(k),k=1,kts)
          csign0 = csign
          iyear0 = iyear
          idate0 = idate
          imon0 = imon
          iday0 = iday
          xp0 = xp
          yp0 = yp
          do k=1,kt
            cnt(k) = 1
            st(k) = tt(k)
            stz(k) = ttz(k)
            sz(k) = zt(k)
          enddo
          do k=kt+1,kax
            cnt(k) = 0
            st(k) = 0.0
            stz(k) = 0.0
            sz(k) = 0.0
          enddo
        endif
        do k=1,kt
          cnt(k) = cnt(k) + 1
          st(k) = st(k) + tt(k)
          stz(k) = stz(k) + ttz(k)
          sz(k) = zt(k)
        enddo
      enddo
      do k=1,kax
        if (cnt(k) .gt. 0) then
          st(k) = st(k) / real(cnt(k))
          stz(k) = stz(k) / real(cnt(k))
          kts = k
        endif
      enddo
      write(83) iyear0,idate0,csign0,yp0,xp0,kts, &
                              &    (st(k),stz(k),sz(k),k=1,kts)
!
      close(82)
      close(80)
      rewind(83)
!
! count profiles in modified file (fort.83)
!
      nprf = 0
      do while (.true.)
        read(83,end=5,err=130) iyear,idate,csign,yp,xp,np
        nprf = nprf + 1
      end do
    5 continue
      rewind (83)
!
      allocate(nprfl(nprf))
      allocate(ndxd(nprf))
      allocate(ndxp(nprf))
!
! open direct access scratch file
!
      nb = 4*(3*kax + 5) + 8
      open (80, status='scratch', form='unformatted', &
                                     & access='direct',recl=nb)
!
! begin loop on profile file
!
      do n=1,nprf
        read (83, err=130) iyear,idate,csign,yp,xp,kt, &
                                    & (tt(k),ttz(k),zt(k),k=1,kt)
        ndxd(n) = (iyear-ayear)*100000000 + idate
!       ndxd(n) = iyear*100000000 + idate
        nprfl(n) = n
        do k=kt+1,kax
          tt(k) = mVal
          ttz(k) = mVal
        enddo
        write (80,rec=n) iyear,idate,csign,yp,xp,kt,tt,ttz,zt
      end do
!
      close (83)
!
! sort on ndxd
!
      ns = 1
      ne = nprf
      call qSort(ndxd,nprfl,ns,ne)
!
! rewrite fort.83 profile file in time sequence
!
      open (83, status='scratch', form='unformatted', &
                                      &  access='sequential')
!
! begin loop on sorted list
!
      do n=1,nprf
        read (80,rec=nprfl(n)) iyear,idate,csign,yp,xp,kt,tt,ttz,zt
        write (83, err=140) iyear,idate,csign,yp,xp,kt,& 
                                      & (tt(k),ttz(k),zt(k),k=1,kt)
      end do
!
      close(80)
      rewind(83)
!
! open assimilation profile file
!
      open (51,form='unformatted',access='sequential',err=510)
!
! open direct access scratch file and a sequential scratch file
!
      nb = 4*(3*kax + 5) + 8
      open (80,status='scratch',form='unformatted', &
                               & access='direct',recl=nb)
!
      open (81,status='scratch',form='unformatted', &
                               & access='sequential')
!
! get reference dates for first week
!
      xi = -1.0
      yj = -1.0
      jday0 = -1
      nout = -1
      okay = .false.
      do while (.not. inbnds(xi,yj) .or. jday0 .lt. jdayS)
        read (83, end=10, err=130) iyear,idate,csign,yp,xp,kt, &
                                        & (tt(k),ttz(k),zt(k),k=1,kt)

        write (cdate,'(i8.8)') idate
        read (cdate,'(4i2)') imon, iday, ihr, imin
        jday0 = jlnDayR(iyear, imon, iday, iyearR)
        xi = real_indx(xp,xt,imx)
        yj = real_indx(yp,yt,jmx)
        nout = nout + 1
      end do
      okay = .true.
   10 continue
      if (.not. okay) then
        write(6,'(a)') 'Reading modified profile file (unit 83) :'
        write(6,'(a)') ' all profiles appear out of bounds.'
        go to 130
      end if
!
      write (cdate,'(i8.8)') idate
      read (cdate,'(4i2)') imon, iday, ihr, imin
      dow0 = weekDay(iyear, imon, iday)
!
      jday0 = jlnDayR(iyear, imon, iday, iyearR)
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
      do k=kt+1,kax
        tt(k) = mVal
        ttz(k) = mVal
      enddo
      write (80,rec=nrc) iyear,idate,csign,xi,yj,kt,tt,ttz,zt
!
! continue loop on profile file
!
      do while (.true.)
        read (83, end=100, err=130) iyear,idate,csign,yp,xp,kt, &
                                         & (tt(k),ttz(k),zt(k),k=1,kt)
        write (cdate,'(i8.8)') idate
        read (cdate,'(4i2)') imon, iday, ihr, imin
        dow = weekDay(iyear, imon, iday)
        jday = jlnDayR(iyear, imon, iday, iyearR)
!
        if (jday .gt. jdayE) then
          exit
        end if
!
        if (jday .lt. jday0) then
          write(6,'(a)') 'Reading modified profile file (unit 83) :'
          write(6,'(a)') ' time sequence error.'
          go to 130
        end if
!
        if (jday .le. jday2) then
          xi = real_indx(xp,xt,imx)
          yj = real_indx(yp,yt,jmx)
          if (inbnds(xi,yj)) then
            nrc = nrc + 1
            ndxp(nrc) = int(yj*1000.0)
            ndxp(nrc) = ndxp(nrc)*1000000 + int(xi*1000.0)
            nprfl(nrc) = nrc
            do k=kt+1,kax
              tt(k) = mVal
              ttz(k) = mVal
            enddo
            write (80,rec=nrc) iyear,idate,csign,xi,yj,kt,tt,ttz,zt
          else
            nout = nout + 1
          end if
        else
!
! save data
!
          iyear0 = iyear
          idate0 = idate
          csign0 = csign
          kts = kt
          do k=1,kt
            st(k) = tt(k)
            stz(k) = ttz(k)
            sz(k) = zt(k)
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
            read (80,rec=nprfl(n)) iyear,idate,csign,xi,yj,kt,tt,ttz,zt
            write (cdate,'(i8.8)') idate
            read (cdate,'(4i2)') imon, iday, ihr, imin
            dow = weekDay(iyear, imon, iday)
            fdow = float(dow-1) + (float(ihr) + float(imin)/60.0)/24.0
            rday = fdow - 3.5
            do k=1,kt
              zk = float(k)
              se = se0 + seF * ttz(k)
              rerr = 1.0 / (se*se)
              write (81) csign,rday,tt(k),xi,yj,zk,rerr
              nobs = nobs + 1
            end do
          end do
!
          rewind (81)
!
          call calDayR(iyr,imo,idy,iyearR,jdayR)
          write (stamp,'(a,i2,a,i2,a,i4,a)') &
                 & 'm/d/y=',imo,'/',idy,'/',iyr,', h:m:s=12: 0: 0'
          write(6,'(a,i8)') stamp,nobs
          write(51,err=520) stamp,nobs
          do n=1,nobs
            read(81) csign,rday,ttk,xi,yj,zk,rerr
            write(51,err=520) csign,rday,ttk,xi,yj,zk,rerr
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
          if (inbnds(xi,yj)) then
            nrc = nrc + 1
            ndxp(nrc) = int(yj*1000.0)
            ndxp(nrc) = ndxp(nrc)*1000000 + int(xi*1000.0)
            nprfl(nrc) = nrc
            write (80,rec=nrc) iyear0,idate0,csign0,xi,yj,kts,st,stz,sz
          else
            nout = nout + 1
          end if
!
        end if
!
      end do
  100 continue
!
      close (83)
!
! sort on ndxp
!
      ns = 1
      ne = nrc
      call qSort(ndxp,nprfl,ns,ne)
!
      nobs = 0
      do n=1,nrc
        read (80,rec=nprfl(n)) iyear,idate,csign,xi,yj,kt,tt,ttz,zt
        write (cdate,'(i8.8)') idate
        read (cdate,'(4i2)') imon, iday, ihr, imin
        dow = weekDay(iyear, imon, iday)
        fdow = float(dow-1) + (float(ihr) + float(imin)/60.0)/24.0
        rday = fdow - 3.5
        do k=1,kt
          zk = float(k)
          se = se0 + seF * ttz(k)
          rerr = 1.0 / (se*se)
          write (81) csign,rday,tt(k),xi,yj,zk,rerr
          nobs = nobs + 1
        end do
      end do
!
      close (80)
      rewind (81)
!
      call calDayR(iyr,imo,idy,iyearR,jdayR)
      write (stamp,'(a,i2,a,i2,a,i4,a)') &
             & 'm/d/y=',imo,'/',idy,'/',iyr,', h:m:s=12: 0: 0'
      write(6,'(a,i8)') stamp,nobs
      write(51,err=520) stamp,nobs
      do n=1,nobs
        read(81) csign,rday,ttk,xi,yj,zk,rerr
        write(51,err=520) csign,rday,ttk,xi,yj,zk,rerr
      end do
!
      close (81)
      close (51)
!
      call w3tage('GODAS_MKASMPRF')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a)') 'Error opening profile file on unit 11'
      call w3tage('GODAS_MKASMPRF')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading profile file on unit 11'
      call w3tage('GODAS_MKASMPRF')
      call errexit(12)
!     call exit(12)
!
  130 write(6,'(a)') 'Error reading modified profile file on unit 83'
      call w3tage('GODAS_MKASMPRF')
      call errexit(13)
!     call exit(13)
!
  140 write(6,'(a)') 'Error writing profile file on unit 83'
      call w3tage('GODAS_MKASMPRF')
      call errexit(14)
!     call exit(14)
!
  210 write(6,'(a)') 'Error opening grid/mask file on unit 12'
      call w3tage('GODAS_MKASMPRF')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a)') 'Error reading grid/mask file on unit 12'
      call w3tage('GODAS_MKASMPRF')
      call errexit(22)
!     call exit(22)
!
  510 write(6,'(a)') 'Error opening profile file on unit 51'
      call w3tage('GODAS_MKASMPRF')
      call errexit(51)
!     call exit(51)
!
  520 write(6,'(a)') 'Error writing profile file on unit 51'
      call w3tage('GODAS_MKASMPRF')
      call errexit(52)
!     call exit(52)
!
  610 write(6,'(a)') 'Error reading date from command line'
      call w3tage('GODAS_MKASMPRF')
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
      function inbnds (x, y)
!
      real x, y
      logical inbnds
!
      if (x .lt. 1.01 .or. y .lt. 0.0) then
        inbnds = .false.
      else
        i = int(x)
        j = int(y)
        if (mskt(i,j) .gt. 0 .or. mskt(i+1,j) .gt. 0 .or. &
            & mskt(i+1,j+1) .gt. 0 .or. mskt(i,j+1) .gt. 0) then
          inbnds = .true.
        else
          inbnds = .false.
        end if
      end if
!
      end function inbnds
!
! -------------------------------------------------------------- 
!
      subroutine intrpZ(kt, ti, tzi, zi, ki, to, tzo, zo, ko)
!
      real ti(*), tzi(*), zi(*), to(*), tzo(*), zo(*)
      integer kt, ki, ko
!
      integer k, kk
      integer km, kp, kv2, cnt
      real zr, tzmn, tzmx
    
      kk = ki
      do k=2,ki
        km = k - 1
        if (zi(km) .le. z0g .and. (zi(k) - zi(km)) .ge. dZg) then
          kk = km
          exit
        end if
      end do
      ki = kk

      if (ki .lt. 3) then
        kt = 0
        return
      end if
    
      do k=2,ki-1
        tzi(k) = (ti(k-1) - ti(k+1)) / (zi(k+1) - zi(k-1))
        if (tzi(k) .lt. 0.0) tzi(k) = 0.0
      end do
      tzi(1) = tzi(2)
      tzi(ki) = tzi(ki-1)
    
      k = 1
      kt = 1
      if (zi(1) .lt. 30.0) then
        do while (kt .lt. ko)
          if (zo(kt) .le. zi(1)) then
            to(kt) = ti(1)
            tzo(kt) = 0.0
          else
            do while (zi(k) .lt. zo(kt) .and. k .lt. ki)
              k = k + 1
            end do
            if (zo(kt) .le. zi(k)) then
              zr = (zo(kt) - zi(k-1)) / (zi(k) - zi(k-1))
              to(kt) = ti(k-1) + (ti(k) - ti(k-1)) * zr
              tzo(kt) = tzi(k-1) + (tzi(k) - tzi(k-1)) * zr
            else
              kt = kt - 1
              exit
            end if
          end if
          kt = kt + 1
        end do
      end if

      kk = kt
      do k=1,kk
        if (to(k) .le. -5.0) kt = 0
      end do
         
      if (kt .gt. 0) then
    
        if (sclDz) then
          do k=1,kt
            tzo(k) = tzo(k) / dz(k)
          end do
        end if
    
        if (srTz) then
          do k=1,kt
            if (tzo(k) > 0.0) then
              tzo(k) = sqrt(tzo(k))
            else
              tzo(k) = 0.0
            end if
          end do
        end if
    
        if (kav .gt. 1) then
          kv2 = kav / 2
          do k=1,kt
            km = k - kv2
            if (km .lt. 1) km = 1
            kp = k + kv2
            if (kp .gt. kt) kp = kt
            cnt = 0
            ti(k) = 0.0
            do kk=km,kp
              ti(k) = ti(k) + tzo(kk)
              cnt = cnt + 1
            end do
            if (cnt .gt. 0) ti(k) = ti(k) / float(cnt)
          end do
          do k=1,kt
            tzo(k) = ti(k)
          end do
        end if
    
        tzmn = tzo(1)
        tzmx = tzo(1)
        do k=1,kt
          if (tzmn .gt. tzo(k)) tzmn = tzo(k)
          if (tzmx .lt. tzo(k)) tzmx = tzo(k)
        end do
        tzmx = tzmx - tzmn
        if (tzmx .lt. tzeps) then
          do k=1,kt
            tzo(k) = 0.0
          end do
        else
          do k=1,kt
            tzo(k) = tzo(k) - tzmn
            tzo(k) = tzo(k) / tzmx
          end do
        end if
      end if
    
      end subroutine intrpZ
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
! -------------------------------------------------------------------
!
      end program mkAsmPrf
