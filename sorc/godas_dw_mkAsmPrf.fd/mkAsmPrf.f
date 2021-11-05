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
! 2007-08-28  David W. Behringer - Extend profiles from 30 levels max
!              to 35 levels, using climatological temperature data.
! 2015-09-29  David W. Behringer - If a platform reports mutltiple
!              profiles on a single day, they are averaged
! 2017-09-15  David W. Behringer - Allows daily averages of massive clusters
!              of profiles, e.g. those air-dropped into hurricanes
! 2018-01-15  David W. Behringer - Fix bug causing very occasional misplacement
!              of profiles. Fix bug causing occasional zero-filled profiles
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - TEMPERATURE PROFILE DATA IN IEEE
!     UNIT 12  - GRID/MASK FOR OCEAN MODEL
!     UNIT 13  - TEMPERATURE CLIMATOLOGY DATA IN IEEE
!
!   OUTPUT FILES:
!     UNIT 51  - TEMPERATURE PROFILE DATA FOR ASSIMILATION IN IEEE
!     UNIT 61  - EXTENDED TEMPERATURE PROFILE DATA (WILL ALSO BE USED
!                 FOR COMPUTING SYNTHETIC SALINITY) IN IEEE
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
!                  jlnDayR, jlnDayC, calDayR, calDayC, weekDay,
!                  intZclmT, getZExt, cmpTz, tBlend
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - real_indx, inbnds, intrpZ, qSort, swap, ssort
!                  jlnDayR, calDayR, weekDay, intZclmT,
!                  getZExt, cmpTz, tBlend
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  13 - ERROR READING UNIT 83
!     COND =  14 - ERROR WRITING UNIT 83
!     COND =  21 - ERROR OPENING UNIT 12
!     COND =  22 - ERROR READING UNIT 12
!     COND =  31 - ERROR OPENING UNIT 13
!     COND =  32 - ERROR READING UNIT 13
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!     COND =  61 - ERROR OPENING UNIT 61
!     COND =  62 - ERROR WRITING UNIT 61
!     COND =  91 - ERROR READING FROM COMMAND LINE
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
      integer imx, jmx, kmx, kax
      integer imc, jmc, kmc
      character csign*8,sid*2,dtyp*2,qkey*1
      character csign0*8,sid0*2,dtyp0*2,qkey0*1
      character str*80, cdate*8
      real, allocatable, dimension(:) :: pt, pz, w
      real, allocatable, dimension(:) :: xt, yt, zt, zw, dz, tt, ttz
      real, allocatable, dimension(:) :: sz, st, stz
      integer, allocatable, dimension(:) :: cnt
      integer, allocatable, dimension(:,:) :: mskt
      real, allocatable, dimension(:,:) :: b
      real, allocatable, dimension(:) :: xc, yc, zc, te
      real, allocatable, dimension(:,:,:) :: tr, tc
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
!     kan = minimum number of model levels for assimilation
!     kax = maximum number of model levels for assimilation
!     nwks = time window (weeks) for profile data
!
      kan = 30
      kax = 35
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
!  to set std errors between 1.5 and 2.5
!    SE0 = 1.5 and SEF = 1.0
!  then 1/e.v. = 1/(se*se)
!  if std error = 0.5, then EVR = 4.0
!                 1.0             1.0
!                 1.5             0.4444
!
      se0 = 1.5
      seF = 1.0
! set missing value
      mVal = 999.999
!
! get date from command line
!
      narg = iargc()
      if (narg .eq. 0) go to 910
      call getarg(1,str)
      read(str,'(i4,2i2)',err=910) ayear, amonth, aday
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
      allocate(w(kmx))
      allocate(st(kmx))
      allocate(stz(kmx))
      allocate(sz(kmx))
!
! read 3d array of climatological temperature
!
      open(unit=13,status='old',form='UNFORMATTED', &
                            & access='sequential', err=310)
!
      read(13,err=320) str
      read(13,err=320) stamp,dpm,imc,jmc,kmc
!
      allocate(xc(imc))
      allocate(yc(jmc))
      allocate(zc(kmc))
      allocate(tr(imc,jmc,kmc))
      allocate(tc(imc,jmc,kmx))
      allocate(b(imc,jmc))
      allocate(te(kmx))
!
      rewind(13)
!
      do k=1,kmc
        read(13,err=320) str
        read(13,err=320) stamp,dpm,i,j,km,k0,month,xc,yc,zc,b
!
        do j=1,jmc
          do i=1,imc
            tr(i,j,k) = b(i,j)
          end do
        end do
!
      end do
!
      close(13)
!
      do k=1,kmc
        zc(k) = 0.01*zc(k)
      enddo
!
! interpolate temperature climatology to model vertical coordinate
!
      call intZclmT
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
! open file for Z-interpolated and Z-extended temperature profiles
! it will then be used to write an assimilation file for temperature
! and by a subsequent program to write a synthetic salinity file
!
      open (61,form='unformatted',access='sequential',err=610)
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
      allocate(cnt(npmx))
!
! open direct access scratch file and a sequential scratch file
!
      nb = 4*(2*kax + 5) + 13
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
          call intrpZ(kt,pt,pz,np,tt,zt,kmx)
        end if                                ! bug fix
        if (kt .gt. 3) then                   ! 3/15/2018
          if (kt .gt. kax) kt = kax
          do k=kt+1,kax
            tt(k) = mVal
          enddo
          nrc = nrc+1
          write(80,rec=nrc) iyear,idate,csign,sid,dtyp, &
                           &  qkey,yp,xp,kt,(tt(k),zt(k),k=1,kax)
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
      read(80,rec=n) iyear,idate,csign,sid, &
                      &  dtyp,qkey,yp,xp,kt,(tt(k),zt(k),k=1,kt)
      csign0 = csign
      sid0 = sid
      dtyp0 = dtyp
      qkey = qkey0
      idate0 = idate
      iyear0 = iyear
      imon0 = imon
      iday0 = iday
      xp0 = xp
      yp0 = yp
      do k=1,kt
        cnt(k) = 1
        st(k) = tt(k)
        sz(k) = zt(k)
      enddo
      do k=kt+1,kax
        cnt(k) = 0
        st(k) = 0.0
        sz(k) = 0.0
      enddo
      do i=2,nrc
        read(82,'(a8,i4,2i2,i10)') csign,iyear,imon,iday,n
        read(80,rec=n) iyear,idate,csign,sid, &
                         &  dtyp,qkey,yp,xp,kt,(tt(k),zt(k),k=1,kt)
        if (csign .ne. csign0 .or. iyear .ne. iyear0  .or. &
                     & imon .ne. imon0  .or. iday .ne. iday0) then
          do k=1,kax
            if (cnt(k) .gt. 0) then
              st(k) = st(k) / real(cnt(k))
              kts = k
            endif
          enddo
          write(83) iyear0,idate0,csign0,sid0, &
                       &  dtyp0,qkey0,yp0,xp0,kts,(st(k),sz(k),k=1,kts)
          csign0 = csign
          sid0 = sid
          dtyp0 = dtyp
          qkey = qkey0
          iyear0 = iyear
          idate0 = idate
          imon0 = imon
          iday0 = iday
          xp0 = xp
          yp0 = yp
          kts = kt
          do k=1,kt
            cnt(k) = 1
            st(k) = tt(k)
            sz(k) = zt(k)
          enddo
          do k=kt+1,kax
            cnt(k) = 0
            st(k) = 0.0
            sz(k) = 0.0
          enddo
        endif
        do k=1,kt
          cnt(k) = cnt(k) + 1
          st(k) = st(k) + tt(k)
          sz(k) = zt(k)
        enddo
      enddo
      do k=1,kax
        if (cnt(k) .gt. 0) then
          st(k) = st(k) / real(cnt(k))
          kts = k
        endif
      enddo
      write(83) iyear0,idate0,csign0,sid0, &
                        &  dtyp0,qkey0,yp0,xp0,kts,(st(k),sz(k),k=1,kts)
!
      close(82)
      close(80)
      rewind(83)
!
! count profiles in modified file (fort.83)
!
      nprf = 0
      do while (.true.)
        read(83,end=6,err=130) iyear,idate,csign,sid,dtyp,qkey,yp,xp,np
        nprf = nprf + 1
      end do
    6 continue
      rewind (83)
!
      allocate(nprfl(nprf))
      allocate(ndxd(nprf))
      allocate(ndxp(nprf))
!
! open direct access scratch file
!
      nb = 4*(2*kax + 5) + 13
      open (80, status='scratch', form='unformatted', &
                                     & access='direct',recl=nb)
!
! begin loop on profile file
!
      do n=1,nprf
        read (83, err=130) iyear,idate,csign,sid, &
                           &  dtyp,qkey,yp,xp,kt,(tt(k),zt(k),k=1,kt)
        ndxd(n) = (iyear-ayear)*100000000 + idate
!       ndxd(n) = iyear*1000000 + idate/100  - bug fix 01/2018
        nprfl(n) = n
        do k=kt+1,kax
          tt(k) = mVal
        enddo
        write (80,rec=n) iyear,idate,csign,sid, &
                                       &  dtyp,qkey,yp,xp,kt,tt,zt
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
                                           & access='sequential')
!
! begin loop on sorted list
!
      do n=1,nprf
        read (80,rec=nprfl(n)) iyear,idate,csign,sid, &
                                       &  dtyp,qkey,yp,xp,kt,tt,zt
        write (83, err=140) iyear,idate,csign,sid, &
                              &  dtyp,qkey,yp,xp,kt,(tt(k),zt(k),k=1,kt)
      end do
!
      close(80)
      rewind(83)
!
! read and extend temperature profiles, keeping only
! those in bounds and within desired dates
!
! open assimilation profile file
!
      open (51,form='unformatted',access='sequential',err=510)
!
      nprf = 0
      do while (.true.)
        read (83, end=10, err=130) iyear,idate,csign,sid, &
                              &  dtyp,qkey,yp,xp,kt,(tt(k),zt(k),k=1,kt)
        write (cdate,'(i8.8)') idate
        read (cdate,'(4i2)') imon, iday, ihr, imin
        jday = jlnDayR(iyear, imon, iday, iyearR)
        xi = real_indx(xp,xt,imx)
        yj = real_indx(yp,yt,jmx)
        if (jday .ge. jdayS .and. jday .le. jdayE .and. inbnds(xi,yj)) then
          i = xi
          j = yj
          km = min(mskt(i,j),mskt(i,j+1),mskt(i+1,j+1),mskt(i+1,j))
          if (kt .gt. 0 .and. km .ge. kan) then
            if (kt .lt. km) then
              kex = getZExt(xp, yp)
              if (kex .gt. 0) then
                if (te(kt+1) .le. tt(kt)) then
                  do k=kt+1,kex
                    tt(k) = te(k)
                  enddo
                else
                  call tBlend(kt,kex)
                endif
                kt = kex
              endif
            endif
            if (kt .gt. 0) then
              write(61,err=620) iyear,idate,csign,sid, &
                             &  dtyp,qkey,yp,xp,kt,(zt(k),tt(k),k=1,kt)
              nprf = nprf + 1
            endif
          endif
        else if (jday .gt. jdayE) then
          exit
        endif
      enddo
   10 continue
!
      close (83)
!
      if (nprf .eq. 0) then
        write(6,'(a)') 'Reading modified profile file (unit 83) :'
        write(6,'(a)') ' all profiles appear out of bounds.'
        go to 130
      end if
!
      rewind (61)
!
! open direct access scratch file and a sequential scratch file
!
      nb = 4*(2*kax + 7) + 13
      open (80,status='scratch',form='unformatted', &
                               & access='direct',recl=nb)
!
      open (81,status='scratch',form='unformatted', &
                               & access='sequential')
!
! get reference dates for first week
!
      read (61, err=630) iyear,idate,csign,sid, &
                            &  dtyp,qkey,yp,xp,kt,(zt(k),tt(k),k=1,kt)
      write (cdate,'(i8.8)') idate
      read (cdate,'(4i2)') imon, iday, ihr, imin
      jday0 = jlnDayR(iyear, imon, iday, iyearR)
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
      write (80,rec=nrc) iyear,idate,csign,sid, &
                       &  dtyp,qkey,xp,yp,xi,yj,kt,(tt(k),zt(k),k=1,kt)
!
! loop on profile file
!
      do while (.true.)
        read (61, end=100, err=630) iyear,idate,csign,sid, &
                            &  dtyp,qkey,yp,xp,kt,(zt(k),tt(k),k=1,kt)
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
          write(6,'(a)') 'Reading profile file (unit 61) :'
          write(6,'(a)') ' time sequence error.'
          go to 630
        end if
!
        if (jday .le. jday2) then
          xi = real_indx(xp,xt,imx)
          yj = real_indx(yp,yt,jmx)
          nrc = nrc + 1
          ndxp(nrc) = int(yj*1000.0)
          ndxp(nrc) = ndxp(nrc)*1000000 + int(xi*1000.0)
          nprfl(nrc) = nrc
          write (80,rec=nrc) iyear,idate,csign,sid, &
                       &  dtyp,qkey,xp,yp,xi,yj,kt,(tt(k),zt(k),k=1,kt)
        else
!
! save data
!
          iyear0 = iyear
          idate0 = idate
          csign0 = csign
          sid0 = sid
          dtyp0 = dtyp
          qkey = qkey0
          xp0 = xp
          yp0 = yp
          kts = kt
          do k=1,kt
            st(k) = tt(k)
            sz(k) = zt(k)
          end do
!
! sort on ndxp
!
          ns = 1
          ne = nrc
          call qSort(ndxp,nprfl,ns,ne)
!
! read data in sorted sequence, compute error estimate and write
! model assimilation file
!
          nobs = 0
          do n=1,nrc
            read (80,rec=nprfl(n)) iyear,idate,csign,sid, &
                      & dtyp,qkey,xp,yp,xi,yj,kt,(tt(k),zt(k),k=1,kt)
            write (cdate,'(i8.8)') idate
            read (cdate,'(4i2)') imon, iday, ihr, imin
            dow = weekDay(iyear, imon, iday)
            fdow = float(dow-1) + (float(ihr) + float(imin)/60.0)/24.0
            rday = fdow - 3.5
!
            if (kt .gt. 2) then
              call cmpTz(kt,tt,zt,ttz,w)
            else
              ttz(1) = 0.0
              ttz(2) = 0.0
            endif
!
            i = xi
            j = yj
            km = min(mskt(i,j),mskt(i,j+1),mskt(i+1,j+1),mskt(i+1,j))
            kt = min(kt,km,kax)
!
            do k=1,kt
              zk = float(k)
              se = se0 + seF * ttz(k)
              rerr = 1.0 / (se*se)
              write (81) csign,rday,tt(k),xi,yj,zk,rerr
              nobs = nobs + 1
            end do
!
          end do
!
          rewind (81)
!
          call calDayR(iyr,imo,idy,iyearR,jdayR)
!         call calDayR(iyr,imo,idy,iyear0,jdayR) - bug fix 01/2018
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
          xi = real_indx(xp0,xt,imx)
          yj = real_indx(yp0,yt,jmx)
!         xi = real_indx(xpS,xt,imx) - bug fix 01/2018
!         yj = real_indx(ypS,yt,jmx)
          nrc = nrc + 1
          ndxp(nrc) = int(yj*1000.0)
          ndxp(nrc) = ndxp(nrc)*1000000 + int(xi*1000.0)
          nprfl(nrc) = nrc
          write (80,rec=nrc) iyear0,idate0,csign0,sid0, &
                  &  dtyp0,qkey0,xp0,yp0,xi,yj,kts,(st(k),sz(k),k=1,kts)
!
        end if
!
      end do
  100 continue
!
      close (61)
!
! sort on ndxp
!
      ns = 1
      ne = nrc
      call qSort(ndxp,nprfl,ns,ne)
!
! read data in sorted sequence, compute error estimate and write
! model assimilation file
!
      nobs = 0
      do n=1,nrc
        read (80,rec=nprfl(n)) iyear,idate,csign,sid, &
                      &  dtyp,qkey,xp,yp,xi,yj,kt,(tt(k),zt(k),k=1,kt)
        write (cdate,'(i8.8)') idate
        read (cdate,'(4i2)') imon, iday, ihr, imin
        dow = weekDay(iyear, imon, iday)
        fdow = float(dow-1) + (float(ihr) + float(imin)/60.0)/24.0
        rday = fdow - 3.5
!
        if (kt .gt. 2) then
          call cmpTz(kt,tt,zt,ttz,w)
        else
          ttz(1) = 0.0
          ttz(2) = 0.0
        endif
!
        i = xi
        j = yj
        km = min(mskt(i,j),mskt(i,j+1),mskt(i+1,j+1),mskt(i+1,j))
        kt = min(kt,km,kax)
!
        do k=1,kt
          zk = float(k)
          se = se0 + seF * ttz(k)
          rerr = 1.0 / (se*se)
          write (81) csign,rday,tt(k),xi,yj,zk,rerr
          nobs = nobs + 1
        end do
!
      end do
!
      close (80)
      rewind (81)
!
      call calDayR(iyr,imo,idy,iyearR,jdayR)
!     call calDayR(iyr,imo,idy,iyear0,jdayR) - bug fix 01/2018
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
  310 write(6,'(a)') 'Error opening temperature file on unit 13'
      call w3tage('GODAS_MKASMPRF')
      call errexit(31)
!     call exit(31)
!
  320 write(6,'(a)') 'Error reading temperature file on unit 13'
      call w3tage('GODAS_MKASMPRF')
      call errexit(32)
!     call exit(32)
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
  610 write(6,'(a)') 'Error opening profile file on unit 61'
      call w3tage('GODAS_MKASMPRF')
      call errexit(61)
!     call exit(61)
!
  620 write(6,'(a)') 'Error writing profile file on unit 61'
      call w3tage('GODAS_MKASMPRF')
      call errexit(62)
!     call exit(62)
!
  630 write(6,'(a)') 'Error reading profile file on unit 61'
      call w3tage('GODAS_MKASMPRF')
      call errexit(63)
!     call exit(63)
!
  910 write(6,'(a)') 'Error reading date from command line'
      call w3tage('GODAS_MKASMPRF')
      call errexit(81)
!     call exit(81)
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
      subroutine intrpZ(kt, ti, zi, ki, to, zo, ko)
!
      real ti(*), zi(*), to(*), zo(*)
      integer kt, ki, ko
!
      integer k, kk
      integer km, kp, kv2, cnt
      real zr

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
      endif

      kt = 0
      if (zi(1) .lt. 30.0) then
        k = 1
        kt = 1
        do while (kt .lt. ko)
          if (zo(kt) .le. zi(1)) then
            to(kt) = ti(1)
          else
            do while (zi(k) .lt. zo(kt) .and. k .lt. ki)
              k = k + 1
            end do
            if (zo(kt) .le. zi(k)) then
              zr = (zo(kt) - zi(k-1)) / (zi(k) - zi(k-1))
              to(kt) = ti(k-1) + (ti(k) - ti(k-1)) * zr
            else
              kt = kt - 1
              exit
            end if
          end if
          kt = kt + 1
        end do
!
        kk = kt
        do k=1,kk
          if (to(k) .le. -5.0) kt = 0
        end do
!
      end if
!
      if (kt .gt. kax) kt = kax
!
      end subroutine intrpZ
!
! --------------------------------------------------------------
!
      subroutine cmpTz(kt, t, z, tz, tw)
!
      integer kt, cnt
      real t(*), z(*), tz(*), tw(*)
!
      do k=2,kt-1
        tz(k) = (t(k-1) - t(k+1)) / (z(k+1) - z(k-1))
        if (tz(k) .lt. 0.0) tz(k) = 0.0
      end do
      tz(1) = tz(2)
      tz(kt) = tz(kt-1)

      if (sclDz) then
        do k=1,kt
          tz(k) = tz(k) / dz(k)
        end do
      end if
!
      if (srTz) then
        do k=1,kt
          if (tz(k) > 0.0) then
            tz(k) = sqrt(tz(k))
          else
            tz(k) = 0.0
          end if
        end do
      end if
!
      if (kav .gt. 1) then
        kv2 = kav / 2
        do k=1,kt
          km = k - kv2
          if (km .lt. 1) km = 1
          kp = k + kv2
          if (kp .gt. kt) kp = kt
          cnt = 0
          tw(k) = 0.0
          do kk=km,kp
            tw(k) = tw(k) + tz(kk)
            cnt = cnt + 1
          end do
          if (cnt .gt. 0) tw(k) = tw(k) / float(cnt)
        end do
        do k=1,kt
          tz(k) = tw(k)
        end do
      end if
!
      tzmn = tz(1)
      tzmx = tz(1)
      do k=1,kt
        if (tzmn .gt. tz(k)) tzmn = tz(k)
        if (tzmx .lt. tz(k)) tzmx = tz(k)
      end do
      tzmx = tzmx - tzmn
      if (tzmx .lt. tzeps) then
        do k=1,kt
          tz(k) = 0.0
        end do
      else
        do k=1,kt
          tz(k) = tz(k) - tzmn
          tz(k) = tz(k) / tzmx
        end do
      end if
!
      end subroutine cmpTz
!
! --------------------------------------------------------------
!
      subroutine intZclmT
!
      integer k, km, kim, kip
      real rzm, rzp, zr
!
      do km=1,kmx
        if (zt(km) .le. zc(1)) then
          kim = 0
          rzm = 0.0
          kip = 0
          rzp = 1.0
        else if (zt(km) .ge. zc(kmc)) then
          kim = kmc
          rzm = 1.0
          kip = kmc
          rzp = 0.0
        else
          do k=2,kmc
            if (zt(km) .ge. zc(k-1) .and. zt(km) .le. zc(k)) then
              kim = k-1
              kip = k
              zr = zc(k) - zc(k-1)
              rzm = (zc(k) - zt(km)) / zr
              rzp = (zt(km) - zc(k-1)) / zr
              exit
            endif
          enddo
        endif
!
        tc(:,:,km) = rzm * tr(:,:,kim) + rzp * tr(:,:,kip)
!
      enddo
!
      end subroutine intZclmT
!
! --------------------------------------------------------------
!
      function getZExt(x,y) result(ke)
!
      real x, y
      integer ke
      integer i, j, k, km, kp
      real dxm, dxp, dym, dyp
!
      ke = -1
!
!  Check for regions where no extension will be made
!
!  Black Sea
      if (x .ge. 27.0 .and. x .le. 42.0 .and. &
                     & y .ge. 41.0 .and. y .le. 48.0) then
        ke = 0
      endif
!  Philippine Sea
      if (x .ge. 116.0 .and. x .le. 123.0 .and. &
                      & y .ge. 4.0 .and. y .le. 10.0) then
        ke = 0
      endif
!  Gibraltar
      if (x .ge. 352.5 .and. x .le. 356.5 .and. &
                      & y .ge. 35.0 .and. y .le. 37.0) then
        ke = 0
      endif
!
!  Check for regions where extension must be set to a constant
!
!  Red Sea
      if (x .ge. 33.0 .and. x .le. 44.0 .and. &
                     & y .ge. 13.0 .and. y .le. 29.0) then
        te = 21.55
        ke = kmx
      endif
!
!  Otherwise interpolate a climatological profile
!
      if (ke .lt. 0) then
        i = -1
        if (x .ge. xc(1) .and. x .le. xc(imc)) then
          i = 1
          do while (xc(i+1) .lt. x)
            i = i + 1
          enddo
        endif
        j = -1
        if (y .ge. yc(1) .and. y .le. yc(jmc)) then
          j = 1
          do while (yc(j+1) .lt. y)
            j = j + 1
          enddo
        endif
!
        if (i .gt. 0 .and. j .gt. 0) then
          dxm = (x - xc(i)) / (xc(i+1) - xc(i))
          dxp = (xc(i+1) - x) / (xc(i+1) - xc(i))
          dym = (y - yc(j)) / (yc(j+1) - yc(j))
          dyp = (yc(j+1) - y) / (yc(j+1) - yc(j))
          do k=1,kmx
            te(k) = tc(i,j,k)*dxp*dyp + tc(i+1,j,k)*dxm*dyp +  &
                  & tc(i,j+1,k)*dxp*dym + tc(i+1,j+1,k)*dxm*dym
          enddo
          ke = kmx
        else
          ke = 0
        endif
      endif
!
      end function getZExt
!
! --------------------------------------------------------------
!
      subroutine tBlend(kt,kex)
!
      integer kt, kex
!
      do k=kt+1,kex
        if (te(k) .lt. tt(kt)) then
          kb = k
          exit
        endif
      enddo
!
      if (kb .gt. 0 .and. kb .le. kex) then
        zr = zt(kb) - zt(kt)
        do k=kt+1,kb-1
          dzt = (zt(kb)-zt(k))/zr
          dzb = (zt(k)-zt(kt))/zr
          tt(k) = tt(kt)*dzt + te(kb)*dzb
        enddo
        do k=kb,kex
          tt(k) = te(k)
        enddo
      else
        kex = kt
      endif
!
      end subroutine tBlend
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
      end program mkAsmPrf
