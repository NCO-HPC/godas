!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  TRKAVECYC
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2006-08-08
!
! ABSTRACT:  Averages data in altimetry cycle files into 1-dgree, 
!   along-track bins. Applies a calibration to bring Jason data from
!   the BUFR tank into line with previously established convention
!   of referencing SSH to the 1993-1999 mean.  A middle step in 
!   preparing altimetry for use in Global Ocean Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2008-08-08  David W. Behringer
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - ALTIMETRY CYCLE FILE IN ASCII
!     UNIT 21  - ALTIMETRY CALIBRATION FILE IN IEEE
!
!   OUTPUT FILES:
!     UNIT 51  - 1-DEGREE AVERAGED ALTIMETRY CYCLE FILE IN ASCII
!     UNIT 61  - A TABLE OF START DATE, END DATE, AND NUMBER OF RECORDS
!                IN THE CYCLE FILE
!
!   WORK FILES:  NONE
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - applyCal, jlnDayR, calDayR
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - applyCal, jlnDayR, calDayR
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  21 - ERROR OPENING UNIT 21
!     COND =  22 - ERROR READING UNIT 21
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!     COND =  53 - ERROR READING UNIT 51
!     COND =  61 - ERROR OPENING UNIT 61
!     COND =  62 - ERROR WRITING UNIT 61
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
      program trkAveCyc
!
!  trkAveCyc averages to 1-dgree bins along track
!   a calibration is also applied
!
      use tmutil_mod

      integer, parameter :: nbin = 129, joff = 65, year0 = 1985
      integer :: dte, hr, mn, cycle, orbit
      integer :: odte, ohr, omn, ocycle, oorbit, orb2
      real :: ylt, oylt, nylt, xlg, oxlg, ssh, ossh, ojd, ad
      real, allocatable, dimension(:) :: sylt, sxlg, sjd, sssh, cnt
      logical :: asc, gap
      real, save :: eps = 0.5, yoff = 65.5
      integer :: year, month, day
!  calibration
      integer :: imx, jmx
      real, allocatable, dimension(:) :: xg, yg
      real, allocatable, dimension(:,:) :: cg
      integer, allocatable, dimension(:,:) :: msk
!  end calibration
!  table
      real :: dy85
      integer :: nlines
      integer :: years, months, days, yeare, monthe, daye, yday
!  end table
      character(len=50) :: cycleFile, calFile, avFile
      character(len=30) :: str
!
      call w3tagb('GODAS_TRKAVECYC',2006,0164,0164,'NP23')
!
      allocate (sylt(nbin))
      allocate (sxlg(nbin))
      allocate (sjd(nbin))
      allocate (sssh(nbin))
      allocate (cnt(nbin))
      sylt = 0.0
      sxlg = 0.0
      sjd = 0.0
      sssh = 0.0
      cnt = 0.0
!
!  open and read the gridded calibration file
!
      open(unit=21,form='UNFORMATTED', status='OLD', &
     &                             access='SEQUENTIAL',err=210)
      read(21,err=220) imx,jmx
      allocate (xg(imx))
      allocate (yg(jmx))
      allocate (cg(imx,jmx))
      allocate (msk(imx,jmx))
      read(21,err=220) xg
      read(21,err=220) yg
      read(21,err=220) cg
      read(21,err=220) msk
      close(21)
!
!  open the unaveraged cycle file for reading on unit 11
!
      open(unit=11,form='FORMATTED', status='OLD', &
     &                             access='SEQUENTIAL',err=110)
!
!  open unit 51 for writing the averaged cycle file
!
      open(unit=51,form='FORMATTED',access='SEQUENTIAL',err=510)
!
!  read the first 2 records from unit 11
!
      read(11,'(i9,2i3,2i5,2f12.5,f8.3)',err=120,end=100)  &
     &          odte, ohr, omn, ocycle, oorbit, oylt, oxlg, ossh
      ossh = 100.0*ossh
      read(11,'(i9,2i3,2i5,2f12.5,f8.3)',err=120,end=100)  &
     &          dte, hr, mn, cycle, orbit, ylt, xlg, ssh
      ssh = 100.0*ssh
!
!  determine whether ascending or descending track and compute a new
!  orbit number that is odd for ascending and even for descending
!
      if ((ylt - oylt) .gt. 0.0) then
        asc = .TRUE.
        ad = 1.0
        orb2 = 2*oorbit - 1
      else
        asc = .FALSE.
        ad = -1.0
        orb2 = 2*oorbit
      endif
!
!  calculate the averaging bin index jlt and do sums if jlt is within
!  specified limits (1 -> 129) corresponding to 64S to 64N
!
      jlt = int(oylt + yoff)
      if (jlt .ge. 1 .and. jlt .le. nbin) then
        year = odte / 10000
        month = mod(odte,10000) / 100
        day = mod(odte,100)
        ojd = float(jlnDayR(year, month, day, year0)) - 1.0
        ojd = ojd + (float(ohr) + float(mn)/60.0) / 24.0
        sylt(jlt) = sylt(jlt) + oylt
        sxlg(jlt) = sumlon(sxlg(jlt),oxlg)
        sjd(jlt) = sjd(jlt) + ojd
        sssh(jlt) = sssh(jlt) + ossh
        cnt(jlt) = cnt(jlt) + 1.0
      endif
      odte = dte
      ohr = hr
      omn = mn
      ocycle = cycle
      oorbit = orbit
      oylt = ylt
      oxlg = xlg
      ossh = ssh
      do while (.TRUE.)
!
!  continue reading records from unit 11 until file is exhausted,
!  then jump to 100
!
        read(11,'(i9,2i3,2i5,2f12.5,f8.3)',end=100,err=120)  &
     &            dte, hr, mn, cycle, orbit, ylt, xlg, ssh
        ssh = 100.0*ssh
!
!  check for a "gap".  a "gap" is defined as a passage over land, during
!  which the direction changes from ascending to descending or vice versa
!
        gap = .FALSE.
        if (abs(ylt-oylt) .gt. eps) then
          if (asc .and. (ylt - oylt) .gt. 0.0) then
            read(11,'(25x,f12.5)',err=120,end=100) nylt
            if ((nylt - ylt) .lt. 0.0) gap = .TRUE.
            backspace(11)
          else if (.not.asc .and. (ylt - oylt) .lt. 0.0) then
            read(11,'(25x,f12.5)',err=120,end=100) nylt
            if ((nylt - ylt) .gt. 0.0) gap = .TRUE.
            backspace(11)
          endif
        endif
        if ((orbit .ne. oorbit) .or.  &
     &       (asc .and. (ylt - oylt) .lt. 0.0) .or.  &
     &                  (.not.asc .and. (ylt - oylt) .gt. 0.0) .or. gap) then
!
!  if the orbit number changes, the direction changes, or there is a "gap",
!  finish summing the previous read, compute the averages and write the
!  averaged data to unit 51
!
          jlt = int(oylt + yoff)
          if (jlt .ge. 1 .and. jlt .le. nbin) then
            year = odte / 10000
            month = mod(odte,10000) / 100
            day = mod(odte,100)
            ojd = float(jlnDayR(year, month, day, year0)) - 1.0
            ojd = ojd + (float(ohr) + float(mn)/60.0) / 24.0
            sylt(jlt) = sylt(jlt) + oylt
            sxlg(jlt) = sumlon(sxlg(jlt),oxlg)
            sjd(jlt) = sjd(jlt) + ojd
            sssh(jlt) = sssh(jlt) + ossh
            cnt(jlt) = cnt(jlt) + 1.0
          endif
          do j=1,nbin
            if (cnt(j) .gt. 0.0) then
              sylt(j) = sylt(j) / cnt(j)
              sxlg(j) = sxlg(j) / cnt(j)
              if (sxlg(j) .lt. 0.0) sxlg(j) = sxlg(j) + 360.0
              sjd(j) = sjd(j) / cnt(j)
              sssh(j) = sssh(j) / cnt(j)
!
              if (y .lt. -62.0) then
                sssh(j) = 0.0
              else
                call applyCal(sylt(j), sxlg(j), sssh(j))
              endif
!
              write(51,'(f9.2,f8.2,f9.2,3f8.2,i4,i5,i4)',err=520) &
     &            sylt(j), sxlg(j), ad, sjd(j), sssh(j), cnt(j),  &
     &            ocycle, orb2, j-joff
!
              sylt(j) = 0.0
              sxlg(j) = 0.0
              sjd(j) = 0.0
              sssh(j) = 0.0
              cnt(j) = 0.0
            endif
          enddo
          odte = dte
          ohr = hr
          omn = mn
          ocycle = cycle
          oorbit = orbit
          oylt = ylt
          oxlg = xlg
          ossh = ssh
          read(11,'(i9,2i3,2i5,2f12.5,f8.3)',end=100,err=120)  &
     &              dte, hr, mn, cycle, orbit, ylt, xlg, ssh
          ssh = 100.0*ssh
          if ((ylt - oylt) .gt. 0.0) then
            asc = .TRUE.
            ad = 1.0
            orb2 = 2*orbit - 1
          else
            asc = .FALSE.
            ad = -1.0
            orb2 = 2*orbit
          endif
          jlt = int(oylt + yoff)
          if (jlt .ge. 1 .and. jlt .le. nbin) then
            year = odte / 10000
            month = mod(odte,10000) / 100
            day = mod(odte,100)
            ojd = float(jlnDayR(year, month, day, year0)) - 1.0
            ojd = ojd + (float(ohr) + float(mn)/60.0) / 24.0
            sylt(jlt) = sylt(jlt) + oylt
            sxlg(jlt) = sumlon(sxlg(jlt),oxlg)
            sjd(jlt) = sjd(jlt) + ojd
            sssh(jlt) = sssh(jlt) + ossh
            cnt(jlt) = cnt(jlt) + 1.0
          endif
        else
!
!  continue summing if there is no orbit or direction change and no "gap"
!
          jlt = int(oylt + yoff)
          if (jlt .ge. 1 .and. jlt .le. nbin) then
            year = odte / 10000
            month = mod(odte,10000) / 100
            day = mod(odte,100)
            ojd = float(jlnDayR(year, month, day, year0)) - 1.0
            ojd = ojd + (float(ohr) + float(mn)/60.0) / 24.0
            sylt(jlt) = sylt(jlt) + oylt
            sxlg(jlt) = sumlon(sxlg(jlt),oxlg)
            sjd(jlt) = sjd(jlt) + ojd
            sssh(jlt) = sssh(jlt) + ossh
            cnt(jlt) = cnt(jlt) + 1.0
          endif
        endif
        odte = dte
        ohr = hr
        omn = mn
        ocycle = cycle
        oorbit = orbit
        oylt = ylt
        oxlg = xlg
        ossh = ssh
      enddo
100   continue
!
!  finish summing the last record, compute the final averages and
!  write the averaged data to unit 51
!
      jlt = int(oylt + yoff)
      if (jlt .ge. 1 .and. jlt .le. nbin) then
        year = odte / 10000
        month = mod(odte,10000) / 100
        day = mod(odte,100)
        ojd = float(jlnDayR(year, month, day, year0)) - 1.0
        ojd = ojd + (float(ohr) + float(mn)/60.0) / 24.0
        sylt(jlt) = sylt(jlt) + oylt
        sxlg(jlt) = sumlon(sxlg(jlt),oxlg)
        sjd(jlt) = sjd(jlt) + ojd
        sssh(jlt) = sssh(jlt) + ossh
        cnt(jlt) = cnt(jlt) + 1.0
      endif
      do j=1,nbin
        if (cnt(j) .gt. 0.0) then
          sylt(j) = sylt(j) / cnt(j)
          sxlg(j) = sxlg(j) / cnt(j)
          if (sxlg(j) .lt. 0.0) sxlg(j) = sxlg(j) + 360.0
          sjd(j) = sjd(j) / cnt(j)
          sssh(j) = sssh(j) / cnt(j)
!
          if (y .lt. -62.0) then
            sssh(j) = 0.0
          else
            call applyCal(sylt(j), sxlg(j), sssh(j))
          endif
!
          write(51,'(f9.2,f8.2,f9.2,3f8.2,i4,i5,i4)',err=520) &
     &        sylt(j), sxlg(j), ad, sjd(j), sssh(j), cnt(j),  &
     &        ocycle, orb2, j-joff
        endif
      enddo
      close(11)
!
!  reread the output file and add its start date, end date, line count
!  and name to a table file
!
      rewind(51)
!
      inquire(unit=51,name=avFile)
!
      read(51,'(26x,f8.2)',err=530) dy85
      yday = int(dy85)
      call calDayR (years, months, days, year0, yday)
      nlines = 1
      do while (.TRUE.)
        read(51,'(26x,f8.2)',end=200,err=530) dy85
        nlines = nlines + 1
      enddo
200   continue
      yday = int(dy85)
      call calDayR (yeare, monthe, daye, year0, yday)
!
      close(51)
!
      open(unit=61,form='FORMATTED', status='UNKNOWN', &
     &                    access='SEQUENTIAL',position='APPEND',err=610)
      write(61,'(i4,a,i2.2,a,i2.2,2x,i4,a,i2.2,a,i2.2,2x,i5,2x,a)',err=620) &
     &  years,'-',months,'-',days, yeare,'-',monthe,'-',daye,  &
     &  nlines, trim(avFile)
      write(6,'(a,i4,a,i2.2,a,i2.2,2x,i4,a,i2.2,a,i2.2,2x,i5,2x,a)')  &
     & 'Write to 61:  ',&
     &  years,'-',months,'-',days, yeare,'-',monthe,'-',daye,  &
     &  nlines, trim(avFile)
!
      close(61)
!
      call w3tage('GODAS_TRKAVECYC')
      call errexit(0)
!     call exit(0)
!
110   write(6,'(a)') 'Error opening altimetry cycle file on unit 11'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(11)
!     call exit(11)
!
120   write(6,'(a)') 'Error reading altimetry cycle file on unit 11'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(12)
!     call exit(12)
!
210   write(6,'(a)') 'Error opening altimetry calibration file on unit 21'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(21)
!     call exit(21)
!
220   write(6,'(a)') 'Error reading altimetry calibration file on unit 21'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(22)
!     call exit(22)
!
510   write(6,'(a)') 'Error opening averaged altimetry cycle file on unit 51'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(51)
!     call exit(51)
!
520   write(6,'(a)') 'Error writing averaged altimetry cycle file on unit 51'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(52)
!     call exit(52)
!
530   write(6,'(a)') 'Error reading averaged altimetry cycle file on unit 51'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(53)
!     call exit(53)
!
610   write(6,'(a)') 'Error opening table file on unit 61'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(61)
!     call exit(61)
!
620   write(6,'(a)') 'Error writing table file on unit 61'
      call w3tage('GODAS_TRKAVECYC')
      call errexit(62)
!     call exit(62)
!
      contains
!
! --------------------------------------------------------------


      function sumlon(sumin,rlon)
      real :: sumin, rlon
! since rlon may range from -180 to 180, want to be sure avg lon will
!  be accurate as tracks cross dateline.  keep sum going in one direction
!  whether it be positive or negative

! do simple sum if lon between -90 and 90.  (to ensure points near greenwich are left as is).
      if(rlon.gt.-90.and.rlon.lt.90) then
!   we are assuming there would not be more than a 90 degree jump in longitude
        sumlon=sumin+rlon                
! for points closer to dateline...
      else if(rlon.gt.0)then     ! between 90 and 180
        if(sumin.ge.0)then         ! both pos 
          sumlon=sumin+rlon       ! simple sum
!         print*,'sumlon=sumin+rlon       ! simple sum'
        else                      !incoming sum neg and new lon pos
          sumlon=sumin+(rlon-360) !make sum increment in neg direction
!         print*,'sumlon=sumin+(rlon-360) !make sum increment in neg direction'
        endif
      else                       ! between -90 and -180
        if(sumin.le.0)then         ! both neg (or sumin 0) 
          sumlon=sumin+rlon       ! simple sum
!         print*,'sumlon=sumin+rlon       ! simple sum'
        else                      !incoming sum pos and new lon neg
          sumlon=sumin+(rlon+360) !make sum continue in pos direction
!         print*,'sumlon=sumin+(rlon+360) !make sum continue in pos direction'
        endif
      endif

      end function sumlon

! --------------------------------------------------------------
!
      subroutine applyCal(lat, lon, hgt)
!
      real :: lat, lon, hgt
!
      if (lat .le. yg(1)) then
        jm = 1
        jp = 2
      else if (lat .ge. yg(jmx)) then
        jm = jmx-1
        jp = jmx
      else
        do jl=2,jmx
          if (lat .ge. yg(jl-1) .and. lat .le. yg(jl)) then
            jm = jl-1
            jp = jl
            exit
          endif
        enddo
      endif
!
      if (lon .le. xg(1)) then
        im = imx
        ip = 1
      else if (lon .ge. xg(imx)) then
        im = imx
        ip = 1
      else
        do il=2,imx
          if (lon .ge. xg(il-1) .and. lon .le. xg(il)) then
            im = il-1
            ip = il
            exit
          endif
        enddo
      endif
!
      n = msk(im,jm) + msk(im,jp) + msk(ip,jm) + msk(ip,jp)
      if (n .eq. 4) then
        dxdy = (xg(ip) - xg(im)) * (yg(jp) - yg(jm))
        dxm = lon - xg(im)
        dxp = xg(ip) - lon
        dym = lat - yg(jm)
        dyp = yg(jp) - lat
        rmm = dxp*dyp/dxdy
        rmp = dxp*dym/dxdy
        rpm = dxm*dyp/dxdy
        rpp = dxm*dym/dxdy
        c = cg(im,jm)*rmm + cg(im,jp)*rmp + cg(ip,jm)*rpm + cg(ip,jp)*rpp
      else if (n .eq. 0) then
        c = 0.0
      else
        n = 0
        c = 0.0
        if (msk(im,jm) .ne. 0) then
          c = c + cg(im,jm)
          n = n + 1
        endif
        if (msk(im,jp) .ne. 0) then
          c = c + cg(im,jp)
          n = n + 1
        endif
        if (msk(ip,jm) .ne. 0) then
          c = c + cg(ip,jm)
          n = n + 1
        endif
        if (msk(ip,jp) .ne. 0) then
          c = c + cg(ip,jp)
          n = n + 1
        endif
        c = c / float(n)
      endif
!
      hgt = hgt - c
!
      end subroutine applyCal
!
      end program trkAveCyc
