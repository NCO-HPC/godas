! MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  NC2GRIB
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2005-08-15
!
! ABSTRACT:  Reads a GODAS time_mean netCDF file and writes fields to a
!             GRIB file
!
! PROGRAM HISTORY LOG:
! 2005-08-15  David W. Behringer
! 2006-02-13  David W. Behringer - Modifications to extILDNc (surface
!             isothermal layer) and extMLDNc (surface mixed layer) to mask
!             out ambiguous regions.  These are essentially at high latitudes
!             where there may be inversions or ice cover.
!
! USAGE:
!   INPUT FILES:
!     INPUT TIME_MEAN FILE IS MANAGED BY NETCDF LIBRARY
!
!   OUTPUT FILES:
!     OUTPUT GRIB FILE IS MANAGED BY NCEP W3 LIBRARY
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     NONE
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - getDatePrd, ext3dFNc, ext2dFNc, extILDNc, extMLDNc,
!                  jt_interp, ju_interp, press, zeta, theta,
!                  atg, density, lnstr, leapYear
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit, putgb
!       BACIO    - baopenw, baclose
!       NETCDF   - nf_open, nf_inq, nf_inq_dim, nf_inq_var,
!                  nf_get_vara_real, nf_get_var_real,
!                  nf_get_att_real, nf_close
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - getDatePrd, ext3dFNc, ext2dFNc, extILDNc, extMLDNc,
!                  jt_interp, ju_interp, lnstr
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit, putgb
!       BACIO    - baopenw, baclose
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  21 - ERROR OPENING NETCDF TIME_MEAN FILE
!     COND =  22 - ERROR READING NETCDF TIME_MEAN FILE
!     COND =  23 - ERROR CLOSING NETCDF TIME_MEAN FILE
!     COND =  24 - ERROR GETTING DATE FROM NETCDF TIME_MEAN FILE
!     COND =  31 - ERROR OPENING OUTPUT GRIB FILE
!     COND =  32 - ERROR WRITING OUTPUT GRIB FILE
!     COND =  33 - ERROR CLOSING OUTPUT GRIB FILE
!     COND =  41 - ERROR READING FROM COMMAND LINE
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     Reads an input netCDF filename from the command line.
!     The input filename is assumed to have the form 
!      time_mean.00yyyy.mm.dd.X.nc
!     The output file will have the form time_mean.00yyyy.mm.dd.X.grb
!     The time stamp on the GRIB file will be one day earlier than the
!      stamp on the netCDF file.  This reflects the MOM3 convention of
!      setting the time stamp with the date of the instant after midnight
!      of the last day of the averaging period vs the GRIB convention of
!      using a time stamp indicating the last day of the averaging period.
!     The 3D fields T, S, U, V, W and the 2D fields Isothermal Layer
!      Depth (ILD), Mixed Layer Depth (MLD), surface height (Eta), 
!      wind stress components (TauX, TauY), heat flux and saltflux
!      are extracted from the netCDF file and written to the GRIB file.
!     The netCDF grid has variable spacing in the NS direction, so all
!      fields are interpolated to a uniform 1/3 degree spacing in the 
!      NS direction before writing them to the GRIB file.
!     No interpolation is done in the EW direction. Fields retain their
!      uniform 1 degree grid spacing in the EW direction.
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program nc2grib
!
      include 'netcdf.inc'
!
      integer i,j,k
      integer, parameter :: im=360, jm=200, km=40
      integer imx, jmx, kmx, nf, iprdm
      integer ncid,status
      integer ndims, nvars, natts, idunlm
      integer dsiz, vtype, nvdm, vdm(10), nvatts
      real period, mV
      integer start(4), count(4)
      character*20 dname, vname, aname
      character*80 str
!
      real, dimension(im) :: xg
      real, dimension(jm) :: yg
      real, dimension(km) :: zg
      real, dimension(km) :: pt, sig
      real, dimension(im,jm,km) :: ta, sa
      real, dimension(im,jm) :: arr
!
      character*10 name3(5)
      data name3 /'temp', 'salinity', 'u', 'v', 'w'/
      character*10 fname
      logical tgrd, ugrd, wgrd
!
      character*10 name2(7)
      data name2 /'ild', 'mld', 'eta', 'taux', 'tauy', 'hflx', 'sflx'/
!
      integer, parameter :: jm3=418
      real, parameter :: deps=0.99
      real, dimension(jm3) :: y3
      real, dimension(im,jm3) ::  arj3
!
      integer, dimension(5) :: pid3, ptbl3, plvl3, dsf3
      data pid3 /13, 88, 49, 50, 40/
      data ptbl3 /2, 2, 2, 2, 2/
      data plvl3 /160, 160, 160, 160, 160/
      data dsf3 /3, 6, 4, 4, 8/
      real, dimension(5) :: fctr3, oset3
      data fctr3/1.0, 1.0, 0.01, 0.01, 0.01/
      data oset3/273.15, 0.035, 0.0, 0.0, 0.0/
!
      integer, dimension(7) :: pid2, ptbl2, plvl2, dsf2
      data pid2 /195, 195, 198, 124, 125, 202, 199/
      data ptbl2 /129, 129, 129, 2, 2, 129, 129/
      data plvl2 /238, 237, 1, 1, 1, 1, 1/
      data dsf2 /3, 3, 4, 4, 4, 4, 9/
      real, dimension(7) :: fctr2, oset2
      data fctr2/1.0, 1.0, 0.01, 0.1, 0.1, 41840.0, 1.0/
      data oset2/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      real, parameter :: delT=0.8
      real, parameter :: grav=980.0
      real delSg, sb, ds, rm, rp, th, dt, tb
      integer kb,kbm,kbp,krf
      integer, dimension(im,jm) :: mask
      character(len=4), parameter :: tname = 'temp'
      character(len=8), parameter :: sname = 'salinity'
!
      integer, dimension(200) :: kpds, kgds
      logical*1, dimension(im,jm3) :: lb
      integer ltmn, ltmx, lgmn, lgmx, iret
      integer, parameter :: lugb=50
      integer year, month, day, yr2, cen, lz
      integer p1, p2
!
      character*80 ncFile, grbFile
!
      call w3tagb('GODAS_NC2GRIB',2008,50,50,'NP23')
!
! get command line arguments
!
      narg = iargc()
      if (narg .lt. 1) go to 410
      call getarg(1,ncFile)
!
! get date and period
!
      call getDatePrd(iret)
      go to (210, 220, 230, 240) iret
!
      ln2 = lnstr(ncFile) - 2
      ln1 = ln2 - 2
      ln0 = ln1 - 11
      write(str,'(i4.4,a,i2.2,a,i2.2)') year,'.',month,'.',day
      grbFile = ncFile(1:ln0)//str(1:10)//ncFile(ln1:ln2)//'grb'
!
      p1 = int(period - 0.99)
      p2 = 0
!
      yr2 = mod(year-1,100) + 1
      cen = (year + 99) / 100
!
      nf = im*jm3
!
! set non-varying KPDS
!
      kpds(:) = 0
      kpds(1) = 7                ! NCEP ID
      kpds(2) = 129              ! GODAS ID
      kpds(3) = 255              ! grid specified in GDS
      kpds(4) = 192              ! flags both GDS and BMS
!     kpds(5) -> kpds(7)   set according to field
      kpds(8) = yr2              ! year         ^
      kpds(9) = month            ! month        |
      kpds(10) = day             ! day          Tr
      kpds(11) = 0               ! hour         |
      kpds(12) = 0               ! minute       v
      kpds(13) = 2               ! time unit ID (day)
      kpds(14) = p1              ! P1 
      kpds(15) = p2              ! P2
      kpds(16) = 7               ! Average, (Tr - P1) -> (Tr + P2)
      kpds(17) = p1 + p2 + 1     ! number in ave.
      kpds(18) = 1               ! GRIB version
!     kpds(19)             set according to field
      kpds(20) = 0               ! number missing from average
      kpds(21) = cen             ! century of reference time
!     kpds(22)             set according to field
      kpds(23) = 4               ! sub-center (EMC; NCO=4)
!
! set non-varying KGDS
!
      kgds(:) = 0
      kgds(1) = 0                ! Lat/Long grid
      kgds(2) = im               ! i-dimension
      kgds(3) = jm3              ! j-dimension
!     kgds(4) -> kgds(5)   set according to field
      kgds(6) = 128              ! res flag
!     kgds(7) -> kgds(8)   set according to field
      kgds(9) = 1000             ! Long increment
      kgds(10) = 333             ! Lat increment
      kgds(11) = 64              ! Scanning flag
      kgds(19) = 0               ! nnumber of vertical coord parameters
      kgds(20) = 255             ! no
!
! open grib file
!
      call baopenw(lugb, grbFile, iret)
      if (iret .ne. 0) go to 310
!
! Begin reading 3D fields (T,S,U,V,W) from netCDF files, converting
!  them to GRIB and writing them out
!
      do n3=1,5
!
        fname = name3(n3)
        if (n3 .eq. 1 .or. n3 .eq. 2) then
          tgrd = .true.
          ugrd = .false.
          wgrd = .false.
        else if (n3 .eq. 3 .or. n3 .eq. 4) then
          tgrd = .false.
          ugrd = .true.
          wgrd = .false.
        else
          tgrd = .true.
          ugrd = .false.
          wgrd = .true.
        endif
!
        call ext3dFNc(iret)
        go to (210, 220, 230) iret
!
        ltmn = int(yg(1)*1000.0)
        ltmx =  int(yg(jmx)*1000.0)
        lgmn = int(xg(1)*1000.0)
        lgmx =  int(xg(imx)*1000.0)
!
! set varying KPDS
!
        kpds(5) = pid3(n3)         ! parameter ID 
        kpds(6) = plvl3(n3)        ! type of level, depth in m
        kpds(19) = ptbl3(n3)       ! parameter table ID 
        kpds(22) = dsf3(n3)        ! scale factor
!
! set varying KGDS
!
        kgds(4) = ltmn             ! Lat origin in milli-degrees
        kgds(5) = lgmn             ! Long origin in milli-degrees
        kgds(7) = ltmx             ! Lat max
        kgds(8) = lgmx             ! Long max
!
! cycle through vertical levels
!
        do k=1,kmx
          do j=1,jmx
            do i=1,imx
              arr(i,j) = ta(i,j,k)
              if (arr(i,j) .ne. mV) then
                arr(i,j) = arr(i,j)*fctr3(n3) + oset3(n3)
              endif
            enddo
          enddo
!
! do interpolation in NS direction
!
          if (tgrd) then
            call jt_interp
          else
            call ju_interp
          endif
!
          do j=1,jm3
            do i=1,imx
              if (arj3(i,j) .ne. mV) then
                lb(i,j) = .true.
              else
                lb(i,j) = .false.
                arj3(i,j) = 0.0
              endif
            enddo
          enddo
!
          lz = int(zg(k)+0.001)
          kpds(7) = lz             ! value of depth
!
! write to grib file
!
          iret = 0
          call putgb(lugb,nf,kpds,kgds,lb,arj3,iret)
          if (iret .ne. 0) go to 320
!
        enddo
!
      enddo
!
! Begin reading 2D fields (ILD,MLD,ETS,Tx,Ty,HF,SF) from netCDF files, 
!  converting them to GRIB and writing them out
!
      do n2=1,7
!
        fname = name2(n2)
        if (n2 .le. 3 .or. n2 .ge. 6) then
          tgrd = .true.
          ugrd = .false.
          wgrd = .false.
        else
          tgrd = .false.
          ugrd = .true.
          wgrd = .false.
        endif
!
        if (n2 .eq. 1) then
          call extILDNc(iret)
          go to (210, 220, 230) iret
        else if (n2 .eq. 2) then
          call extMLDNc(iret)
          go to (210, 220, 230) iret
        else
          call ext2dFNc(iret)
          go to (210, 220, 230) iret
        endif
!
        ltmn = int(yg(1)*1000.0)
        ltmx =  int(yg(jmx)*1000.0)
        lgmn = int(xg(1)*1000.0)
        lgmx =  int(xg(imx)*1000.0)
!
! set varying KPDS
!
        kpds(5) = pid2(n2)         ! parameter ID
        kpds(6) = plvl2(n2)        ! type of level, depth in m
        kpds(19) = ptbl2(n2)       ! parameter table ID
        kpds(22) = dsf2(n2)        ! scale factor
!
! set varying KGDS
!
        kgds(4) = ltmn             ! Lat origin in milli-degrees
        kgds(5) = lgmn             ! Long origin in milli-degrees
        kgds(7) = ltmx             ! Lat max
        kgds(8) = lgmx             ! Long max
!
! rescale
!
        do j=1,jmx
          do i=1,imx
            if (arr(i,j) .ne. mV) then
              arr(i,j) = arr(i,j)*fctr2(n2)
            endif
          enddo
        enddo
!
! do interpolation in NS direction
!
        if (tgrd) then
          call jt_interp
        else
          call ju_interp
        endif
!
        do j=1,jm3
          do i=1,imx
            if (arj3(i,j) .ne. mV) then
              lb(i,j) = .true.
            else
              lb(i,j) = .false.
              arj3(i,j) = 0.0
            endif
          enddo
        enddo
!
        kpds(7) = 0              ! value of depth
!
! write to grib file
!
        iret = 0
        call putgb(lugb,nf,kpds,kgds,lb,arj3,iret)
        if (iret .ne. 0) go to 320
!
      enddo
!
! close grib file
!
      call baclose(lugb, iret)
      if (iret .ne. 0) go to 330
!
!
      call w3tage('GODAS_NC2GRIB')
      call errexit(0)
!     call exit(0)
!
  210 write(6,'(a,a)') 'Error opening netCDF file ', ncFile
      call w3tage('GODAS_NC2GRIB')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a,a)') 'Error reading netCDF file ', ncFile
      call w3tage('GODAS_NC2GRIB')
      call errexit(22)
!     call exit(22)
!
  230 write(6,'(a,a)') 'Error closing netCDF file ', ncFile
      call w3tage('GODAS_NC2GRIB')
      call errexit(23)
!     call exit(23)
!
  240 write(6,'(a,a)') 'Error getting date from netCDF file ', ncFile
      call w3tage('GODAS_NC2GRIB')
      call errexit(24)
!     call exit(24)
!
  310 write(6,'(a,a)') 'Error opening GRIB file ', grbFile
      call w3tage('GODAS_NC2GRIB')
      call errexit(31)
!     call exit(31)
!
  320 write(6,'(a,a)') 'Error writing GRIB file ', grbFile
      call w3tage('GODAS_NC2GRIB')
      call errexit(32)
!     call exit(32)
!
  330 write(6,'(a,a)') 'Error closing GRIB file ', grbFile
      call w3tage('GODAS_NC2GRIB')
      call errexit(33)
!     call exit(33)
!
  410 write(6,'(a)') 'Error reading command line'
      call w3tage('GODAS_NC2GRIB')
      call errexit(41)
!     call exit(41)
!
      contains
!
! ---------------------------------------------------------------------
!
      subroutine getDatePrd(ret)
!
      integer ret
!
      integer, dimension(12), save :: dpm
      data dpm/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
      ret = 0
!
! open netCDF dataset
!
      status = nf_open(ncFile, 0, ncid)
      if (status .ne. NF_NOERR) then
        ret = 1
        return
      endif
!
! inquire about the file
!
      status = nf_inq(ncid, ndims, nvars, natts, idunlm)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
!
! inquire about the variables
!
      do n = 1,nvars
        status = nf_inq_var(ncid, n, vname, vtype, nvdm, vdm, nvatts)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (vname .eq. 'period') then
          status = nf_get_var_real(ncid, n, period)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        endif
      enddo
!
! inquire about global attributes
!
      status = nf_inq_varnatts(ncid, NF_GLOBAL, nvatts)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
      do m = 1,nvatts
        status = nf_inq_attname(ncid, NF_GLOBAL, m, aname)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (aname(1:8) .eq. 'filename') then
          status = nf_inq_attlen(ncid, NF_GLOBAL, aname, ln1)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
          status = nf_get_att_text(ncid, NF_GLOBAL, aname, str)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        endif
      enddo
!
! close netCDF file
!
      status = nf_close(ncid)
      if (status .ne. NF_NOERR) then
        ret = 3 
        return
      endif
!
      read(str(11:22),'(i6,1x,i2,1x,i2)',iostat=ios) year, month, day
      if (ios .ne. 0) then
        ret = 4                                                                 
        return
      endif
!
      if (day .eq. 1) then
        if (month .eq. 1) then
          month = 12
          day = 31
          year = year - 1
        else
          month = month - 1
          day = dpm(month)
          if (month .eq. 2) then
            if (leapYear(year)) day = day + 1
          endif
        endif
      else
        day = day - 1
      endif
!
      end subroutine getDatePrd
!
! ---------------------------------------------------------------------
!
      subroutine ext3dFNc(ret)
!
      integer ret
!
      ret = 0
!
! open netCDF dataset
!
      status = nf_open(ncFile, 0, ncid)
      if (status .ne. NF_NOERR) then
        ret = 1
        return
      endif
!
! inquire about the file
!
      status = nf_inq(ncid, ndims, nvars, natts, idunlm)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
!
! inquire about the dimensions
!
      do n = 1,ndims
        status = nf_inq_dim(ncid, n, dname, dsiz)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (dname .eq. 'xt_i') then
          imx = dsiz
        else if (dname .eq. 'yt_j') then
          jmx = dsiz
        else if (dname .eq. 'zt_k') then
          kmx = dsiz
        endif
      enddo
!
! inquire about the variables
! inquire about the variable attributes
!
      do n = 1,nvars
        status = nf_inq_var(ncid, n, vname, vtype, nvdm, vdm, nvatts)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (vname .eq. 'xt_i' .and. tgrd) then
          status = nf_get_var_real(ncid, n, xg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'yt_j' .and. tgrd) then
          status = nf_get_var_real(ncid, n, yg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'xu_i' .and. ugrd) then
          status = nf_get_var_real(ncid, n, xg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'yu_j' .and. ugrd) then
          status = nf_get_var_real(ncid, n, yg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'zt_k' .and. .not. wgrd) then
          status = nf_get_var_real(ncid, n, zg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'zw_k' .and. wgrd) then
          status = nf_get_var_real(ncid, n, zg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. fname) then
          start(1) = 1
          count(1) = imx
          start(2) = 1
          count(2) = jmx
          start(3) = 1
          count(3) = kmx
          start(4) = 1
          count(4) = 1
          status = nf_get_vara_real(ncid, n, start, count, ta)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
          status = nf_get_att_real(ncid, n, 'missing_value', mV)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        endif
      enddo
!
      status = nf_close(ncid)
      if (status .ne. NF_NOERR) then
        ret = 3
        return
      endif
!
      end subroutine ext3dFNc
!
! ---------------------------------------------------------------------
!
      subroutine ext2dFNc(ret)
!
      integer ret
!
! open netCDF dataset
!
      status = nf_open(ncFile, 0, ncid)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
!
! inquire about the file
!
      status = nf_inq(ncid, ndims, nvars, natts, idunlm)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
!
! inquire about the dimensions
!
      do n = 1,ndims
        status = nf_inq_dim(ncid, n, dname, dsiz)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (dname .eq. 'xt_i') then
          imx = dsiz
        else if (dname .eq. 'yt_j') then
          jmx = dsiz
        endif
      enddo
!
! inquire about the variables
! inquire about the variable attributes
!
      do n = 1,nvars
        status = nf_inq_var(ncid, n, vname, vtype, nvdm, vdm, nvatts)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (vname .eq. 'xt_i' .and. tgrd) then
          status = nf_get_var_real(ncid, n, xg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'yt_j' .and. tgrd) then
          status = nf_get_var_real(ncid, n, yg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'xu_i' .and. ugrd) then
          status = nf_get_var_real(ncid, n, xg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'yu_j' .and. ugrd) then
          status = nf_get_var_real(ncid, n, yg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. fname) then
          start(1) = 1
          count(1) = imx
          start(2) = 1
          count(2) = jmx
          start(3) = 1
          count(3) = 1
          if (nvdm == 4) then
            start(3) = 1
            count(3) = 1
            start(4) = 1
            count(4) = 1
          endif
          status = nf_get_vara_real(ncid, n, start, count, arr)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
          status = nf_get_att_real(ncid, n, 'missing_value', mV)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        endif
      enddo
!
      status = nf_close(ncid)
      if (status .ne. NF_NOERR) then
        ret = 3
        return
      endif
!
      end subroutine ext2dFNc
!
! ---------------------------------------------------------------------
!
      subroutine extILDNc(ret)
!
      integer ret
!
! open netCDF dataset
!
      status = nf_open(ncFile, 0, ncid)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
!
! inquire about the file
!
      status = nf_inq(ncid, ndims, nvars, natts, idunlm)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
!
! inquire about the dimensions
!
      do n = 1,ndims
        status = nf_inq_dim(ncid, n, dname, dsiz)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (dname .eq. 'xt_i') then
          imx = dsiz
        else if (dname .eq. 'yt_j') then
          jmx = dsiz
        else if (dname .eq. 'zt_k') then
          kmx = dsiz
        endif
      enddo
!
! inquire about the variables
! inquire about the variable attributes
!
      do n = 1,nvars
        status = nf_inq_var(ncid, n, vname, vtype, nvdm, vdm, nvatts)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (vname .eq. 'xt_i') then
          status = nf_get_var_real(ncid, n, xg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'yt_j') then
          status = nf_get_var_real(ncid, n, yg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'zt_k') then
          status = nf_get_var_real(ncid, n, zg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. tname) then
          start(1) = 1
          count(1) = imx
          start(2) = 1
          count(2) = jmx
          start(3) = 1
          count(3) = kmx
          start(4) = 1
          count(4) = 1
          status = nf_get_vara_real(ncid, n, start, count, ta)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
          status = nf_get_att_real(ncid, n, 'missing_value', mV)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        endif
      enddo
!
      status = nf_close(ncid)
      if (status .ne. NF_NOERR) then
        ret = 3
        return
      endif
!
! make a mask
!
      do j = 1,jmx
        do i = 1,imx
          mask(i,j) = 0
          do k = 1,kmx
            if (ta(i,j,k) .eq. mV) exit
            mask(i,j) = k
          enddo
        enddo
      enddo
!
! compute isothermal layer depth (ILD)
!
      do j = 1,jmx
        do i = 1,imx
          if (mask(i,j) .eq. 0 .or. ta(i,j,1) .lt. 0.0) then
            arr(i,j) = mV
          else
            krf = 1
            kbm = 0
            kbp = 0
            do k = krf,mask(i,j)
              if ((ta(i,j,krf) - ta(i,j,k)) .ge. delT) then
                kbp = k
                exit
              endif
            enddo
            if (kbp .le. 1) then
              arr(i,j) = mV
            else
              kbm = kbp - 1
              tb = ta(i,j,krf) - delT
              dt = ta(i,j,kbp) - ta(i,j,kbm)
              rm = (ta(i,j,kbp) - tb) / dt
              rp = (tb - ta(i,j,kbm)) / dt
              arr(i,j) = zg(kbm)*rm + zg(kbp)*rp
            endif
          endif
        enddo
      enddo
!
      end subroutine extILDNc
!
! ---------------------------------------------------------------------
!
      subroutine extMLDNc(ret)
!
      integer ret
!
! open netCDF dataset
!
      status = nf_open(ncFile, 0, ncid)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
!
! inquire about the file
!
      status = nf_inq(ncid, ndims, nvars, natts, idunlm)
      if (status .ne. NF_NOERR) then
        ret = 2
        return
      endif
!
! inquire about the dimensions
!
      do n = 1,ndims
        status = nf_inq_dim(ncid, n, dname, dsiz)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (dname .eq. 'xt_i') then
          imx = dsiz
        else if (dname .eq. 'yt_j') then
          jmx = dsiz
        else if (dname .eq. 'zt_k') then
          kmx = dsiz
        endif
      enddo
!
! inquire about the variables
! inquire about the variable attributes
!
      do n = 1,nvars
        status = nf_inq_var(ncid, n, vname, vtype, nvdm, vdm, nvatts)
        if (status .ne. NF_NOERR) then
          ret = 2
          return
        endif
        if (vname .eq. 'xt_i') then
          status = nf_get_var_real(ncid, n, xg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'yt_j') then
          status = nf_get_var_real(ncid, n, yg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. 'zt_k') then
          status = nf_get_var_real(ncid, n, zg)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. tname) then
          start(1) = 1
          count(1) = imx
          start(2) = 1
          count(2) = jmx
          start(3) = 1
          count(3) = kmx
          start(4) = 1
          count(4) = 1
          status = nf_get_vara_real(ncid, n, start, count, ta)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
          status = nf_get_att_real(ncid, n, 'missing_value', mV)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        else if (vname .eq. sname) then
          start(1) = 1
          count(1) = imx
          start(2) = 1
          count(2) = jmx
          start(3) = 1
          count(3) = kmx
          start(4) = 1
          count(4) = 1
          status = nf_get_vara_real(ncid, n, start, count, sa)
          if (status .ne. NF_NOERR) then
            ret = 2
            return
          endif
        endif
      enddo
!
      status = nf_close(ncid)
      if (status .ne. NF_NOERR) then
        ret = 3
        return
      endif
!
! compute pressure, ignore latitude
!
      do k = 1,kmx
        pt(k) = press(zg(k),980.0)
      enddo
!
! make a mask
!
      do j = 1,jmx
        do i = 1,imx
          mask(i,j) = 0
          do k = 1,kmx
            if (ta(i,j,k) .eq. mV) exit
            mask(i,j) = k
          enddo
        enddo
      enddo
!
! rescale salinity
!
      do j = 1,jmx
        do i = 1,imx
          do k = 1,mask(i,j)
            sa(i,j,k) = 1000.0*sa(i,j,k) + 35.0
          enddo
        enddo
      enddo
!
! compute mixed layer depth (MLD)
!
      do j = 1,jmx
        do i = 1,imx
          if (mask(i,j) .eq. 0 .or. ta(i,j,1) .lt. 0.0) then
            arr(i,j) = mV
          else
            delSg = (density(0.0,ta(i,j,1)-delT,sa(i,j,1)) -
     &                     density(0.0,ta(i,j,1),sa(i,j,1)))
            do k = 1,mask(i,j)
              th = theta(pt(k),ta(i,j,k),sa(i,j,k),0.0)
              sig(k) = (density(0.0,th,sa(i,j,k)) - 1000.0)
            enddo
            krf = 1
            kbm = 0
            kbp = 0
            do k = krf,mask(i,j)
              if (sig(k) - sig(krf) .ge. delSg) then
                kbp = k
                exit
              endif
            enddo
            if (kbp .le. 1) then
              arr(i,j) = mV
            else
              kbm = kbp - 1
              sb = sig(krf) + delSg
              ds = sig(kbp) - sig(kbm)
              rm = (sig(kbp) - sb) / ds
              rp = (sb - sig(kbm)) / ds
              arr(i,j) = zg(kbm)*rm + zg(kbp)*rp
            endif
          endif
        enddo
      enddo
!
      end subroutine extMLDNc
!
! ---------------------------------------------------------------------
!
      subroutine jt_interp
!
      integer i, j, j0, j1
      real, dimension(jm3), save :: yt3
      data yt3 /-74.500, -74.167, -73.833, -73.500, -73.167,
     &          -72.833, -72.500, -72.167, -71.833, -71.500,
     &          -71.167, -70.833, -70.500, -70.167, -69.833,
     &          -69.500, -69.167, -68.833, -68.500, -68.167,
     &          -67.833, -67.500, -67.167, -66.833, -66.500,
     &          -66.167, -65.833, -65.500, -65.167, -64.833,
     &          -64.500, -64.167, -63.833, -63.500, -63.167,
     &          -62.833, -62.500, -62.167, -61.833, -61.500,
     &          -61.167, -60.833, -60.500, -60.167, -59.833,
     &          -59.500, -59.167, -58.833, -58.500, -58.167,
     &          -57.833, -57.500, -57.167, -56.833, -56.500,
     &          -56.167, -55.833, -55.500, -55.167, -54.833,
     &          -54.500, -54.167, -53.833, -53.500, -53.167,
     &          -52.833, -52.500, -52.167, -51.833, -51.500,
     &          -51.167, -50.833, -50.500, -50.167, -49.833,
     &          -49.500, -49.167, -48.833, -48.500, -48.167,
     &          -47.833, -47.500, -47.167, -46.833, -46.500,
     &          -46.167, -45.833, -45.500, -45.167, -44.833,
     &          -44.500, -44.167, -43.833, -43.500, -43.167,
     &          -42.833, -42.500, -42.167, -41.833, -41.500,
     &          -41.167, -40.833, -40.500, -40.167, -39.833,
     &          -39.500, -39.167, -38.833, -38.500, -38.167,
     &          -37.833, -37.500, -37.167, -36.833, -36.500,
     &          -36.167, -35.833, -35.500, -35.167, -34.833,
     &          -34.500, -34.167, -33.833, -33.500, -33.167,
     &          -32.833, -32.500, -32.167, -31.833, -31.500,
     &          -31.167, -30.833, -30.500, -30.167, -29.833,
     &          -29.500, -29.167, -28.833, -28.500, -28.167,
     &          -27.833, -27.500, -27.167, -26.833, -26.500,
     &          -26.167, -25.833, -25.500, -25.167, -24.833,
     &          -24.500, -24.167, -23.833, -23.500, -23.167,
     &          -22.833, -22.500, -22.167, -21.833, -21.500,
     &          -21.167, -20.833, -20.500, -20.167, -19.833,
     &          -19.500, -19.167, -18.833, -18.500, -18.167,
     &          -17.833, -17.500, -17.167, -16.833, -16.500,
     &          -16.167, -15.833, -15.500, -15.167, -14.833,
     &          -14.500, -14.167, -13.833, -13.500, -13.167,
     &          -12.833, -12.500, -12.167, -11.833, -11.500,
     &          -11.167, -10.833, -10.500, -10.167,  -9.833,
     &           -9.500,  -9.167,  -8.833,  -8.500,  -8.167,
     &           -7.833,  -7.500,  -7.167,  -6.833,  -6.500,
     &           -6.167,  -5.833,  -5.500,  -5.167,  -4.833,
     &           -4.500,  -4.167,  -3.833,  -3.500,  -3.167,
     &           -2.833,  -2.500,  -2.167,  -1.833,  -1.500,
     &           -1.167,  -0.833,  -0.500,  -0.167,   0.167,
     &            0.500,   0.833,   1.167,   1.500,   1.833,
     &            2.167,   2.500,   2.833,   3.167,   3.500,
     &            3.833,   4.167,   4.500,   4.833,   5.167,
     &            5.500,   5.833,   6.167,   6.500,   6.833,
     &            7.167,   7.500,   7.833,   8.167,   8.500,
     &            8.833,   9.167,   9.500,   9.833,  10.167,
     &           10.500,  10.833,  11.167,  11.500,  11.833,
     &           12.167,  12.500,  12.833,  13.167,  13.500,
     &           13.833,  14.167,  14.500,  14.833,  15.167,
     &           15.500,  15.833,  16.167,  16.500,  16.833,
     &           17.167,  17.500,  17.833,  18.167,  18.500,
     &           18.833,  19.167,  19.500,  19.833,  20.167,
     &           20.500,  20.833,  21.167,  21.500,  21.833,
     &           22.167,  22.500,  22.833,  23.167,  23.500,
     &           23.833,  24.167,  24.500,  24.833,  25.167,
     &           25.500,  25.833,  26.167,  26.500,  26.833,
     &           27.167,  27.500,  27.833,  28.167,  28.500,
     &           28.833,  29.167,  29.500,  29.833,  30.167,
     &           30.500,  30.833,  31.167,  31.500,  31.833,
     &           32.167,  32.500,  32.833,  33.167,  33.500,
     &           33.833,  34.167,  34.500,  34.833,  35.167,
     &           35.500,  35.833,  36.167,  36.500,  36.833,
     &           37.167,  37.500,  37.833,  38.167,  38.500,
     &           38.833,  39.167,  39.500,  39.833,  40.167,
     &           40.500,  40.833,  41.167,  41.500,  41.833,
     &           42.167,  42.500,  42.833,  43.167,  43.500,
     &           43.833,  44.167,  44.500,  44.833,  45.167,
     &           45.500,  45.833,  46.167,  46.500,  46.833,
     &           47.167,  47.500,  47.833,  48.167,  48.500,
     &           48.833,  49.167,  49.500,  49.833,  50.167,
     &           50.500,  50.833,  51.167,  51.500,  51.833,
     &           52.167,  52.500,  52.833,  53.167,  53.500,
     &           53.833,  54.167,  54.500,  54.833,  55.167,
     &           55.500,  55.833,  56.167,  56.500,  56.833,
     &           57.167,  57.500,  57.833,  58.167,  58.500,
     &           58.833,  59.167,  59.500,  59.833,  60.167,
     &           60.500,  60.833,  61.167,  61.500,  61.833,
     &           62.167,  62.500,  62.833,  63.167,  63.500,
     &           63.833,  64.167,  64.500/
      integer, dimension(jm3), save :: jt0
      data jt0 /  1,   1,   1,   2,   2,   2,   3,   3,   3,
     &            4,   4,   4,   5,   5,   5,   6,   6,   6,
     &            7,   7,   7,   8,   8,   8,   9,   9,   9,
     &           10,  10,  10,  11,  11,  11,  12,  12,  12,
     &           13,  13,  13,  14,  14,  14,  15,  15,  15,
     &           16,  16,  16,  17,  17,  17,  18,  18,  18,
     &           19,  19,  19,  20,  20,  20,  21,  21,  21,
     &           22,  22,  22,  23,  23,  23,  24,  24,  24,
     &           25,  25,  25,  26,  26,  26,  27,  27,  27,
     &           28,  28,  28,  29,  29,  29,  30,  30,  30,
     &           31,  31,  31,  32,  32,  32,  33,  33,  33,
     &           34,  34,  34,  35,  35,  35,  36,  36,  36,
     &           37,  37,  37,  38,  38,  38,  39,  39,  39,
     &           40,  40,  40,  41,  41,  41,  42,  42,  42,
     &           43,  43,  43,  44,  44,  44,  45,  45,  45,
     &           46,  46,  46,  47,  47,  47,  48,  48,  48,
     &           49,  49,  49,  50,  50,  50,  51,  51,  51,
     &           52,  52,  52,  53,  53,  54,  54,  54,  55,
     &           55,  55,  56,  56,  57,  57,  58,  58,  59,
     &           59,  59,  60,  60,  61,  62,  62,  63,  63,
     &           64,  64,  65,  66,  67,  67,  68,  69,  70,
     &           71,  72,  73,  74,  75,  76,  77,  78,  79,
     &           80,  81,  82,  83,  84,  85,  86,  87,  88,
     &           89,  90,  91,  92,  93,  94,  95,  96,  97,
     &           98,  99, 100, 101, 102, 103, 104, 105, 106,
     &          107, 108, 109, 110, 111, 112, 113, 114, 115,
     &          116, 117, 118, 119, 120, 121, 122, 123, 124,
     &          125, 126, 127, 128, 129, 130, 131, 132, 133,
     &          134, 135, 136, 136, 137, 138, 139, 140, 141,
     &          142, 143, 143, 144, 145, 146, 146, 147, 147,
     &          148, 148, 149, 150, 150, 151, 151, 151, 152,
     &          152, 153, 153, 154, 154, 155, 155, 155, 156,
     &          156, 156, 157, 157, 158, 158, 158, 159, 159,
     &          159, 160, 160, 160, 161, 161, 161, 162, 162,
     &          162, 163, 163, 163, 164, 164, 165, 165, 165,
     &          166, 166, 166, 167, 167, 167, 168, 168, 168,
     &          169, 169, 169, 170, 170, 170, 171, 171, 171,
     &          172, 172, 172, 173, 173, 173, 174, 174, 174,
     &          175, 175, 175, 176, 176, 176, 177, 177, 177,
     &          178, 178, 178, 179, 179, 179, 180, 180, 180,
     &          181, 181, 181, 182, 182, 182, 183, 183, 183,
     &          184, 184, 184, 185, 185, 185, 186, 186, 186,
     &          187, 187, 187, 188, 188, 188, 189, 189, 189,
     &          190, 190, 190, 191, 191, 191, 192, 192, 192,
     &          193, 193, 193, 194, 194, 194, 195, 195, 195,
     &          196, 196, 196, 197, 197, 197, 198, 198, 198,
     &          199, 199, 199, 199/
      integer, dimension(jm3), save :: jt1
      data jt1 /  2,   2,   2,   3,   3,   3,   4,   4,   4,
     &            5,   5,   5,   6,   6,   6,   7,   7,   7,
     &            8,   8,   8,   9,   9,   9,  10,  10,  10,
     &           11,  11,  11,  12,  12,  12,  13,  13,  13,
     &           14,  14,  14,  15,  15,  15,  16,  16,  16,
     &           17,  17,  17,  18,  18,  18,  19,  19,  19,
     &           20,  20,  20,  21,  21,  21,  22,  22,  22,
     &           23,  23,  23,  24,  24,  24,  25,  25,  25,
     &           26,  26,  26,  27,  27,  27,  28,  28,  28,
     &           29,  29,  29,  30,  30,  30,  31,  31,  31,
     &           32,  32,  32,  33,  33,  33,  34,  34,  34,
     &           35,  35,  35,  36,  36,  36,  37,  37,  37,
     &           38,  38,  38,  39,  39,  39,  40,  40,  40,
     &           41,  41,  41,  42,  42,  42,  43,  43,  43,
     &           44,  44,  44,  45,  45,  45,  46,  46,  46,
     &           47,  47,  47,  48,  48,  48,  49,  49,  49,
     &           50,  50,  50,  51,  51,  51,  52,  52,  52,
     &           53,  53,  53,  54,  54,  55,  55,  55,  56,
     &           56,  56,  57,  57,  58,  58,  59,  59,  60,
     &           60,  60,  61,  61,  62,  63,  63,  64,  64,
     &           65,  65,  66,  67,  68,  68,  69,  70,  71,
     &           72,  73,  74,  75,  76,  77,  78,  79,  80,
     &           81,  82,  83,  84,  85,  86,  87,  88,  89,
     &           90,  91,  92,  93,  94,  95,  96,  97,  98,
     &           99, 100, 101, 102, 103, 104, 105, 106, 107,
     &          108, 109, 110, 111, 112, 113, 114, 115, 116,
     &          117, 118, 119, 120, 121, 122, 123, 124, 125,
     &          126, 127, 128, 129, 130, 131, 132, 133, 134,
     &          135, 136, 137, 137, 138, 139, 140, 141, 142,
     &          143, 144, 144, 145, 146, 147, 147, 148, 148,
     &          149, 149, 150, 151, 151, 152, 152, 152, 153,
     &          153, 154, 154, 155, 155, 156, 156, 156, 157,
     &          157, 157, 158, 158, 159, 159, 159, 160, 160,
     &          160, 161, 161, 161, 162, 162, 162, 163, 163,
     &          163, 164, 164, 164, 165, 165, 166, 166, 166,
     &          167, 167, 167, 168, 168, 168, 169, 169, 169,
     &          170, 170, 170, 171, 171, 171, 172, 172, 172,
     &          173, 173, 173, 174, 174, 174, 175, 175, 175,
     &          176, 176, 176, 177, 177, 177, 178, 178, 178,
     &          179, 179, 179, 180, 180, 180, 181, 181, 181,
     &          182, 182, 182, 183, 183, 183, 184, 184, 184,
     &          185, 185, 185, 186, 186, 186, 187, 187, 187,
     &          188, 188, 188, 189, 189, 189, 190, 190, 190,
     &          191, 191, 191, 192, 192, 192, 193, 193, 193,
     &          194, 194, 194, 195, 195, 195, 196, 196, 196,
     &          197, 197, 197, 198, 198, 198, 199, 199, 199,
     &          200, 200, 200, 200/
      real, dimension(jm3), save :: dyt0
      data dyt0/1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.666, 0.332,
     &          0.998, 0.662, 0.327, 0.991, 0.652, 0.313,
     &          0.974, 0.631, 0.287, 0.943, 0.594, 0.245,
     &          0.894, 0.538, 0.182, 0.822, 0.458, 0.093,
     &          0.721, 0.346, 0.971, 0.584, 0.198, 0.805,
     &          0.405, 0.005, 0.590, 0.174, 0.748, 0.315,
     &          0.877, 0.424, 0.970, 0.494, 0.019, 0.520,
     &          0.020, 0.494, 0.964, 0.406, 0.839, 0.248,
     &          0.635, 0.008, 0.341, 0.655, 0.943, 0.192,
     &          0.407, 0.590, 0.738, 0.850, 0.927, 0.973,
     &          0.994, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 0.005, 0.026, 0.072,
     &          0.149, 0.261, 0.409, 0.592, 0.807, 0.056,
     &          0.345, 0.658, 0.992, 0.364, 0.752, 0.160,
     &          0.593, 0.035, 0.506, 0.979, 0.479, 0.980,
     &          0.505, 0.030, 0.576, 0.123, 0.684, 0.251,
     &          0.825, 0.410, 0.995, 0.595, 0.195, 0.802,
     &          0.415, 0.029, 0.653, 0.279, 0.907, 0.542,
     &          0.177, 0.817, 0.461, 0.105, 0.754, 0.405,
     &          0.056, 0.712, 0.369, 0.026, 0.687, 0.348,
     &          0.009, 0.673, 0.337, 0.002, 0.668, 0.334,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 0.000/
      real, dimension(jm3), save :: dyt1
      data dyt1/0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.334, 0.668,
     &          0.002, 0.338, 0.673, 0.009, 0.348, 0.687,
     &          0.026, 0.369, 0.713, 0.057, 0.406, 0.755,
     &          0.106, 0.462, 0.818, 0.178, 0.542, 0.907,
     &          0.279, 0.654, 0.029, 0.416, 0.802, 0.195,
     &          0.595, 0.995, 0.410, 0.826, 0.252, 0.685,
     &          0.123, 0.576, 0.030, 0.506, 0.981, 0.480,
     &          0.980, 0.506, 0.036, 0.594, 0.161, 0.752,
     &          0.365, 0.992, 0.659, 0.345, 0.057, 0.808,
     &          0.593, 0.410, 0.262, 0.150, 0.073, 0.027,
     &          0.006, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.995, 0.974, 0.928,
     &          0.851, 0.739, 0.591, 0.408, 0.193, 0.944,
     &          0.655, 0.342, 0.008, 0.636, 0.248, 0.840,
     &          0.407, 0.965, 0.494, 0.021, 0.521, 0.020,
     &          0.495, 0.970, 0.424, 0.877, 0.316, 0.749,
     &          0.175, 0.590, 0.005, 0.405, 0.805, 0.198,
     &          0.585, 0.971, 0.347, 0.721, 0.093, 0.458,
     &          0.823, 0.183, 0.539, 0.895, 0.246, 0.595,
     &          0.944, 0.288, 0.631, 0.974, 0.313, 0.652,
     &          0.991, 0.327, 0.663, 0.998, 0.332, 0.666,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 1.000/
!
      do j = 1,jm3
        y3(j) = yt3(j)
      enddo
!
      do j = 1,jm3
        j0 = jt0(j)
        j1 = jt1(j)
        if (dyt0(j) .gt. deps) then
          do i = 1,imx
            arj3(i,j) = arr(i,j0)
          enddo
        else if (dyt1(j) .gt. deps) then
          do i = 1,imx
            arj3(i,j) = arr(i,j1)
          enddo
        else
          do i = 1,imx
            if (arr(i,j0) .ne. mV .and. arr(i,j1) .ne. mV) then
              arj3(i,j) = arr(i,j0) * dyt0(j) + arr(i,j1) * dyt1(j)
            else
              arj3(i,j) = mV
            endif
          enddo
        endif
      enddo
!
      end subroutine jt_interp
!
! ---------------------------------------------------------------------
!
      subroutine ju_interp
!
      integer i, j, j0, j1
      real, dimension(jm3), save :: yu3
      data yu3 /-74.000, -73.667, -73.333, -73.000, -72.667,
     &          -72.333, -72.000, -71.667, -71.333, -71.000,
     &          -70.667, -70.333, -70.000, -69.667, -69.333,
     &          -69.000, -68.667, -68.333, -68.000, -67.667,
     &          -67.333, -67.000, -66.667, -66.333, -66.000,
     &          -65.667, -65.333, -65.000, -64.667, -64.333,
     &          -64.000, -63.667, -63.333, -63.000, -62.667,
     &          -62.333, -62.000, -61.667, -61.333, -61.000,
     &          -60.667, -60.333, -60.000, -59.667, -59.333,
     &          -59.000, -58.667, -58.333, -58.000, -57.667,
     &          -57.333, -57.000, -56.667, -56.333, -56.000,
     &          -55.667, -55.333, -55.000, -54.667, -54.333,
     &          -54.000, -53.667, -53.333, -53.000, -52.667,
     &          -52.333, -52.000, -51.667, -51.333, -51.000,
     &          -50.667, -50.333, -50.000, -49.667, -49.333,
     &          -49.000, -48.667, -48.333, -48.000, -47.667,
     &          -47.333, -47.000, -46.667, -46.333, -46.000,
     &          -45.667, -45.333, -45.000, -44.667, -44.333,
     &          -44.000, -43.667, -43.333, -43.000, -42.667,
     &          -42.333, -42.000, -41.667, -41.333, -41.000,
     &          -40.667, -40.333, -40.000, -39.667, -39.333,
     &          -39.000, -38.667, -38.333, -38.000, -37.667,
     &          -37.333, -37.000, -36.667, -36.333, -36.000,
     &          -35.667, -35.333, -35.000, -34.667, -34.333,
     &          -34.000, -33.667, -33.333, -33.000, -32.667,
     &          -32.333, -32.000, -31.667, -31.333, -31.000,
     &          -30.667, -30.333, -30.000, -29.667, -29.333,
     &          -29.000, -28.667, -28.333, -28.000, -27.667,
     &          -27.333, -27.000, -26.667, -26.333, -26.000,
     &          -25.667, -25.333, -25.000, -24.667, -24.333,
     &          -24.000, -23.667, -23.333, -23.000, -22.667,
     &          -22.333, -22.000, -21.667, -21.333, -21.000,
     &          -20.667, -20.333, -20.000, -19.667, -19.333,
     &          -19.000, -18.667, -18.333, -18.000, -17.667,
     &          -17.333, -17.000, -16.667, -16.333, -16.000,
     &          -15.667, -15.333, -15.000, -14.667, -14.333,
     &          -14.000, -13.667, -13.333, -13.000, -12.667,
     &          -12.333, -12.000, -11.667, -11.333, -11.000,
     &          -10.667, -10.333, -10.000,  -9.667,  -9.333,
     &           -9.000,  -8.667,  -8.333,  -8.000,  -7.667,
     &           -7.333,  -7.000,  -6.667,  -6.333,  -6.000,
     &           -5.667,  -5.333,  -5.000,  -4.667,  -4.333,
     &           -4.000,  -3.667,  -3.333,  -3.000,  -2.667,
     &           -2.333,  -2.000,  -1.667,  -1.333,  -1.000,
     &           -0.667,  -0.333,   0.000,   0.333,   0.667,
     &            1.000,   1.333,   1.667,   2.000,   2.333,
     &            2.667,   3.000,   3.333,   3.667,   4.000,
     &            4.333,   4.667,   5.000,   5.333,   5.667,
     &            6.000,   6.333,   6.667,   7.000,   7.333,
     &            7.667,   8.000,   8.333,   8.667,   9.000,
     &            9.333,   9.667,  10.000,  10.333,  10.667,
     &           11.000,  11.333,  11.667,  12.000,  12.333,
     &           12.667,  13.000,  13.333,  13.667,  14.000,
     &           14.333,  14.667,  15.000,  15.333,  15.667,
     &           16.000,  16.333,  16.667,  17.000,  17.333,
     &           17.667,  18.000,  18.333,  18.667,  19.000,
     &           19.333,  19.667,  20.000,  20.333,  20.667,
     &           21.000,  21.333,  21.667,  22.000,  22.333,
     &           22.667,  23.000,  23.333,  23.667,  24.000,
     &           24.333,  24.667,  25.000,  25.333,  25.667,
     &           26.000,  26.333,  26.667,  27.000,  27.333,
     &           27.667,  28.000,  28.333,  28.667,  29.000,
     &           29.333,  29.667,  30.000,  30.333,  30.667,
     &           31.000,  31.333,  31.667,  32.000,  32.333,
     &           32.667,  33.000,  33.333,  33.667,  34.000,
     &           34.333,  34.667,  35.000,  35.333,  35.667,
     &           36.000,  36.333,  36.667,  37.000,  37.333,
     &           37.667,  38.000,  38.333,  38.667,  39.000,
     &           39.333,  39.667,  40.000,  40.333,  40.667,
     &           41.000,  41.333,  41.667,  42.000,  42.333,
     &           42.667,  43.000,  43.333,  43.667,  44.000,
     &           44.333,  44.667,  45.000,  45.333,  45.667,
     &           46.000,  46.333,  46.667,  47.000,  47.333,
     &           47.667,  48.000,  48.333,  48.667,  49.000,
     &           49.333,  49.667,  50.000,  50.333,  50.667,
     &           51.000,  51.333,  51.667,  52.000,  52.333,
     &           52.667,  53.000,  53.333,  53.667,  54.000,
     &           54.333,  54.667,  55.000,  55.333,  55.667,
     &           56.000,  56.333,  56.667,  57.000,  57.333,
     &           57.667,  58.000,  58.333,  58.667,  59.000,
     &           59.333,  59.667,  60.000,  60.333,  60.667,
     &           61.000,  61.333,  61.667,  62.000,  62.333,
     &           62.667,  63.000,  63.333,  63.667,  64.000,
     &           64.333,  64.667,  65.000/
      integer, dimension(jm3), save :: ju0
      data ju0 /  1,   1,   1,   2,   2,   2,   3,   3,   3,
     &            4,   4,   4,   5,   5,   5,   6,   6,   6,
     &            7,   7,   7,   8,   8,   8,   9,   9,   9,
     &           10,  10,  10,  11,  11,  11,  12,  12,  12,
     &           13,  13,  13,  14,  14,  14,  15,  15,  15,
     &           16,  16,  16,  17,  17,  17,  18,  18,  18,
     &           19,  19,  19,  20,  20,  20,  21,  21,  21,
     &           22,  22,  22,  23,  23,  23,  24,  24,  24,
     &           25,  25,  25,  26,  26,  26,  27,  27,  27,
     &           28,  28,  28,  29,  29,  29,  30,  30,  30,
     &           31,  31,  31,  32,  32,  32,  33,  33,  33,
     &           34,  34,  34,  35,  35,  35,  36,  36,  36,
     &           37,  37,  37,  38,  38,  38,  39,  39,  39,
     &           40,  40,  40,  41,  41,  41,  42,  42,  42,
     &           43,  43,  43,  44,  44,  44,  45,  45,  45,
     &           46,  46,  46,  47,  47,  47,  48,  48,  48,
     &           49,  49,  49,  50,  50,  50,  51,  51,  51,
     &           52,  52,  52,  53,  53,  54,  54,  54,  55,
     &           55,  56,  56,  56,  57,  57,  58,  58,  59,
     &           59,  60,  60,  61,  61,  62,  62,  63,  64,
     &           64,  65,  66,  66,  67,  68,  69,  70,  71,
     &           72,  73,  74,  75,  76,  77,  78,  79,  80,
     &           81,  82,  83,  84,  85,  86,  87,  88,  89,
     &           90,  91,  92,  93,  94,  95,  96,  97,  98,
     &           99, 100, 101, 102, 103, 104, 105, 106, 107,
     &          108, 109, 110, 111, 112, 113, 114, 115, 116,
     &          117, 118, 119, 120, 121, 122, 123, 124, 125,
     &          126, 127, 128, 129, 130, 131, 132, 133, 134,
     &          135, 136, 136, 137, 138, 139, 140, 141, 142,
     &          143, 143, 144, 145, 145, 146, 147, 147, 148,
     &          148, 149, 149, 150, 150, 151, 151, 152, 152,
     &          153, 153, 153, 154, 154, 155, 155, 155, 156,
     &          156, 157, 157, 157, 158, 158, 158, 159, 159,
     &          159, 160, 160, 160, 161, 161, 161, 162, 162,
     &          162, 163, 163, 164, 164, 164, 165, 165, 165,
     &          166, 166, 166, 167, 167, 167, 168, 168, 168,
     &          169, 169, 169, 170, 170, 170, 171, 171, 171,
     &          172, 172, 172, 173, 173, 173, 174, 174, 174,
     &          175, 175, 175, 176, 176, 176, 177, 177, 177,
     &          178, 178, 178, 179, 179, 179, 180, 180, 180,
     &          181, 181, 181, 182, 182, 182, 183, 183, 183,
     &          184, 184, 184, 185, 185, 185, 186, 186, 186,
     &          187, 187, 187, 188, 188, 188, 189, 189, 189,
     &          190, 190, 190, 191, 191, 191, 192, 192, 192,
     &          193, 193, 193, 194, 194, 194, 195, 195, 195,
     &          196, 196, 196, 197, 197, 197, 198, 198, 198,
     &          199, 199, 199, 199/
      integer, dimension(jm3), save :: ju1
      data ju1 /  2,   2,   2,   3,   3,   3,   4,   4,   4,
     &            5,   5,   5,   6,   6,   6,   7,   7,   7,
     &            8,   8,   8,   9,   9,   9,  10,  10,  10,
     &           11,  11,  11,  12,  12,  12,  13,  13,  13,
     &           14,  14,  14,  15,  15,  15,  16,  16,  16,
     &           17,  17,  17,  18,  18,  18,  19,  19,  19,
     &           20,  20,  20,  21,  21,  21,  22,  22,  22,
     &           23,  23,  23,  24,  24,  24,  25,  25,  25,
     &           26,  26,  26,  27,  27,  27,  28,  28,  28,
     &           29,  29,  29,  30,  30,  30,  31,  31,  31,
     &           32,  32,  32,  33,  33,  33,  34,  34,  34,
     &           35,  35,  35,  36,  36,  36,  37,  37,  37,
     &           38,  38,  38,  39,  39,  39,  40,  40,  40,
     &           41,  41,  41,  42,  42,  42,  43,  43,  43,
     &           44,  44,  44,  45,  45,  45,  46,  46,  46,
     &           47,  47,  47,  48,  48,  48,  49,  49,  49,
     &           50,  50,  50,  51,  51,  51,  52,  52,  52,
     &           53,  53,  53,  54,  54,  55,  55,  55,  56,
     &           56,  57,  57,  57,  58,  58,  59,  59,  60,
     &           60,  61,  61,  62,  62,  63,  63,  64,  65,
     &           65,  66,  67,  67,  68,  69,  70,  71,  72,
     &           73,  74,  75,  76,  77,  78,  79,  80,  81,
     &           82,  83,  84,  85,  86,  87,  88,  89,  90,
     &           91,  92,  93,  94,  95,  96,  97,  98,  99,
     &          100, 101, 102, 103, 104, 105, 106, 107, 108,
     &          109, 110, 111, 112, 113, 114, 115, 116, 117,
     &          118, 119, 120, 121, 122, 123, 124, 125, 126,
     &          127, 128, 129, 130, 131, 132, 133, 134, 135,
     &          136, 137, 137, 138, 139, 140, 141, 142, 143,
     &          144, 144, 145, 146, 146, 147, 148, 148, 149,
     &          149, 150, 150, 151, 151, 152, 152, 153, 153,
     &          154, 154, 154, 155, 155, 156, 156, 156, 157,
     &          157, 158, 158, 158, 159, 159, 159, 160, 160,
     &          160, 161, 161, 161, 162, 162, 162, 163, 163,
     &          163, 164, 164, 165, 165, 165, 166, 166, 166,
     &          167, 167, 167, 168, 168, 168, 169, 169, 169,
     &          170, 170, 170, 171, 171, 171, 172, 172, 172,
     &          173, 173, 173, 174, 174, 174, 175, 175, 175,
     &          176, 176, 176, 177, 177, 177, 178, 178, 178,
     &          179, 179, 179, 180, 180, 180, 181, 181, 181,
     &          182, 182, 182, 183, 183, 183, 184, 184, 184,
     &          185, 185, 185, 186, 186, 186, 187, 187, 187,
     &          188, 188, 188, 189, 189, 189, 190, 190, 190,
     &          191, 191, 191, 192, 192, 192, 193, 193, 193,
     &          194, 194, 194, 195, 195, 195, 196, 196, 196,
     &          197, 197, 197, 198, 198, 198, 199, 199, 199,
     &          200, 200, 200, 200/
      real, dimension(jm3), save :: dyu0
      data dyu0/1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.667, 0.333, 1.000, 0.667, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.664, 0.329,
     &          0.994, 0.657, 0.320, 0.982, 0.641, 0.300,
     &          0.959, 0.613, 0.267, 0.919, 0.567, 0.214,
     &          0.859, 0.499, 0.138, 0.772, 0.403, 0.033,
     &          0.654, 0.273, 0.889, 0.496, 0.103, 0.699,
     &          0.292, 0.879, 0.455, 0.031, 0.590, 0.147,
     &          0.690, 0.226, 0.751, 0.263, 0.764, 0.251,
     &          0.722, 0.180, 0.616, 0.042, 0.436, 0.816,
     &          0.170, 0.493, 0.795, 0.066, 0.297, 0.497,
     &          0.662, 0.793, 0.888, 0.950, 0.984, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     &          1.000, 1.000, 0.015, 0.049, 0.111, 0.206,
     &          0.337, 0.503, 0.703, 0.934, 0.205, 0.506,
     &          0.830, 0.183, 0.564, 0.957, 0.383, 0.819,
     &          0.277, 0.749, 0.235, 0.736, 0.249, 0.773,
     &          0.309, 0.853, 0.410, 0.968, 0.544, 0.120,
     &          0.708, 0.300, 0.896, 0.503, 0.110, 0.726,
     &          0.346, 0.966, 0.597, 0.227, 0.861, 0.501,
     &          0.141, 0.785, 0.433, 0.081, 0.733, 0.387,
     &          0.041, 0.699, 0.358, 0.017, 0.680, 0.343,
     &          0.005, 0.670, 0.335, 1.000, 0.667, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 1.000, 0.666, 0.333,
     &          1.000, 0.666, 0.333, 0.000/
      real, dimension(jm3), save :: dyu1
      data dyu1/0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.333, 0.667, 0.000, 0.333, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.336, 0.671,
     &          0.006, 0.343, 0.680, 0.018, 0.359, 0.700,
     &          0.041, 0.387, 0.733, 0.081, 0.433, 0.786,
     &          0.141, 0.501, 0.862, 0.228, 0.597, 0.967,
     &          0.346, 0.727, 0.111, 0.504, 0.897, 0.301,
     &          0.708, 0.121, 0.545, 0.969, 0.410, 0.853,
     &          0.310, 0.774, 0.249, 0.737, 0.236, 0.749,
     &          0.278, 0.820, 0.384, 0.958, 0.564, 0.184,
     &          0.830, 0.507, 0.205, 0.934, 0.703, 0.503,
     &          0.338, 0.207, 0.112, 0.050, 0.016, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     &          0.000, 0.000, 0.985, 0.951, 0.889, 0.794,
     &          0.663, 0.497, 0.297, 0.066, 0.795, 0.494,
     &          0.170, 0.817, 0.436, 0.043, 0.617, 0.181,
     &          0.723, 0.251, 0.765, 0.264, 0.751, 0.227,
     &          0.691, 0.147, 0.590, 0.032, 0.456, 0.880,
     &          0.292, 0.700, 0.104, 0.497, 0.890, 0.274,
     &          0.654, 0.034, 0.403, 0.773, 0.139, 0.499,
     &          0.859, 0.215, 0.567, 0.919, 0.267, 0.613,
     &          0.959, 0.301, 0.642, 0.983, 0.320, 0.657,
     &          0.995, 0.330, 0.665, 0.000, 0.333, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 0.000, 0.334, 0.667,
     &          0.000, 0.334, 0.667, 1.000/
!
      do j = 1,jm3
        y3(j) = yu3(j)
      enddo
!
      do j = 1,jm3
        j0 = ju0(j)
        j1 = ju1(j)
        if (dyu0(j) .gt. deps) then
          do i = 1,imx
            arj3(i,j) = arr(i,j0)
          enddo
        else if (dyu1(j) .gt. deps) then
          do i = 1,imx
            arj3(i,j) = arr(i,j1)
          enddo
        else
          do i = 1,imx
            if (arr(i,j0) .ne. mV .and. arr(i,j1) .ne. mV) then
              arj3(i,j) = arr(i,j0) * dyu0(j) + arr(i,j1) * dyu1(j)
            else
              arj3(i,j) = mV
            endif
          enddo
        endif
      enddo
!
      end subroutine ju_interp
!
! ---------------------------------------------------------------------
!
!   Supply depth (zz) in meters and grav acc'l (g) in cm/sec**2
!
      function press(z, g)
!
      integer, parameter :: itr=20
      real p, a0
      real(kind=8) :: e, ae, es
!
      p = z*(1.0076+z*(2.3487e-6-z*1.2887e-11))
      e = zeta(p,g)-z
      ae = abs(e)
      es = ae*2.
      do i = 1,itr
        a0 = 0.972643+p*(1.32696e-5-p*(6.228e-12+p*1.885e-16))
        a0 = a0/(1.0+1.83e-5*p)
        p = p-((g+1.113e-4*p)/a0)*e*0.001
        es = ae
        e = zeta(p,g)-z
        ae = abs(e)
        if (ae .le. 0.01 .or. ae .ge. es) exit
      enddo
!
      if (ae .ge. es) print *, '   PRESS DIVERGENCE'
!
      press = p
!
      end function press
!
! ---------------------------------------------------------------------
!
      function zeta(p, glat)
!
      real p, glat, z
                                                                                
      z = ((-3.434e-12*p+1.113e-7)*p+0.712953)*p
     &                            +14190.7*log(1.0+1.83e-5*p)
      z = (z/(glat+1.113e-4*p))*1000.
                
      zeta = z
!
      end function zeta
!
! ---------------------------------------------------------------------
!
      function theta(p, t, s, pref)
!
      real(kind=8), parameter :: sqrt2 = 0.7071067811865475
      real del_p, del_t1, del_t2, del_t3, del_t4, tp, th
                
      del_p = pref-p
      del_t1 = del_p*atg(p,t,s)
      tp = t+0.5*del_t1
      del_t2 = del_p*atg((p+0.5*del_p),tp,s)
      tp = t+(-0.5+sqrt2)*del_t1+(1.0-sqrt2)*del_t2
      del_t3 = del_p*atg((p+0.5*del_p),tp,s)
      tp = t-sqrt2*del_t2+(1.0+sqrt2)*del_t3
      del_t4 = del_p*atg(pref,tp,s)
      th = (t+(del_t1+(1.0-sqrt2)*del_t2*2.0 +
     &                 (1.0+sqrt2)*del_t3*2.0+del_t4)/6.0)
      theta = th
!
      end function theta
!
! ---------------------------------------------------------------------
!
      function atg(p, t, s)
!
      real ds, a
                
      ds = s-35.0
      a = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p
     &       +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t
     &       +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p
     &       +(-4.2393e-8*t+1.8932e-6)*ds
     &       +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5
!
      atg = a
!
      end function atg
!
! ---------------------------------------------------------------------
!
!     Density is in units of kg/m**3  (1 g/cm**3 = 1000 kg/m**3)
!
      function density(prs, tmp, sal)
!
      real p, t, s, kstp, k0, kw, d0, dw
!
      s = sal
      t = tmp
      p = prs/10.00
!
      kw = 19652.21+(148.4206
     &       -(2.327105-(1.360477e-2-5.155288e-5*t)*t)*t)*t
!
      k0 = kw+s*(54.6746-(0.603459-(1.09987e-2-6.1670e-5*t)*t)*t)
     &       +sqrt(s*s*s)*(7.944e-2+(1.6483e-2-5.3009e-4*t)*t)
!
      kstp = k0+p*((3.239908+(1.43713e-3+(1.16092e-4-5.77905e-7*t)*t)*t)
     &       +s*(2.2838e-3-(1.0981e-5+1.6078e-6*t)*t)
     &       +sqrt(s*s*s)*1.91075e-4
     &       +p*((8.50935e-5-(6.12293e-6-5.2787e-8*t)*t)
     &       -s*(9.9348e-7-(2.0816e-8+9.1697e-10*t)*t)))
!
      dw = 999.842594+(6.793952e-2-(9.095290e-3-(1.001685e-4
     &       -(1.120083e-6-6.536332e-9*t)*t)*t)*t)*t
!
      d0 = dw+s*(0.824493-(4.0899e-3-(7.6438e-5-(8.2467e-7
     &       -5.3875e-9*t)*t)*t)*t)
     &       -sqrt(s*s*s)*(5.72466e-3-(1.0227e-4-1.6546e-6*t)*t)
     &       +s*s*4.8314e-4
!
      density = d0/(1.00-p/kstp)
!
      end function density
!
! ---------------------------------------------------------------------
!
      logical function leapYear(iyr)
!
      integer :: iyr
      logical :: lpyr
!
      if (mod(iyr,400) .eq. 0) then
        lpyr = .true.
      else if (mod(iyr,4) .eq. 0 .and. mod(iyr,100) .ne. 0) then
        lpyr = .true.
      else
        lpyr = .false.
      end if
!
      leapYear = lpyr
!
      end function leapYear
!
! ---------------------------------------------------------------------
!
      function lnstr(s)
!
      character*(*) s
      integer n, nmax
!
      nmax = len(s)
      n = 0
      do while (s(n+1:n+1) .ne. ' ' .and. n .lt. nmax)
        n = n+1
      enddo
!
      lnstr = n
!
      end function lnstr
!
! ---------------------------------------------------------------------
!
      end program nc2grib
