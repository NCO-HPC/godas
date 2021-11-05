!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  DLY2PNTDNC
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2005-08-15
!
! ABSTRACT:  Average daily time_mean files from GODAS into pentad (5-day) 
!             time_mean
!
! PROGRAM HISTORY LOG:
! 2005-08-15  David W. Behringer
!
! USAGE:
!   INPUT FILES:
!     INPUT DAILY TIME_MEAN FILES ARE MANAGED BY NETCDF LIBRARY
!
!   OUTPUT FILES:
!     OUTPUT PENTAD TIME_MEAN FILES ARE MANAGED BY NETCDF LIBRARY
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     NONE
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - lnstr, leapYear, pntdDay
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit, iw3jdn, w3fs26
!       NETCDF   - nf_open, nf_inq, nf_inq_dim, nf_inq_var, nf_inq_attname,
!                  nf_copy_att, nf_inq_varnatts, nf_inq_attlen,
!                  nf_get_att_text, nf_get_vara_real, nf_get_vara_double,
!                  nf_get_att_real,
!                  nf_create, nf_def_dim, nf_def_var, nf_put_att_text,
!                  nf_enddef, nf_put_vara_real,
!                  nf_put_vara_double, nf_close
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - lnstr, leapYear, pntdDay
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit, iw3jdn, w3fs26
!       NETCDF   - nf_open, nf_inq, nf_inq_dim, nf_inq_var, nf_inq_attname,
!                  nf_copy_att, nf_inq_varnatts, nf_inq_attlen,
!                  nf_get_att_text, nf_get_vara_real, nf_get_vara_double,
!                  nf_get_att_real,
!                  nf_create, nf_def_dim, nf_def_var, nf_put_att_text,
!                  nf_enddef, nf_put_vara_real,
!                  nf_put_vara_double, nf_close
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING RUNNING SUM FILE
!     COND =  12 - ERROR WRITING RUNNING SUM FILE
!     COND =  13 - ERROR READING RUNNING SUM FILE
!     COND =  14 - ERROR CLOSING RUNNING SUM FILE
!     COND =  15 - ERROR OPENING SCRATCH FILE
!     COND =  16 - ERROR WRITING SCRATCH FILE
!     COND =  17 - ERROR READING SCRATCH FILE
!     COND =  18 - ERROR CLOSING SCRATCH FILE
!     COND =  21 - ERROR OPENING DAILY TIME_MEAN FILE
!     COND =  22 - ERROR READING DAILY TIME_MEAN FILE
!     COND =  23 - ERROR CLOSING DAILY TIME_MEAN FILE
!     COND =  31 - ERROR OPENING PENTAD TIME_MEAN FILE
!     COND =  32 - ERROR WRITING PENTAD TIME_MEAN FILE
!     COND =  33 - ERROR CLOSING PENTAD TIME_MEAN FILE
!     COND =  41 - ERROR READING FROM COMMAND LINE
!     COND =  51 - ERROR: DATE OF DAILY TIME_MEAN FILE OUT OF SEQUENCE
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     Reads a date YYYYMMDD from command line
!     A sequential flat file is used to maintain a running sum of the
!      daily analyses of the current pentad (SUM_time_mean.00yyyy.pp.P.dta)
!      as well as a log in the form of a list of Julian days.
!     When the input day is the last of a pentad, the sum is completed 
!      for the previous 5 days and an average netCDF file is created.
!     When the input day is the first of a pentad, a new flat summing file
!      is started.
!     When the the input day falls within a pentad, a new daily analysis is
!      added to the running sum flat file.
!     In a leap year, the average file dated March 2, will be the average 
!      6 days rather than the usual 5 days.
!     The program assumes daily files are named according to the
!      convention:  time_mean.00yyyy.mm.dd.nc
!     The output pentad time_mean file will be named according to the
!      convention: time_mean.00yyyy.mm.dd.P.nc
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program Dly2PntdNc
!
      include 'netcdf.inc'
!
      logical leapYear
!
      integer i,j,k,t,nmx
      parameter (i=363,j=203,k=41,t=1)
      parameter (nmx=i*j*k*t)
      integer is,js,ks,nm2,nm3
      parameter (is=360,js=200,ks=40)
      parameter (nm2=is*js)
      parameter (nm3=is*js*ks)
      integer ii,jj,kk,ncid,status
      integer ndims, nvars, natts, idunlm
      integer dsiz(20), vtype, nvdm, vdm(10), nvatts, indx(4)
      integer atype, len
      integer start(4), count(4)
      character*20 dname, vname, aname
      character*80 str
      real xVal, aval, prd
      real*8 tdbl
!
      character*80 dlyFile, pntdFile, sumFile
      integer year, month, day, yrm, nDays, iprd, pDay, nPntd
      integer yr, mn, dy, dow, doy, npt
      integer jdLog(6)
      logical first, last, doSum, newSum, ex
!
      integer ncout
!
      real wght
!
      real temp(nmx), a(nmx)
      real*8 smf2(nm2), smf3(nm3)
!
      call w3tagb('GODAS_DLY2PNTDNC',2008,50,50,'NP23')
!
! get date from command line
!
      narg = iargc()
      if (narg .eq. 0) go to 410
      call getarg(1,str)
      read(str,'(i4,2i2)',err=410) year, month, day
!
      call pntdDay(year,month,day,pDay,nPntd,first,last)
      if (last) then
        doSum = .false.
      else
        doSum = .true.
        if (first) then
          newSum = .true.
        else
          newSum = .false.
        endif
      endif
!
      write(str,'(i4,a1,i2.2,a1,i2.2)') year, '.',month, '.', day
      ln2 = lnstr(str)
      dlyFile = 'time_mean.00'//str(1:ln2)//'.nc'
!
      if (doSum) then
        if (newSum) then
!
! begin a new pentad sum
!
! first check whether previous average was concluded properly
!
          npt = nPntd - 1
          yr = year
          if (npt .lt. 1) then
            npt = 73
            yr = yr - 1
          endif
          write(str,'(i4,a1,i2.2)') yr, '.',npt
          ln2 = lnstr(str)
*          sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.P.dta'
          sumFile = 'SUM_time_mean.P.dta'
          inquire(file=sumFile, exist=ex, iostat = ios)
          if (ex) then
            open(11,file=sumFile,form='UNFORMATTED',access='SEQUENTIAL',
     &          status='OLD',err=110)
            read(11,err=130) nDays
            read(11,err=130) jdLog
            call w3fs26(jdLog(nDays),yr,mn,dy,dow,doy)
            write(6,'(a,i4,a,i2.2,a,i2.2)') 'Previous date: ',
     &                                        yr,'-',mn,'-',dy
            write(6,'(a,i4,a,i2.2,a,i2.2)') 'Current date:  ',
     &                                        year,'-',month,'-',day
            jd1 = iw3jdn(year, month, day)
            jday = jdLog(nDays) + 1
            do while (jday .lt. jd1)
              call w3fs26(jday,yr,mn,dy,dow,doy)
              write(6,'(a,i4,a,i2.2,a,i2.2)') 'Missing date:  ',
     &                                        yr,'-',mn,'-',dy
              jday = jday + 1
            enddo
            go to 510
          endif
!
! create flat file for summing daily fields
!
          write(str,'(i4,a1,i2.2)') year, '.',nPntd
          ln2 = lnstr(str)
*          sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.P.dta'
          sumFile = 'SUM_time_mean.P.dta'
          open(11,file=sumFile,form='UNFORMATTED',access='SEQUENTIAL',
     &          status='NEW',err=110)
          nDays = 1
          write(11, err=120) nDays
          do n = 1,6
            jdLog(n) = 0
          enddo
          jdLog(nDays) = iw3jdn(year, month, day)
          write(11, err=120) jdLog
!
! open daily netCDF dataset
!
          status = nf_open(dlyFile, 0, ncid)
          if (status .ne. NF_NOERR) go to 210
!
! inquire about the file
!
          status = nf_inq(ncid, ndims, nvars, natts, idunlm)
          if (status .ne. NF_NOERR) go to 220
!
! inquire about the dimensions
!
          do n = 1,ndims
            status = nf_inq_dim(ncid, n, dname, dsiz(n))
            if (status .ne. NF_NOERR) go to 220
          enddo
!
! read only variable fields for summing
!
          do nv = 15,31
            status =
     &         nf_inq_var(ncid, nv, vname, vtype, nvdm, vdm, nvatts)
            if (status .ne. NF_NOERR) go to 220
            nd = 1
            do m = 1,nvdm
              start(m) = 1
              count(m) = dsiz(vdm(m))
              nd = nd * dsiz(vdm(m))
            enddo
            status = nf_get_vara_real(ncid, nv, start, count, temp)
            if (status .ne. NF_NOERR) go to 220
            status = nf_get_att_real(ncid, nv, 'missing_value', xVal)
            if (status .ne. NF_NOERR) go to 220
!
            if (nd .eq. nm2) then
              do m=1,nd
                smf2(m) = temp(m)
              enddo
              write(11, err=120) smf2
            else
              do m=1,nd
                smf3(m) = temp(m)
              enddo
              write(11, err=120) smf3
            endif
!
          enddo
          status = nf_close(ncid)
          if (status .ne. NF_NOERR) go to 230
          close(11, err=140)
!
! first day of pentad sum saved
!
        else
!
! continue by adding a day to pentad sum
! open flat file for summing daily fields
!
          write(str,'(i4,a1,i2.2)') year, '.',nPntd
          ln2 = lnstr(str)
*          sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.P.dta'
          sumFile = 'SUM_time_mean.P.dta'
          inquire(file=sumFile, exist=ex, iostat = ios)
          if (.not.ex) then
!
! current summing file missing
! check for previous summing file
!
            npt = nPntd - 1
            yr = year
            if (npt .lt. 1) then
              npt = 73
              yr = yr - 1
            endif
            write(str,'(i4,a1,i2.2)') yr, '.',npt
            ln2 = lnstr(str)
*           sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.P.dta'
            sumFile = 'SUM_time_mean.P.dta'
            inquire(file=sumFile, exist=ex, iostat = ios)
            if (ex) then
              open(11,file=sumFile,form='UNFORMATTED',
     &          access='SEQUENTIAL', status='OLD',err=110)
              read(11,err=130) nDays
              read(11,err=130) jdLog
              call w3fs26(jdLog(nDays),yr,mn,dy,dow,doy)
              write(6,'(a,i4,a,i2.2,a,i2.2)') 'Previous date: ',
     &                                        yr,'-',mn,'-',dy
              write(6,'(a,i4,a,i2.2,a,i2.2)') 'Current date:  ',
     &                                        year,'-',month,'-',day
              jday = jdLog(nDays) + 1
              jd1 = iw3jdn(year, month, day)
              do while (jday .lt. jd1)
                call w3fs26(jday,yr,mn,dy,dow,doy)
                write(6,'(a,i4,a,i2.2,a,i2.2)') 'Missing date:  ',
     &                                        yr,'-',mn,'-',dy
                jday = jday + 1
              enddo
            else
              write(6,'(a,i4,a,i2.2,a,i2.2)') 'Current date:  ',
     &                                        year,'-',month,'-',day
              jd1 = iw3jdn(year, month, day)
              jday = jd1 - pDay + 1
              do while (jday .lt. jd1)
                call w3fs26(jday,yr,mn,dy,dow,doy)
                write(6,'(a,i4,a,i2.2,a,i2.2)') 'Missing date:  ',
     &                                        yr,'-',mn,'-',dy
                jday = jday + 1
              enddo
            endif
            go to 510
          endif
!
          open(11,file=sumFile,form='UNFORMATTED',access='SEQUENTIAL',
     &          status='OLD',err=110)
          read(11,err=130) nDays
          nDays = nDays + 1
          read(11,err=130) jdLog
          jdLog(nDays) = iw3jdn(year, month, day)
!
          if (jdLog(nDays) - jdLog(nDays-1) .gt. 1) then
            call w3fs26(jdLog(nDays-1),yr,mn,dy,dow,doy)
            write(6,'(a,i4,a,i2.2,a,i2.2)') 'Previous date: ',
     &                                        yr,'-',mn,'-',dy
            write(6,'(a,i4,a,i2.2,a,i2.2)') 'Current date:  ',
     &                                        year,'-',month,'-',day
            jday = jdLog(nDays-1) + 1
            do while (jday .lt. jdLog(nDays))
              call w3fs26(jday,yr,mn,dy,dow,doy)
              write(6,'(a,i4,a,i2.2,a,i2.2)') 'Missing date:  ',
     &                                        yr,'-',mn,'-',dy
              jday = jday + 1
            enddo
            go to 510
          else if (jdLog(nDays) .le. jdLog(nDays-1)) then
            call w3fs26(jdLog(nDays-1),yr,mn,dy,dow,doy)
            write(6,'(a,i4,a,i2.2,a,i2.2)') 'Current date:  ',
     &                                        year,'-',month,'-',day
            write(6,'(a)') ' precedes or matches'
            write(6,'(a,i4,a,i2.2,a,i2.2)') 'Previous date: ',
     &                                        yr,'-',mn,'-',dy
            go to 510
          endif
!
          open(12,file='SCRATCH',form='UNFORMATTED',access='SEQUENTIAL',
     &          status='NEW',err=150)
          write(12,err=160) nDays
          write(12,err=160) jdLog
!
! open daily netCDF dataset
!
          status = nf_open(dlyFile, 0, ncid)
          if (status .ne. NF_NOERR) go to 210
!
! inquire about the file
!
          status = nf_inq(ncid, ndims, nvars, natts, idunlm)
          if (status .ne. NF_NOERR) go to 220
!
! inquire about the dimensions
!
          do n = 1,ndims
            status = nf_inq_dim(ncid, n, dname, dsiz(n))
            if (status .ne. NF_NOERR) go to 220
          enddo
!
! read variable fields and sum
!
          do nv = 15,31
            status =
     &         nf_inq_var(ncid, nv, vname, vtype, nvdm, vdm, nvatts)
            if (status .ne. NF_NOERR) go to 220
            nd = 1
            do m = 1,nvdm
              start(m) = 1
              count(m) = dsiz(vdm(m))
              nd = nd * dsiz(vdm(m))
            enddo
            status = nf_get_vara_real(ncid, nv, start, count, temp)
            if (status .ne. NF_NOERR) go to 220
            status = nf_get_att_real(ncid, nv, 'missing_value', xVal)
            if (status .ne. NF_NOERR) go to 220
!
            if (nd .eq. nm2) then
              read(11,err=130) smf2
              do m=1,nd
                if (temp(m) .ne. xVal) then
                  smf2(m) = smf2(m) + temp(m)
                endif
              enddo
              write(12, err=120) smf2
            else
              read(11,err=130) smf3
              do m=1,nd
                if (temp(m) .ne. xVal) then
                  smf3(m) = smf3(m) + temp(m)
                endif
              enddo
              write(12, err=120) smf3
            endif
!
          enddo
          status = nf_close(ncid)
          if (status .ne. NF_NOERR) go to 230
!
          close(11, err=140)
          close(12, err=180)
!
          call system('mv SCRATCH '//sumFile)
!
        endif
!
! one day added to pentad sum and saved
!
      else
!
! finish summing and write a pentad average netCDF file
!
        yrm = year
        if (nPntd .eq. 73) yrm = year - 1
!
        if (leapYear(year) .and. month .eq. 3 .and. day .lt. 5) then
          prd = 6.0
          iprd = 6
        else
          prd = 5.0
          iprd = 5
        endif
!
! open flat file for summing daily fields
!
        write(str,'(i4,a1,i2.2)') yrm, '.',nPntd
        ln2 = lnstr(str)
*        sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.P.dta'
          sumFile = 'SUM_time_mean.P.dta'
        open(11,file=sumFile,form='UNFORMATTED',access='SEQUENTIAL',
     &          status='OLD',err=110)
        read(11,err=130) nDays
        nDays = nDays + 1
        read(11,err=130) jdLog
        jdLog(nDays) = iw3jdn(year, month, day)
!
        if (jdLog(nDays) - jdLog(nDays-1) .gt. 1) then
          call w3fs26(jdLog(nDays-1),yr,mn,dy,dow,doy)
          write(6,'(a,i4,a,i2.2,a,i2.2)') 'Previous date: ',
     &                                        yr,'-',mn,'-',dy
          write(6,'(a,i4,a,i2.2,a,i2.2)') 'Current date:  ',
     &                                        year,'-',month,'-',day
          jday = jdLog(nDays-1) + 1
          do while (jday .lt. jdLog(nDays))
            call w3fs26(jday,yr,mn,dy,dow,doy)
            write(6,'(a,i4,a,i2.2,a,i2.2)') 'Missing date:  ',
     &                                        yr,'-',mn,'-',dy
            jday = jday + 1
          enddo
          go to 510
        endif
!
        if (nDays .lt. iprd) then
          write(6,'(a,i2,a)') 'Warning: ', iprd-nDays,
     &                ' days are missing from pentad average'
        endif
!
        wght = 1.0 / float(nDays)
!
        write(str,'(i4,a1,i2.2,a1,i2.2)') year, '.',month, '.', day
        ln2 = lnstr(str)
        pntdFile = 'time_mean.00'//str(1:ln2)//'.P.nc'
!
! open daily netCDF dataset
!
        status = nf_open(dlyFile, 0, ncid)
        if (status .ne. NF_NOERR) go to 210
!
! create new pentad netCDF dataset
!
        status = nf_create(pntdFile, NF_NOCLOBBER, ncout)
        if (status .ne. NF_NOERR) go to 310
!
! inquire about the file
!
        status = nf_inq(ncid, ndims, nvars, natts, idunlm)
        if (status .ne. NF_NOERR) go to 220
!
! inquire about the dimensions and create dimensions in new dataset
!
        do n = 1,ndims
          status = nf_inq_dim(ncid, n, dname, dsiz(n))
          if (status .ne. NF_NOERR) go to 220
          status = nf_def_dim(ncout, dname, dsiz(n), m)
          if (status .ne. NF_NOERR) go to 320
        enddo
!
! inquire about the variables and create variables in new dataset
! inquire about the variable attributes and create them in new dataset
!
        do n = 1,nvars
          status = nf_inq_var(ncid, n, vname, vtype, nvdm, vdm, nvatts)
          if (status .ne. NF_NOERR) go to 220
          status = nf_def_var(ncout, vname, vtype, nvdm, vdm, m)
          if (status .ne. NF_NOERR) go to 320
          if (nvatts .gt. 0) then
            do m = 1,nvatts
              status = nf_inq_attname(ncid, n, m, aname)
              if (status .ne. NF_NOERR) go to 220
              status = nf_copy_att(ncid, n, aname, ncout, n)
              if (status .ne. NF_NOERR) go to 320
            enddo
          endif
        enddo
!
! inquire about global attributes and create them in new dataset
!
        status = nf_inq_varnatts(ncid, NF_GLOBAL, nvatts)
        if (status .ne. NF_NOERR) go to 220
        do m = 1,nvatts
          status = nf_inq_attname(ncid, NF_GLOBAL, m, aname)
          if (status .ne. NF_NOERR) go to 220
          if (aname(1:8) .eq. 'filename') then
            status = nf_inq_attlen(ncid, NF_GLOBAL, aname, ln1)
            if (status .ne. NF_NOERR) go to 220
            status = nf_get_att_text(ncid, NF_GLOBAL, aname, str)
            if (status .ne. NF_NOERR) go to 220
            ln2 = lnstr(pntdFile)
            status =
     &           nf_put_att_text(ncout, NF_GLOBAL, aname, ln2, pntdFile)
            if (status .ne. NF_NOERR) go to 320
          else
            status =
     &           nf_copy_att(ncid, NF_GLOBAL, aname, ncout, NF_GLOBAL)
            if (status .ne. NF_NOERR) go to 320
          endif
        enddo
!
        status = nf_enddef(ncout)
        if (status .ne. NF_NOERR) go to 320
!
        do n = 1,14
          status = nf_inq_var(ncid, n, vname, vtype, nvdm, vdm, nvatts)
          if (status .ne. NF_NOERR) go to 220
          do m = 1,nvdm
            start(m) = 1
            count(m) = dsiz(vdm(m))
          enddo
          if (vtype .eq. NF_FLOAT) then
            status = nf_get_vara_real(ncid, n, start, count, temp)
            if (status .ne. NF_NOERR) go to 220
            status = nf_put_vara_real(ncout, n, start, count, temp)
            if (status .ne. NF_NOERR) go to 320
          else if (vtype .eq. NF_DOUBLE) then
            status = nf_get_vara_double(ncid, n, start, count, tdbl)
            if (status .ne. NF_NOERR) go to 220
            if (vname(1:4) .eq. 'Time') then
              tdbl = tdbl + 7
            endif
            status = nf_put_vara_double(ncout, n, start, count, tdbl)
            if (status .ne. NF_NOERR) go to 320
          endif
        enddo
!
        n = 32
        status = nf_inq_var(ncid, n, vname, vtype, nvdm, vdm, nvatts)
        if (status .ne. NF_NOERR) go to 220
        start(1) = 1
        count(1) = dsiz(vdm(1))
        status = nf_put_vara_real(ncout, n, start, count, prd)
        if (status .ne. NF_NOERR) go to 320
!
        do nv = 15,31
!
          status = nf_inq_var(ncid, nv, vname, vtype, nvdm, vdm, nvatts)
          if (status .ne. NF_NOERR) go to 220
          nd = 1
          do m = 1,nvdm
            start(m) = 1
            count(m) = dsiz(vdm(m))
            nd = nd * dsiz(vdm(m))
          enddo
          status = nf_get_vara_real(ncid, nv, start, count, temp)
          if (status .ne. NF_NOERR) go to 220
          status = nf_get_att_real(ncid, nv, 'missing_value', xVal)
          if (status .ne. NF_NOERR) go to 220
!
          if (nd .eq. nm2) then
            read(11,err=130) smf2
            do m=1,nd
              if (temp(m) .ne. xVal) then
                a(m) = wght * (temp(m) + smf2(m))
              else
                a(m) = temp(m)
              endif
            enddo
          else
            read(11,err=130) smf3
            do m=1,nd
              if (temp(m) .ne. xVal) then
                a(m) = wght * (temp(m) + smf3(m))
              else
                a(m) = temp(m)
              endif
            enddo
          endif
!
          status = nf_put_vara_real(ncout, nv, start, count, a)
          if (status .ne. NF_NOERR) go to 320
!
        enddo
!
        status = nf_close(ncid)
        if (status .ne. NF_NOERR) go to 230
!
        status = nf_close(ncout)
        if (status .ne. NF_NOERR) go to 330
!
        close(11, err=140)
        call system('rm -f '//sumFile)
!
      endif
!
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a,a)') 'Error opening running sum file ', sumFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a,a)') 'Error writing running sum file ', sumFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(12)
!     call exit(12)
!
  130 write(6,'(a,a)') 'Error reading running sum file ', sumFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(13)
!     call exit(13)
!
  140 write(6,'(a,a)') 'Error closing running sum file ', sumFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(14)
!     call exit(14)
!
  150 write(6,'(a)') 'Error opening scratch file '
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(15)
!     call exit(15)
!
  160 write(6,'(a)') 'Error writing scratch file '
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(16)
!     call exit(16)
!
  170 write(6,'(a)') 'Error reading scratch file '
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(17)
!     call exit(17)
!
  180 write(6,'(a)') 'Error closing scratch file '
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(18)
!     call exit(18)
!
  210 write(6,'(a,a)') 'Error opening time_mean file ', dlyFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a,a)') 'Error reading time_mean file ', dlyFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(22)
!     call exit(22)
!
  230 write(6,'(a,a)') 'Error closing time_mean file ', dlyFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(23)
!     call exit(23)
!
  310 write(6,'(a,a)') 'Error opening time_mean file ', pntdFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(31)
!     call exit(31)
!
  320 write(6,'(a,a)') 'Error writing time_mean file ', pntdFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(32)
!     call exit(32)
!
  330 write(6,'(a,a)') 'Error closing time_mean file ', pntdFile
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(33)
!     call exit(33)
!
  410 write(6,'(a)') 'Error reading date from command line'
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(41)
!     call exit(41)
!
  510 write(6,'(a)') 'Error: date of daily file out of sequence'
      call w3tage('GODAS_DLY2PNTDNC')
      call errexit(51)
!     call exit(51)
!
      end
!
! --------------------------------------------------------------------
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
      end
!
! ---------------------------------------------------------------------
!
      function leapYear (year)
!
      logical :: leapYear
      logical :: lpyr
      integer :: year
!
      if (mod(year,400) .eq. 0) then
        lpyr = .true.
      else if (mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) then
        lpyr = .true.
      else
        lpyr = .false.
      end if
!
      leapYear = lpyr
!
      end
!
! ---------------------------------------------------------------------
!
      subroutine pntdDay (year, month, day, pday, np, f1, l5)
!
      integer :: year, month, day, pday, np
      logical f1, l5
      logical leapYear
      integer yday
!
      nlast = 5
      f1 = .false.
      l5 = .false.
!
      yday = iw3jdn(year, month, day) - iw3jdn(year-1, 12, 31)
!
      if (leapYear(year)) then
        if (yday .eq. 1) then
          pday = 5
          np = 73
        else if (yday .ge. 62) then
          pday = mod(yday-3,5) + 1
          np = (yday-3)/5 + 1
          if (yday .eq. 62) pday = pday + 1
        else
          pday = mod(yday-2,5)+1
          np = (yday-2)/5 + 1
        endif
        if (np .eq. 12) nlast = 6
      else
        if (yday .eq. 1) then
          pday = 5
          np = 73
        else
          pday = mod(yday-2,5)+1
          np = (yday-2)/5 + 1
        endif
      endif
!
      if (pday .eq. 1) then
        f1 = .true.
      else if (pday .eq. nlast) then
        l5 = .true.
      endif
!
      end
