!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  DLY2MNTHNC
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2005-08-15
!
! ABSTRACT:  Average daily time_mean files from GODAS into monthly time_mean
!
! PROGRAM HISTORY LOG:
! 2005-08-15  David W. Behringer
!
! USAGE:
!   INPUT FILES:
!     INPUT DAILY TIME_MEAN FILES ARE MANAGED BY NETCDF LIBRARY
!
!   OUTPUT FILES:
!     OUTPUT MONTHLY TIME_MEAN FILES ARE MANAGED BY NETCDF LIBRARY
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     A SEQUENTIAL FLAT FILE IS USED TO MAINTAIN A RUNNING SUM OF THE
!       THE DAILY INPUT FILES
!     A TEMPORARY SCRATCH FILE IS USED TO UPDATE THE RUNNING SUM FILE
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - lnstr, leapYear
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
!     UNIQUE:    - lnstr, leapYear
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
!     COND =  31 - ERROR OPENING MONTHLY TIME_MEAN FILE
!     COND =  32 - ERROR WRITING MONTHLY TIME_MEAN FILE
!     COND =  33 - ERROR CLOSING MONTHLY TIME_MEAN FILE
!     COND =  41 - ERROR READING FROM COMMAND LINE
!     COND =  51 - ERROR: DATE OF DAILY TIME_MEAN FILE OUT OF SEQUENCE 
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     Reads a date YYYYMMDD from command line
!     A sequential flat file is used to maintain a running sum of the
!      daily analyses of the current month (SUM_time_mean.00yyyy.mm.M.dta)
!      as well as a log in the form of a list of Julian days.
!     When the input day equals 1, the sum is completed for the previous
!      month and an average netCDF file is created for the previous month.
!     When the input day equals 2, a new flat summing file is started.
!     When the the input day is greater than 2, a new daily analysis is
!      added to the running sum flat file.
!     The program assumes daily files are named according to the
!      convention:  time_mean.00yyyy.mm.dd.nc
!     The output monthly time_mean file will be named according to the
!      convention: time_mean.00yyyy.mm.dd.M.nc
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program Dly2MnthNc
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
      real dpm(12)
      data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0,
     &          31.0, 31.0, 30.0, 31.0, 30.0, 31.0/
!
      character*80 dlyFile, mnthFile, sumFile
      integer year, month, yrm, monm, day, nDays, iprd
      integer yr, mn, dy, dow, doy
      integer jdLog(31)
      logical doSum, newSum, ex
!
      integer ncout
!
      real wght
!
      real temp(nmx), a(nmx)
      real*8 smf2(nm2), smf3(nm3)
!
      call w3tagb('GODAS_DLY2MNTHNC',2008,50,50,'NP23')
!
! get date from command line
!
      narg = iargc()
      if (narg .eq. 0) go to 410
      call getarg(1,str)
      read(str,'(i4,2i2)',err=410) year, month, day
!
      if (day .eq. 1) then
        doSum = .false.
      else
        doSum = .true.
        if (day .eq. 2) then
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
! begin a new monthly sum
!
! first check whether previous average was concluded properly
!
          mn = month - 1
          yr = year
          if (mn .lt. 1) then
            mn = 12
            yr = yr - 1
          endif
          write(str,'(i4,a1,i2.2)') yr, '.',mn
          ln2 = lnstr(str)
*         sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.M.dta'
          sumFile = 'SUM_time_mean.M.dta'
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
          write(str,'(i4,a1,i2.2)') year, '.',month
          ln2 = lnstr(str)
*          sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.M.dta'
          sumFile = 'SUM_time_mean.M.dta'
          open(11,file=sumFile,form='UNFORMATTED',access='SEQUENTIAL',
     &          status='NEW',err=110)
          nDays = 1
          write(11, err=120) nDays
          do n = 1,31
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
! first day of monthly sum saved
!
        else
!
! continue by adding a day to monthly sum
! open flat file for summing daily fields
!
          write(str,'(i4,a1,i2.2)') year, '.',month
          ln2 = lnstr(str)
*          sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.M.dta'
          sumFile = 'SUM_time_mean.M.dta'
          inquire(file=sumFile, exist=ex, iostat = ios)
          if (.not.ex) then
!
! current summing file missing
! check for previous summing file
!
            mn = month - 1
            yr = year
            if (mn .lt. 1) then
              mn = 12
              yr = yr - 1
            endif
            write(str,'(i4,a1,i2.2)') yr, '.',mn
            ln2 = lnstr(str)
*            sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.M.dta'
             sumFile = 'SUM_time_mean.M.dta'
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
              jday = iw3jdn(year, month, 2)
              jd1 = iw3jdn(year, month, day)
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
! one day added to monthly sum and saved
!
      else
!
! finish summing and write a monthly average netCDF file
!
        monm = month - 1
        yrm = year
        if (monm < 1) then
          monm = 12
          yrm = year - 1
        endif
        prd = dpm(monm)
        if (leapYear(yrm) .and. month .eq. 3) then
          prd = prd + 1
        endif
        iprd = int(prd + 0.001)
!
! open flat file for summing daily fields
!
        write(str,'(i4,a1,i2.2)') yrm, '.',monm
        ln2 = lnstr(str)
*        sumFile = 'SUM_time_mean.00'//str(1:ln2)//'.M.dta'
        sumFile = 'SUM_time_mean.M.dta'
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
     &                ' days are missing from monthly average'
        endif
!
        wght = 1.0 / float(nDays)
!
        write(str,'(i4,a1,i2.2,a1,i2.2)') year, '.',month, '.', day
        ln2 = lnstr(str)
        mnthFile = 'time_mean.00'//str(1:ln2)//'.M.nc'
!
! open daily netCDF dataset
!
        status = nf_open(dlyFile, 0, ncid)
        if (status .ne. NF_NOERR) go to 210
!
! create new monthly netCDF dataset
!
        status = nf_create(mnthFile, NF_NOCLOBBER, ncout)
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
            ln2 = lnstr(mnthFile)
            status = 
     &           nf_put_att_text(ncout, NF_GLOBAL, aname, ln2, mnthFile)
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
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a,a)') 'Error opening running sum file ', sumFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a,a)') 'Error writing running sum file ', sumFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(12)
!     call exit(12)
!
  130 write(6,'(a,a)') 'Error reading running sum file ', sumFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(13)
!     call exit(13)
!
  140 write(6,'(a,a)') 'Error closing running sum file ', sumFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(14)
!     call exit(14)
!
  150 write(6,'(a)') 'Error opening scratch file '
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(15)
!     call exit(15)
!
  160 write(6,'(a)') 'Error writing scratch file '
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(16)
!     call exit(16)
!
  170 write(6,'(a)') 'Error reading scratch file '
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(17)
!     call exit(17)
!
  180 write(6,'(a)') 'Error closing scratch file '
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(18)
!     call exit(18)
!
  210 write(6,'(a,a)') 'Error opening time_mean file ', dlyFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a,a)') 'Error reading time_mean file ', dlyFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(22)
!     call exit(22)
!
  230 write(6,'(a,a)') 'Error closing time_mean file ', dlyFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(23)
!     call exit(23)
!
  310 write(6,'(a,a)') 'Error opening time_mean file ', mnthFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(31)
!     call exit(31)
!
  320 write(6,'(a,a)') 'Error writing time_mean file ', mnthFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(32)
!     call exit(32)
!
  330 write(6,'(a,a)') 'Error closing time_mean file ', mnthFile
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(33)
!     call exit(33)
!
  410 write(6,'(a)') 'Error reading date from command line'
      call w3tage('GODAS_DLY2MNTHNC')
      call errexit(41)
!     call exit(41)
!
  510 write(6,'(a)') 'Error: date of daily file out of sequence'
      call w3tage('GODAS_DLY2MNTHNC')
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
