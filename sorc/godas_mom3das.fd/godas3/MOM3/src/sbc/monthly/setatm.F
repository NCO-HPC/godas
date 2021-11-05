#ifdef time_varying_sbc_data
      subroutine setatm (dtatm)
c
c-----------------------------------------------------------------------
c     set up surface forcing data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c     modified by:   m.j. harrison       e-mail-> mjh@gfdl.gov
c-----------------------------------------------------------------------
c
      use mom3_mpp_mod
      use xctdbc_mod
      use sbc_info
      use sbc_ocean
      use timeinterp_mod

      include 'mpif.h'
# include "param.h"
# include "ctdbc.h"
# include "tmngr.h"
# include "switch.h"
# include "taskrows.h"
      character*32 seqrecl
      dimension xxu(imt), yyu(jmt), obc_data(imt,jmt)
      real*4 xxu4(imt),yyu4(jmt),obc_data4(imt,jmt),zz4,aprec4,k4
      integer :: error
#ifdef opr_day
      namelist /ctlatm/ ndysbc
#endif

c
c     names of sequential unformatted input data (S.B.C.)
c
      n=1
      namein(n) = 'taux.mom'
      n=n+1
      namein(n) = 'tauy.mom'
      n=n+1
      namein(n) = 'sst.mom'
      n=n+1
      namein(n) = 'sss.mom'
# if defined monthly_sbc || defined daily_sbc
      n=n+1
      namein(n) = 'fsr.mom'
      n=n+1
      namein(n) = 'flh.mom'
      n=n+1
      namein(n) = 'fsh.mom'
      n=n+1
      namein(n) = 'ful.mom'
      n=n+1
      namein(n) = 'fslt.mom'
#endif
#ifdef esku_flux
      n=n+1
      namein(n) = 'fsr.clim.mom'
      n=n+1
      namein(n) = 'flh.clim.mom'
      n=n+1
      namein(n) = 'fsh.clim.mom'
      n=n+1
      namein(n) = 'ful.clim.mom'
# ifdef fluxcorr
      n=n+1
      namein(n) = 'fcorr.clim.mom'
# endif
#endif

c
c     names of direct access data  (S.B.C.)
c
      m = 1
      nameout(m) = 'taux.ud'
      m=m+1
      nameout(m) = 'tauy.ud'
      m=m+1
      nameout(m) = 'sst.ud'
      m=m+1
      nameout(m) = 'sss.ud'
#if defined monthly_sbc || defined daily_sbc
      m=m+1
      nameout(m) = 'fsr.ud'
      m=m+1
      nameout(m) = 'flh.ud'
      m=m+1
      nameout(m) = 'fsh.ud'
      m=m+1
      nameout(m)= 'ful.ud'
      m=m+1
      nameout(m)= 'fslt.ud'
#endif
#ifdef esku_flux
      m=m+1
      nameout(m) = 'fsr.ud'
      m=m+1
      nameout(m) = 'flh.ud'
      m=m+1
      nameout(m) = 'fsh.ud'
      m=m+1
      nameout(m)= 'ful.ud'
# ifdef fluxcorr
      m=m+1
      nameout(m) = 'fcorr.ud'
# endif
#endif

c
      if (n.ne.ntdbc.or.m.ne.n) then
        write(stdout,*) 'Error in defining S.B.C. files ntdbc= '
     &,                 ntdbc,' n= ',n,'m= ',m
        call abort
      endif
c


      write (stdout,'(/a/)')
     &   ' ==> S.B.C. for MOM are comming from time varying S.B.C. data'
c
c-----------------------------------------------------------------------
c     transfer the  surface boundary conditions
c     from the sequential unformatted files (prepared by PREP_DATA)
c     to direct access files.
c     disk resource factor of 13 is for 12 months + 1 annual mean
c-----------------------------------------------------------------------
c
c     define how dataset is to be used: (periodic or not) and number
c     of records for each dataset
c

#ifdef daily_sbc
c
# ifdef pentad
c    all SBC but sss have the number of records equal to 5 days + 2
c    for the last February pentad in leap years the number is 6 days + 2
# else
#  ifdef opr_day
c    all SBC but sss have the number of records equal to (ndysbc + 2) days
#  else
c    all SBC but sss have the number of records equal to days in month + 2
c    it is assumed the model is run in monthly chunks, i.e. not spanning
c    the end of one month and the start of the next
#  endif
# endif
c    (sst is NCEP weekly analysis, interpolated to daily)
c    sss is climatological 12 months
c
      write (stdout,'(/a,a)') ' Current time: ',stamp
      call get_stamp (stamp, iyr, imon, iday, ihr, imin, isec)
      a_time = set_date (iyr, imon, iday, ihr, imin, isec)
# ifdef pentad
c     ndymo  = days_in_month(a_time)
c     if (ndymo .eq. 29 .and. iday .ge. 25) then
c       ndim = 6
c     else
c       ndim = 5
c     endif
      ndim = 5
      write (stdout,'(a,i2,2x,32a)') ' Pentad length:',ndim, stamp
# else
#  ifdef opr_day
      ndysbc = 1
      call getunit (ioun, 'namelist.ctlatm'
     &,             'formatted sequential rewind')
      read  (ioun, ctlatm)
      call relunit (ioun)
      ndim = ndysbc
      write (stdout,'(a,i2,a)') ' Run period:', ndim, ' day'
#  else
      ndim  = days_in_month(a_time)
      write (stdout,'(a,i3/)') ' Days in month:', ndim
#  endif
# endif
      do n=1,3
        ntdrec(n) = ndim + 2
        period(n) = .false.
      enddo
      n = 4
      ntdrec(n) = 12
      period(n) = .true.
      do n=5,ntdbc
        ntdrec(n) = ndim + 2
        period(n) = .false.
      enddo
#else
# ifdef monthly_sbc
c
c    all SBC but sss have 14 records: Dec of Y-1 + 12 m of Y + Jan of Y+1
c    sss is climatological 12 months
c
      do n=1,3
        ntdrec(n) = 14
        period(n) = .false.
      enddo
      n = 4
      ntdrec(n) = 12
      period(n) = .true.
      do n=5,ntdbc
        ntdrec(n) = 14
        period(n) = .false.
      enddo
# else
      do n=1,4
c
c climatological tau,sst,sss
c
        ntdrec(n) = 12
	period(n) = .true.

      enddo

#  ifdef esku_flux
      do n=5,8
c
c climatological heat fluxes
c
        ntdrec(n) = 12
	period(n) = .true.

      enddo
#   ifdef fluxcorr
      n=9
      ntdrec(n) = 12
      period(n) = .true.
#   endif
#  endif
# endif
#endif


#ifdef reynolds_sst
c
c  reynolds are monthly averages not from climatology. Files 
c  should be padded on the ends (i.e. dec and jan of prior
c  and subsequent years respectively) and "period" is .false.
c
      n=3
      period(n)= .false.
      ntdrec(n)= 14
#endif

#if defined daily_winds 

      write (stdout,'(/a,a)') ' Current time:',stamp
      call get_stamp (stamp, iyr, imon, iday, ihr, imin, isec)
      a_time = set_date (iyr, imon, iday, ihr, imin, isec)
      ndim  = days_in_month(a_time)
      write (stdout,'(a,i3/)') ' Days in month:', ndim
      n=1
      period(n)= .false.
      ntdrec(n)= ndim + 2
      n=2
      period(n)= .false.
      ntdrec(n)= ndim + 2

#endif

#if defined monthly_winds
c
c Winds are monthly averages
c
      n=1
      ntdrec(n) = 14
      n=2
      ntdrec(n) = 14
#endif

c
# if defined cray_t90 
c
c     lrec needs room for character*32 "stamp" + "period" + jrow
c     + latitude row of data = (4 + 1 + 1 + imt) words
c     8 characters per 64bit word = 4 words per character*32 "stamp" 
c
      lrec = 6 + imt
      write (opt_tdbc,'(a,i10,a,i10,a)')
     &  ' unformatted direct words =',lrec,' sds maxsize words='
     &,  lrec*jmt*ntotrec,' novfl'
# else
c
c     lrec needs room for character*32 "stamp" + "period" + jrow
c     + latitude row of data = 2*(4 + 1 + 1 + imt) words
c     4 characters per 32bit word = 8 words per character*32 "stamp"
c     the factor of "2" is to account for 32bit platforms using
c     implicit double precision (real*8) 
c
      lrec = 2*(6 + imt)
      write (opt_tdbc,'(a,i10)') ' unformatted direct words =',lrec
# endif
c

      do n=1,ntdbc

        write(stdout,*) 'NTDREC, N= ',ntdrec(n), n
        if (ntdrec(n) .gt. maxbcr) then
          write (stdout,'(/,1x,a)') 'Error: "ntdrec" exceeds "maxbcr"'
          stop '=>setatm'
        endif

c
        ntotrec = ntdrec(n)

        write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Sequential access disk resource for file ', namein(n)(1:13)
     &,' = ', lrec*jmt*ntotrec*1.e-6,' (MW)'
        write (seqrecl,'(a,i12)') ' words=',imt*jmt + imt + jmt + 4 + 8
        call getunit (nu, namein(n)
     &,     'unformatted sequential rewind ieee'//seqrecl)
# ifdef ssread_sswrite
#  ifdef minimize_sbc_memory
        call ostart (iotdbc(n), nameout(n), imt*jmt*ntotrec, imt, 1)
#  else
        call ostart (iotdbc(n), nameout(n), imt*jmt*ntotrec, imt*jmt, 1)
#  endif
# else
        if (pn .eq. 1) then
          write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Direct access disk resource for file ', nameout(n)(1:13),' = '
     &,  lrec*jmt*ntdrec(n)*1.e-6,' (MW)'
          call getunit (iotdbc(n), nameout(n), opt_tdbc)
        endif
# endif
c
c       transfer the input data to direct access files for each month
c
        write (stdout,'(/a,a)') 'Initializing ',nameout(n)
        cksum = 0.0


        do m=1,ntdrec(n)
          if (n .le. 2 ) then
c
c           read taux & tauy 
c
	    read (nu)
            read (nu) dstamp(m,n),aprec4,iimt,jjmt,nn
     &,            xxu4,yyu4,obc_data4
            aprec(m,n)=aprec4
            xxu(:)=xxu4(:)
            yyu(:)=yyu4(:)
            obc_data(:,:)=obc_data4(:,:)

!         else if (namein(n).eq.'sst.mom') then
          else if (n.eq.3) then
c
c           read sst 
c
#if defined reynolds_sst
	    read (nu)
            read (nu) dstamp(m,n),aprec4,iimt,jjmt,nn
     &,            xxu4,yyu4,obc_data4
            aprec(m,n)=aprec4
            xxu(:)=xxu4(:)
            yyu(:)=yyu4(:)
            obc_data(:,:)=obc_data4(:,:)
#else
	    read(nu)
            read (nu) dstamp(m,n),aprec4,iimt,jjmt,k4,zz4,nn
     &,            xxu4,yyu4,obc_data4
            aprec(m,n)=aprec4
            xxu(:)=xxu4(:)
            yyu(:)=yyu4(:)
            obc_data(:,:)=obc_data4(:,:)
#endif
          else if (n.eq.4) then
c
c	    read sss
c
	    read(nu)
            read (nu) dstamp(m,n),aprec4,iimt,jjmt,k4,zz4,nn
     &,            xxu4,yyu4,obc_data4
            aprec(m,n)=aprec4
            xxu(:)=xxu4(:)
            yyu(:)=yyu4(:)
            obc_data(:,:)=obc_data4(:,:)

#ifdef salinity_psu            
            do i=1,imt
              do j=1,jmt
                if (obc_data(i,j) .ne. 0.0) then
                  obc_data(i,j) = obc_data(i,j)*1.e3 + 35.0
                endif
              enddo
            enddo
#endif
          else
c
            read (nu)
            read (nu) dstamp(m,n),aprec4,iimt,jjmt,nn
     &,            xxu4,yyu4,obc_data4
            aprec(m,n)=aprec4
            xxu(:)=xxu4(:)
            yyu(:)=yyu4(:)
            obc_data(:,:)=obc_data4(:,:)

          endif
c
          cksum = cksum + checksum (obc_data, imt, jmt)
c
c
c         adjust february length when using data with leap year calendar
c
          if (.not. eqyear .and. period(n)) then
            if (nint(aprec(m,n)) .eq. 28) then
              aprec(m,n) = aprec(m,n) + 0.2425
              write (stdout,*)
     &        '=>Warning: adding 0.2425 days to feb for leap year'
     &,       ' calendar for ',nameout(n) 
            endif
	  endif
c
c         time manager cannot work with year 0 so add one to the year
c         (add 2 years to the average for december)
c
          call get_stamp (dstamp(m,n), iyr, imon, iday, ihr, imin, isec)
          if (iyr .eq. 0) then
            if (imon .eq. 1 .and. iday .eq. 1) then
              iyr = iyr + 2
            else
              iyr = iyr + 1
            endif
            call set_stamp(dstamp(m,n),iyr,imon,iday,ihr,imin,isec)
          endif
c
c       all processors read the data files and all write the data
c       this is redundant, but error occured when only one wrote and
c       all tried subsequently to read.
c
          if (pn .eq. 1) then
            call write_obc_data (iotdbc(n), m, dstamp(m,n)
     &,                                 aprec(m,n), obc_data,n) 
          endif
	enddo

c       call mpp_sync()
        call MPI_BARRIER(MPI_COMM_WORLD, error)

        call relunit (nu)

        if (pn .eq. 1) then
          write (stdout,*)
     & ' checksum for s.b.c. = ', cksum, ' on ', nameout(n)
          call relunit (iotdbc(n))
        endif

c
c       construct time structures to mark boundaries of dataset
c
#ifdef daily_sbc
        if (n .eq. 4) then
          call get_stamp (dstamp(1,n), iyr, imon, iday, ihr, imin, isec)
          if (imon == 1) then
            a_time = set_date (iyr-1, 12, iday, ihr, imin, isec)
          else
            a_time = set_date (iyr, imon-1, iday, ihr, imin, isec)
          endif
          ndim  = days_in_month(a_time)
          a_time = set_date (iyr, imon, iday, ihr, imin, isec)
          data_start_time(n) = decrement_date (a_time, 0,0, ndim, 0,0,0)
        else
          call get_stamp (dstamp(1,n), iyr, imon, iday, ihr, imin, isec)
          a_time = set_date (iyr, imon, iday, ihr, imin, isec)
          data_start_time(n) = decrement_date (a_time, 0, 0, 1, 0,0,0)
        endif
#else
        call get_stamp (dstamp(1,n), iyr, imon, iday, ihr, imin, isec)
        if (imon == 1) then
          a_time = set_date (iyr-1, 12, iday, ihr, imin, isec)
        else
          a_time = set_date (iyr, imon-1, iday, ihr, imin, isec)
        endif
        ndim  = days_in_month(a_time)
        a_time = set_date (iyr, imon, iday, ihr, imin, isec)
	data_start_time(n) = decrement_date (a_time, 0, 0, ndim, 0,0,0)
#endif
        call get_date (data_start_time(n),iyr, imon, iday,ihr,imin,isec)
        call set_stamp (rstamp1, iyr, imon, iday, ihr, imin, isec)
c
        k = ntdrec(n)
        call get_stamp (dstamp(k,n), iyr, imon, iday, ihr, imin, isec) 
        data_end_time(n) = set_date (iyr, imon, iday, ihr, imin, isec)
c
c       construct center of data records in days and show data bounds
c
        do m=1,ntdrec(n)
	  if (m .eq. 1) then
	    tdrec(1,n) = 0.5*aprec(1,n)
	  else
	    tdrec(m,n) = tdrec(m-1,n) + 0.5*(aprec(m-1,n) + aprec(m,n))
	    rstamp1 = dstamp(m-1,n)
	  endif
	  write (stdout,'(a,i3,a,f10.4,4a)') ' data record # ',m
     &,   ' is centered at ',tdrec(m,n), 'days. defined from '
     &,   rstamp1,' to ',dstamp(m,n)
	enddo

      enddo


# if defined minimize_sbc_memory
c
      do n=1,ntdbc 
c
c       express model time as number of days into dataset
c
        daysbc(n) = position_within_data (model_time, data_start_time(n)
     &,                              data_end_time(n), period(n))
#ifdef debug_atmos
        write (stdout,*) ' n=',n,' daysbc(n)=',daysbc(n)
#endif

c
c       initialize the pointers (indices) to memory buffers and choose
c       interpolation method
c
        iprevm(n) = 1
        inextm(n) = 2
	method    = 3
c
c       find the disk indices and interpolation factor
c
        indxtd(n) = -1
        call timeinterp (daysbc(n), indxtd(n), tdrec(:,n), aprec(:,n)
     &,      ntdrec(n), period(n), method, inextd(n), iprevd(n)
     &,      wprev(n), rdtdbc(n), inextm(n), iprevm(n))
      enddo
# else
      write (stdout,'(/a,a,a/)') '==> Note: to save memory when using'
     &,' "time_varying_sbc_data", enable "minimize_sbc_memory"'
     &,' which will use disk instead'
c
c-----------------------------------------------------------------------
c     determine time interpolation factors and disk record indices for 
c     interpolating time dependent surface boundary conditions to the
c     model time. the S.B.C. will be updated at the mid points of each 
c     month. For generality, each S.B.C. is treated separately.
c     This (for example) allows for a mixture of monthly heatflux with
c     daily winds.
c
c     note: to understand how the time interpolator works, use
c           run_timeinterp.
c-----------------------------------------------------------------------
c
      do n=1,ntdbc 
c
c
c       express model time as number of days into dataset
c
        daysbc(n) = position_within_data (model_time, data_start_time(n)
     &,                              data_end_time(n), period(n))
#ifdef debug_atmos
        write (stdout,*) ' n=',n,' daysbc(n)=',daysbc(n)
#endif
c
c       initialize the pointers (indices) to memory buffers and choose
c       interpolation method
c
        iprevm(n) = 1
        inextm(n) = 2
	method    = 3
c
c       find the disk indices and interpolation factor
c
        indxtd(n) = -1
        call timeinterp (daysbc(n), indxtd(n), tdrec(:,n), aprec(:,n)
     &,      ntdrec(n), period(n), method, inextd(n), iprevd(n)
     &,      wprev(n), rdtdbc(n), inextm(n), iprevm(n))
c
c       read data straddling (just before and after) current model time
c
        write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Direct access disk resource for file ', nameout(n)(1:13),' = '
     &,  lrec*jmt*ntdrec(n)*1.e-6,' (MW)'
        call getunit (iotdbc(n), nameout(n), opt_tdbc)
c
        js = 1
	je = jmt
	joff = 0
        call read_obc_data (iotdbc(n), iprevd(n), rstamp1
     &,                     obc(1,1,n,iprevm(n)), js, je, joff)
        call read_obc_data (iotdbc(n), inextd(n), rstamp2
     &,                     obc(1,1,n,inextm(n)), js, je, joff)
c
        write (stdout,8800) n, iprevd(n), rstamp1, iprevm(n) 
     &,                        inextd(n), rstamp2, inextm(n)
     &,                        nameout(n), stamp
        write (stdout,'(/24x,a12,a,i2,a,i2/)')
     &  trim(nameout(n)),' is dataset index ',n
     &, ' for time interpolation using method #', method
        write (stdout,*)
     &   '  ...1st cksum =', checksum (obc(1,1,n,iprevm(n)), imt*jmt, 1)
     &,  '  ...2nd cksum =', checksum (obc(1,1,n,inextm(n)), imt*jmt, 1)
      enddo
8800  format (/1x,'==> Reading S.B.C. #',i2,' from rec ',i6, ' at ', a
     &,' into index ',i3,/24x,' and rec ',i6,' at ', a, ' into index '
     &,i3,/,24x, ' for ',a12, ' at MOM time ',a) 
# endif
c    
c
c-----------------------------------------------------------------------
c     zero all fluxes in case a zero flux S.B.C. is desired
c     "m" relates how S.B.C. "n" is ordered within "sbcocn"
c-----------------------------------------------------------------------
c
      do n=1,numosbc
        do jrow=1,jmt
          do i=1,imt
            sbcocn(i,jrow,n) = 0.0
          enddo
        enddo
      enddo

# ifdef coupled
c
c     verify that "coupled" is not enabled
c
      write (stdout,*) '==> Error: "time_varying_sbc_data" and '
     &, '"coupled" are incompatible'
      stop '=>setatm'
# else
c
c-----------------------------------------------------------------------
c     since this is not an atmospheric model, the atmospheric 
c     time step "dtatm" does not have to be set. mom.F will set it.
c-----------------------------------------------------------------------
c
# endif
# ifdef simple_sbc
      write (stdout,'(/,a,a,/)')
     &'==> Error: "time_varying_sbc_data" and '
     &,' "simple_sbc" are incompatible'
      stop '=>setatm'
# endif
c
c
c-----------------------------------------------------------------------
c     initialize euler2 because it is tested within atmos but not set
c     until the first call to subroutine mom
c-----------------------------------------------------------------------
c
      euler2 = .false.
      return
      end

      subroutine write_obc_data (iou, m, stamp, period, obc, n)
# if defined sponges && !defined sponges_old
      use xsponge_mod
# endif
# include "size.h"
# include "stdunits.h"
# include "taskrows.h"
      dimension obc(imt,jmt), buf(imt)
      character*(*) stamp
# if defined sponges && !defined sponges_old
      if (n .eq. 1 .or. n .eq. 2) then
        do jrow = 2,j_south_spng-2
	  do i=1,imt
	    obc(i,jrow) = obc(i,j_south_spng-1)
	  enddo
          if (n .eq. 1 .and. m .eq. 1) then
	    write (stdout,*) '=>Note: resetting stress on row'
     &,     jrow,' to stress on row', j_south_spng-1
          endif
	enddo
        do jrow = j_north_spng+1,jmt-1
	  do i=1,imt
	    obc(i,jrow) = obc(i,j_north_spng)
	  enddo
          if (n .eq. 1 .and. m .eq. 1) then
	    write (stdout,*) '=>Note: resetting stress on row'
     &,     jrow,' to stress on row', j_north_spng
          endif
	enddo
      else
        do jrow = 2,j_south_spng-1
	  do i=1,imt
	    obc(i,jrow) = obc(i,j_south_spng)
	  enddo
          if (n .eq. 3 .and. m .eq. 1) then
	    write (stdout,*) '=>Note: resetting SST & SSS & fluxes on row'
     &,     jrow,' to SST & SSS on row', j_south_spng
          endif
	enddo
        do jrow = j_north_spng+1,jmt-1
	  do i=1,imt
	    obc(i,jrow) = obc(i,j_north_spng)
	  enddo
          if (n .eq. 3 .and. m .eq. 1) then
	    write (stdout,*) '=>Note: resetting SST & SSS on row'
     &,     jrow,' to SST & SSS on row', j_north_spng
          endif
	enddo
      endif
# endif
# ifdef ssread_sswrite
      call oput (iou, imt*jmt, m, obc)
# else
      do j=1,jmt
        nrec = (m-1)*jmt + j
        do i=1,imt
	  buf(i) = obc(i,j)
        enddo
        write (iou, rec=nrec) stamp, period, j, buf
      enddo
# endif
      write (stdout,'(a,i2,1x,a,a,g14.7,a)') ' record=',m, stamp
     &,  '.  averaged over preceeding ', period, 'days'
      return
      end


      subroutine read_obc_data (iou, m, stamp, obc, js, je, joff)
# include "size.h"
      dimension obc(imt,js:je), buf(imt)
      character*(*) stamp
# ifdef ssread_sswrite
#  ifdef minimize_sbc_memory
      do j=js,je
        jrow = j + joff
        nrec = (m-1)*jmt + jrow
        call oget (iou, imt, nrec, obc(1,j))
      enddo
#  else
      call oget (iou, imt*jmt, m, obc)
#  endif
      stamp = "(not avail with ssread_sswrite) "
# else
      do j=js,je
        jrow = j + joff
        nrec = (m-1)*jmt + jrow
        read (iou, rec=nrec) stamp, period, jj, buf
        do i=1,imt
	  obc(i,j) = buf(i)
        enddo
      enddo
# endif
      return
      end


#else
      subroutine nodatasetatm
      return
      end
#endif
