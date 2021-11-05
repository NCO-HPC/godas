      subroutine setassim
#ifdef assm_dta
!
!-----------------------------------------------------------------------
!     This version of the assimilation code has been written to use 
!     only MPI message passing.
!-----------------------------------------------------------------------
!     Set up the weekly observational data for assimilation.
!
!     author:   d.w.behringer   david.behringer@noaa.gov
!
!-----------------------------------------------------------------------
!
      use mom3_mpp_mod
      use arrays_2d
      use assim_mod
      use asm_x_mod
      use time_manager_mod
      use timeassim_mod
      use timeinterp_mod
!
      include 'mpif.h'
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "tmngr.h"
# include "switch.h"
# include "taskrows.h"
!
      integer, dimension(:), allocatable :: len, off
      integer :: pnw(4), error
      real(kind=4) :: rt, obs, ai, aj, ak, err
      real(kind=4) :: vprec4
      type(time_type) :: t_time
      character(len=8) :: cll
      character(len=20) :: sst_fle, tmp_fle, sal_fle, ssh_fle
      real, dimension(:), allocatable :: wgns
      real, dimension(:,:), allocatable :: s1
      real(kind=4), dimension(:,:,:), allocatable :: bufijk
      real(kind=4), dimension(:,:), allocatable :: bufij
      character(len=128) :: opt_vv
!
      logical assrestrt
      namelist /asmcntrl/ acoef, vcvs, hrscl, rassint, nassdur
     &,                   assrestrt, mitrDB, pnvw
!
      imtka = imt*kass2
!
      do n = 1,4
        pnw(n) = min(n,num_processors)
      enddo
!----------------------------------------------------------------------
!     set number of levels per processor.  setassim must be called after
!     the 1D horizontal decomposition has been done.
!----------------------------------------------------------------------
      num_processors_v = num_processors
      pnv = pn
      nk = kass / num_processors_v
!
      ks = 1
      do n = 1, num_processors_v
        ke = ks + CEILING(float(kass-ks+1)/(num_processors_v-n+1)) - 1
        kscomp3(n) = ks
        kecomp3(n) = ke
        ks = ke + 1
      enddo
      kscomp = kscomp3(pnv)
      kecomp = kecomp3(pnv)
!
!----------------------------------------------------------------------
!     set intervals for assimilation switch
!----------------------------------------------------------------------
      rassint = 0.5
      nassdur = 3
!----------------------------------------------------------------------
!     set global scaling for model error variance
!----------------------------------------------------------------------
      acoef = 0.01
!----------------------------------------------------------------------
!     define constant horizontal scale of model error covariance (deg)
!----------------------------------------------------------------------
      hrscl = 3.99
!----------------------------------------------------------------------
!     set model vertical covariance scale
!----------------------------------------------------------------------
      vcvs = 25.0e2
!----------------------------------------------------------------------
!     set default for assimilation at restart
!----------------------------------------------------------------------
      assrestrt = .false.
!
      mitrDB = maxits
      pnvw = num_processors / 2
!----------------------------------------------------------------------
!     allow for resetting parameters by "namelist"
!-----------------------------------------------------------------------
!
      call getunit (ioun, 'namelist.asmcntrl'
     &,             'formatted sequential rewind')
      read  (ioun, asmcntrl)
      write (stdout, asmcntrl)
      call relunit (ioun)
      call getunit (iodoc, 'document.dta'
     &,             'formatted sequential append')
      write (iodoc, asmcntrl)
      call relunit (iodoc)
      if (assrestrt) then
        nassim = 2
      else
        nassim = nassdur + 1
      endif

!----------------------------------------------------------------------
!     set vertical covariance matrix
!----------------------------------------------------------------------
!
      r2 = vcvs * vcvs
      rsr = sqrt(0.5*vcvs)
      do kk = 1,kass
        do k = kk,kass
          dz2 = (zt(k) - zt(kk)) * (zt(k) - zt(kk))
          cvn(k,kk) = 0.85*exp(-(sqrt(rsr*dz2/r2)))
          cvn(kk,k) = cvn(k,kk)
#if defined cor_sal || defined asm_sal
!         cvnsalt(k,kk) = 0.85*exp(-(sqrt(rsr*dz2/r2)))
          cvnsalt(k,kk) = cvn(k,kk)
          cvnsalt(kk,k) = cvnsalt(k,kk)
#endif
        enddo
      enddo
!
#ifdef fix_vv
!----------------------------------------------------------------------
!     read in fixed in time but geo-varying error variance
!----------------------------------------------------------------------
!
      allocate( bufijk(1:imt,1:jmt,1:kass) )
      allocate( vtmp(kass,imt,jscomp:jecomp) )
      call getunit (nu,'tvv.mom','unformatted sequential rewind ieee')
      read (nu) vstamp0,vprec4,iimt,jjmt,kkss,bufijk
      do j = jscomp,jecomp
        do i = 1,imt
          do k = 1,kass
            vtmp(k,i,j) = bufijk(i,j,k)
          enddo
        enddo
      enddo
      call relunit (nu)
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
# if defined cor_sal || defined asm_sal
      allocate( vsal(kass,imt,jscomp:jecomp) )
      call getunit (nu,'svv.mom','unformatted sequential rewind ieee')
      read (nu) vstamp0,vprec4,iimt,jjmt,kkss,bufijk
      do j = jscomp,jecomp
        do i = 1,imt
          do k = 1,kass
            vsal(k,i,j) = bufijk(i,j,k)
          enddo
        enddo
      enddo
      call relunit (nu)
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
# endif
      deallocate( bufijk )
#else

!----------------------------------------------------------------------
!     write the geo-varying error variance to direct access files
!----------------------------------------------------------------------
!
      allocate( vtmp(kass,imt,jscomp:jecomp,2) )
!
      nvvrec = 14
      vvprd = .false.
      allocate( bufijk(1:imt,1:jmt,1:kass) )
      allocate( bufik(1:imt,1:kass) )
      call getunit (nu,'tvv.mom','unformatted sequential rewind ieee')
      lrec = imt*kass
      write (opt_vv,'(a,i10)') ' unformatted direct words =',lrec
      call getunit (iotvv,'tvv.ud',opt_vv)
      do m = 1,nvvrec
        read (nu) vstamp(m),vprec4,iimt,jjmt,kkss,bufijk
        vprec(m) = vprec4
        do j = 1,jmt
          do i = 1,imt
            do k = 1,kass
              bufik(i,k) = bufijk(i,j,k)
            enddo
          enddo
          nrec = j + (m-1)*jmt
          write (iotvv,rec=nrec) bufik
        enddo
      enddo
      call relunit (nu)
c
      call MPI_BARRIER(MPI_COMM_WORLD, error)
c
#if defined cor_sal || defined asm_sal
!
!   it is assumed that svv and tvv have the same time structures
!
      allocate( vsal(kass,imt,jscomp:jecomp,2) )
!
      call getunit (nu,'svv.mom','unformatted sequential rewind ieee')
      lrec = imt*kass
      write (opt_vv,'(a,i10)') ' unformatted direct words =',lrec
      call getunit (iosvv,'svv.ud',opt_vv)
      do m = 1,nvvrec
        read (nu) vstamp(m),vprec4,iimt,jjmt,kkss,bufijk
        vprec(m) = vprec4
        do j = 1,jmt
          do i = 1,imt
            do k = 1,kass
              bufik(i,k) = bufijk(i,j,k)
            enddo
          enddo
          nrec = j + (m-1)*jmt
          write (iosvv,rec=nrec) bufik
        enddo
      enddo
      call relunit (nu)
c
      call MPI_BARRIER(MPI_COMM_WORLD, error)
#endif
c
      deallocate( bufijk )
#endif
c
#ifdef asm_ssh
!
!  read in mean climatology for model sea surface height
!
      allocate( sshc(1:imt,jstask:jetask) )
      sshc = 0.0
      allocate( bufij(imt,jmt) )
      if (pn == 1) then
        call getunit(nu,'sshc.mom','unformatted sequential rewind ieee')
        read (nu) bufij
        call relunit (nu)
      endif
!
      call MPI_BCAST(bufij, nwds, MPI_REAL4, 0, MPI_COMM_WORLD, error)
      do j=jstask,jetask
        do i=1,imt
          sshc(i,j) = bufij(i,j)
        enddo
      enddo
      deallocate( bufij )
!
!  read in coefficients for linear computation of dynamic height
!  and broadcast them to all processors
!
      allocate( cdnz(km) )
# if defined cor_sal || defined asm_sal
      allocate( cdnzs(km) )
# endif
      if (pn == 1) then
        call getunit(nu,'cdnz.mom','unformatted sequential rewind ieee')
        read (nu) cdnz
# if defined cor_sal || defined asm_sal
        read (nu) cdnzs
# endif
        call relunit (nu)
      else
        cdnz = 0
# if defined cor_sal || defined asm_sal
        cdnzs = 0
# endif
      endif
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      call MPI_BCAST(cdnz, km, MPI_REAL8, 0, MPI_COMM_WORLD, error)
# if defined cor_sal || defined asm_sal
      call MPI_BCAST(cdnzs, km, MPI_REAL8, 0, MPI_COMM_WORLD, error)
# endif
      write (stdout,'(a)') 'Coefficients for linear Dynamic Height'
      write (stdout,'(a)') ' Temperature'
      write (stdout,'(5f13.6)') cdnz
# if defined cor_sal || defined asm_sal
      write (stdout,'(a)') ' Salinity'
      write (stdout,'(5f13.4)') cdnzs
# endif
#endif
c
#ifndef fix_vv
c     construct time structures to mark boundaries of dataset
c
      call get_stamp (vstamp(1), iyr, imon, iday, ihr, imin, isec)
      if (imon == 1) then
        v_time = set_date (iyr-1, 12, iday, ihr, imin, isec)
      else
        v_time = set_date (iyr, imon-1, iday, ihr, imin, isec)
      endif
      ndim  = days_in_month(v_time)
      v_time = set_date (iyr, imon, iday, ihr, imin, isec)
      vv_start_time = decrement_date (v_time, 0, 0, ndim, 0,0,0)
      call get_date (vv_start_time,iyr, imon, iday,ihr,imin,isec)
      call set_stamp (vstamp0, iyr, imon, iday, ihr, imin, isec)
c
      k = nvvrec
      call get_stamp (vstamp(k), iyr, imon, iday, ihr, imin, isec)
      vv_end_time = set_date (iyr, imon, iday, ihr, imin, isec)
c
c       construct center of data records in days and show data bounds
c
      do m=1,nvvrec
        if (m .eq. 1) then
          vvrec(1) = 0.5*vprec(1)
        else
          vvrec(m) = vvrec(m-1) + 0.5*(vprec(m-1) + vprec(m))
          vstamp0 = vstamp(m-1)
        endif
        write (stdout,'(a,i3,a,f10.4,4a)') ' vv record # ',m
     &,   ' is centered at ',vvrec(m), 'days. defined from '
     &,   vstamp0,' to ',vstamp(m)
      enddo
c
c-----------------------------------------------------------------------
c     determine time interpolation factors and disk record indices for
c     interpolating time dependent error variances to the model time. 
c     the variances will be updated at the mid points of each month.
c-----------------------------------------------------------------------
c
c     express model time as number of days into dataset
c
      dayvv = position_within_data (model_time, vv_start_time
     &,                              vv_end_time, vvprd)
c
c     initialize the pointers (indices) to memory buffers and choose
c     interpolation method
c
      iprvvm = 1
      inxtvm = 2
      mthdvv = 3
c
c     find the disk indices and interpolation factor
c
      indxvv = -1
      call timeinterp (dayvv, indxvv, vvrec, vprec, nvvrec, vvprd
     &,        mthdvv, inxtvd, iprvvd, wprvv, rdvv, inxtvm, iprvvm)
c
c     read data just before and after the current model time
c
      do j = jscomp,jecomp
        nrec = j + (iprvvd-1)*jmt
        read (iotvv,rec=nrec) bufik
        do i = 1,imt
          do k = 1,kass
            vtmp(k,i,j,iprvvm) = bufik(i,k)
          enddo
        enddo
      enddo
      do j = jscomp,jecomp
        nrec = j + (inxtvd-1)*jmt
        read (iotvv,rec=nrec) bufik
        do i = 1,imt
          do k = 1,kass
            vtmp(k,i,j,inxtvm) = bufik(i,k)
          enddo
        enddo
      enddo
c
      write (stdout,9000) iprvvd, iprvvm, inxtvd, inxtvm, stamp
9000  format (/1x,'==> Reading temperature variance from rec ',i3
     &, ' to index ',i3,' and from rec ',i3,' to index ',i3,/
     &, 5x,' at MOM time ',a)
!
#if defined cor_sal || defined asm_sal
!
      do j = jscomp,jecomp
        nrec = j + (iprvvd-1)*jmt
        read (iosvv,rec=nrec) bufik
        do i = 1,imt
          do k = 1,kass
            vsal(k,i,j,iprvvm) = bufik(i,k)
          enddo
        enddo
      enddo
      do j = jscomp,jecomp
        nrec = j + (inxtvd-1)*jmt
        read (iosvv,rec=nrec) bufik
        do i = 1,imt
          do k = 1,kass
            vsal(k,i,j,inxtvm) = bufik(i,k)
          enddo
        enddo
      enddo
c
      write (stdout,9001) iprvvd, iprvvm, inxtvd, inxtvm, stamp
9001  format (/1x,'==> Reading salinity variance from rec ',i3
     &, ' to index ',i3,' and from rec ',i3,' to index ',i3,/
     &, 5x,' at MOM time ',a)
!
#endif
#endif
#ifdef exRes
!
!----------------------------------------------------------------------
!     open units for saving residuals
!----------------------------------------------------------------------
!
      write (tresd,'(a,i2.2)') 'TResd.',pn
      call getunit (ioextr, tresd, 'formatted sequential rewind')
# ifdef asm_sal
      write (sresd,'(a,i2.2)') 'SResd.',pn
      call getunit (ioexsr, sresd, 'formatted sequential rewind')
# endif
# ifdef asm_ssh                                                                 
      write (eresd,'(a,i2.2)') 'EResd.',pn
      call getunit (ioexer, eresd, 'formatted sequential rewind')
# endif
#endif
!
!----------------------------------------------------------------------
!     allocate various arrays for assimilation
!----------------------------------------------------------------------
!
      allocate( ares(1:imt,1:kass2,jstask:jetask) )
      allocate( arex(1:imt,1:kass2,jstask:jetask) )
      ares = 0.0
      arex = 0.0
# ifdef asm_ssh
      allocate( etax(1:imt,jstask:jetask) )
      etax = 0.0
      allocate( dssh(1:imt,jstask:jetask) )
      dssh = 0.0
# endif
!
!     arrays to be used by the laplacian smoother. It is decomposed
!     in the vertical, so these arrays are allocated with the global 
!     N-S dimension (jmt). The exception is s2v which facillitates
!     multiplication with the vertical error covariance matrix. The
!     buffer s2buf is used to transpose data between the horizontal
!     and vertical versions of s2.
!
      allocate( wgta(jmt) )
      wgta = 0.0
      allocate( elipt(jmt) )
      elipt = 0.0
      allocate( wew(imt,jmt) )
      allocate( wno(imt,jmt) )
      allocate( wso(imt,jmt) )
      wew = 0.0
      wno = 0.0
      wso = 0.0
      allocate( scra(imt,jmt) )
      allocate( temp(imt,jmt) )
      scra = 0.0
      temp = 0.0
      allocate( dpth(imt,jmt) )
      allocate( kma(imt,jmt) )
      dpth = 0.0
      kma = 0
      allocate( s2h(imt,jmt,kscomp:kecomp) )
      allocate( s2v(imt,kass2,jscomp:jecomp) )
      s2h = 0.0
      s2v = 0.0
      jsz = jecomp3(1) - jscomp3(1) + 1
      ksz = kecomp3(1) - kscomp3(1) + 1
      do n=2,num_processors
        if ( jsz < jecomp3(n)-jscomp3(n)+1 )
     &                           jsz = jecomp3(n) - jscomp3(n) + 1
        if ( ksz < kecomp3(n)-kscomp3(n)+1 )
     &                           ksz = kecomp3(n) - kscomp3(n) + 1
      enddo
      allocate( s2buf(imt,jsz,ksz) )
      s2buf = 0.0
      ns2b = imt*jsz*ksz
!
!     arrays for computation
!
      allocate( t_cg(imt,kass2,jstask:jetask) )
      allocate( d_cg(imt,kass2,jstask:jetask) )
      allocate( e_cg(imt,kass2,jstask:jetask) )
      allocate( f_cg(imt,kass2,jstask:jetask) )
      allocate( g_cg(imt,kass2,jstask:jetask) )
      allocate( h_cg(imt,kass2,jstask:jetask) )
      allocate( r_cg(imt,kass2) )
!
!     arrays for internal management of observations
!
      allocate( rtm(mobs) )
      allocate( val(mobs) )
      allocate( aip(mobs) )
      allocate( ajp(mobs) )
      allocate( akp(mobs) )
      allocate( aerr(mobs) )
      allocate( ov(mobs) )
      allocate( a111(mobs) )
      allocate( a211(mobs) )
      allocate( a221(mobs) )
      allocate( a121(mobs) )
      allocate( ils(mobs) )
      allocate( nobs(jscomp:jecomp) )
      allocate( nsobs(jscomp:jecomp+1) )
      nr = jecomp - jscomp + 1
      allocate( srtm(nr*mobs) )
      allocate( sval(nr*mobs) )
      allocate( sip(nr*mobs) )
      allocate( sjp(nr*mobs) )
      allocate( skp(nr*mobs) )
      allocate( serr(nr*mobs) )
!
!----------------------------------------------------------------------
!     gather kmt into a global array (kma) for each processor
!----------------------------------------------------------------------
      kma(:,:) = 0
      allocate( len(num_processors) )
      allocate( off(num_processors) )
      do n = 1,num_processors
        len(n) = (jecomp3(n) - jscomp3(n) + 1) * imt
        off(n) = (jscomp3(n) - 1) * imt
      enddo
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      call MPI_ALLGATHERV(kmt(:,jscomp),len(pn),MPI_INTEGER,
     &                    kma(:,:),len,off,MPI_INTEGER,
     &                    MPI_COMM_WORLD,error)
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      deallocate( len )
      deallocate( off )
!
!----------------------------------------------------------------------
!  Begin reading in observations. Each processor keeps only the data
!  within its range.  The observations are written to direct access
!  disk files.  They are organized by model row and arrays containing
!  the number of obs per row (*tobs(:,1,:)) and offsets to the start
!  of each row's data (*tobs(:,2,:)) are computed here.
!
!  Time stamps (obstmp) are retained on all processors.  They
!  are used to identify the weekly data periods currently being
!  assimilated.  For now there is no flexibility.  All data sets are
!  organized in weekly periods, all have the same number of periods
!  and all have the same starting and ending dates.  Thus, the time
!  stamps must be the same far all data sets.
!----------------------------------------------------------------------
!
#ifdef asm_sst
!
!----------------------------------------------------------------------
!  read and save sst data to direct access file
!----------------------------------------------------------------------
!
      do n=1,num_processors
        if (pn == n) then
          allocate( itobs(jscomp:jecomp,2,nwkobs) )
          itobs(:,:,:) = 0
!
!         write(stdout,'(/,a,i3)')
          write(6,'(a,i3)')
     & ' starting to read and sort SST data on pn=', pn
          call flush(6)
!
          write (sst_fle,'(a,i2.2)') 'ssta.ud.', pn
          write (opt_dasm,'(a)') ' unformatted direct words = 6'
          call getunit (iosstao, sst_fle, opt_dasm)
          call getunit (nu, 'ssta.mom'
     &,                      'unformatted sequential rewind ieee')
          ak = 1.0
          iak = 1
          joff = 0
          do nwk = 1,nwkobs
            read (nu) obstmp(nwk), icnt
!           write(stdout,'(a10,i3,a7,i5)') ' week no. ',nwk,' icnt= ',icnt
!
!       write(6,'(a10,i3,a7,i5)') ' week no. ',nwk,' icnt= ',icnt
!       call flush(6)
            if (icnt.gt.0) then
              jobs = 0
              jrow = jscomp
              do ic = 1,icnt
                read (nu) cll,rt,obs,ai,aj,err
                iai=ai
                iaj=aj
!
                ierr = 0
!
! Retain only data needed by each processor *********************
                if (iaj < jscomp .or. iaj > jecomp) ierr = ierr + 1
!
! Drop some data close to topography ****************************
                if (iak .gt. kass) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj-1)) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj-1)) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj-1)) ierr = ierr + 1
! End of topography check****************************************
!
                if (ierr == 0) then
                  if (iaj /= jrow) then
                    itobs(jrow,1,nwk) = jobs
                    itobs(jrow,2,nwk) = joff
                    if (jobs > 0) then
                      do j=1,jobs
                        write(iosstao,REC=joff+j)(aobs(i,j),i=1,6)
                      enddo
                    endif
                    joff = joff + jobs
                    do j = jrow+1,iaj-1
                      itobs(j,2,nwk) = joff
                    enddo
                    jobs = 0
                    jrow = iaj
                  endif
                  jobs = jobs + 1
      if (jobs.gt.mobs) then
      write(6,'(a)') 'MOBS exceeded'
      call flush(6)
      endif
                  aobs(1,jobs) = rt + 3.5 + 7.0 * float(nwk-1)
                  aobs(2,jobs) = obs
                  aobs(3,jobs) = ai
                  aobs(4,jobs) = aj
                  aobs(5,jobs) = ak
                  aobs(6,jobs) = err
                endif
              enddo
              itobs(jrow,1,nwk) = jobs
              itobs(jrow,2,nwk) = joff
              if (jobs > 0) then
                do j=1,jobs
                  write(iosstao,REC=joff+j)(aobs(i,j),i=1,6)
                enddo
              endif
              joff = joff + jobs
              do j = jrow+1,jecomp
                itobs(j,2,nwk) = joff
              enddo
            endif
!
          enddo
!
          call relunit (nu)
!         write(stdout,'(i3,a)')
!    &      nwkobs,' weeks of sst obs. processed and saved'
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, error)
      enddo
#endif
#ifdef asm_tmp
!----------------------------------------------------------------------
!  read and save T(z) data to direct access file
!----------------------------------------------------------------------
!
      do n=1,num_processors
        if (pn == n) then
          allocate( jtobs(jscomp:jecomp,2,nwkobs) )
          jtobs(:,:,:) = 0
!
!         write(stdout,'(/,a,i3)') 
          write(6,'(a,i3)') 
     & ' starting to read and sort T(z) data on pn=', pn
          call flush(6)
!
          write (tmp_fle,'(a,i2.2)') 'tmpa.ud.', pn
          write (opt_dasm,'(a)') ' unformatted direct words = 6'
          call getunit (iotmpao, tmp_fle, opt_dasm)
          call getunit (nu, 'tmpa.mom'
     &,                    'unformatted sequential rewind ieee')
          joff = 0
          do nwk = 1,nwkobs
            read (nu) obstmp(nwk), icnt
!           write(stdout,'(a10,i3,a7,i5)') ' week no. ',nwk,' icnt= ',icnt
            if (icnt.gt.0) then
              jobs = 0
              jrow = jscomp
              do ic = 1,icnt
                read (nu) cll,rt,obs,ai,aj,ak,err
                iai=ai
                iaj=aj
                iak=ak
!
                ierr = 0
!
! Retain only data needed by each processor *********************
                if (iaj < jscomp .or. iaj > jecomp) ierr = ierr + 1
!
! Drop some data close to topography ****************************
                if (iak .gt. kass) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj-1)) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj-1)) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj-1)) ierr = ierr + 1
! End of topography check****************************************
!
                if (ierr.eq.0) then
                  if (iaj.ne.jrow) then
                    jtobs(jrow,1,nwk) = jobs
                    jtobs(jrow,2,nwk) = joff
                    if (jobs > 0) then
                      do j=1,jobs
                        write(iotmpao,REC=joff+j)(aobs(i,j),i=1,6) 
                      enddo
                    endif
                    joff = joff + jobs
                    do j = jrow+1,iaj-1
                      jtobs(j,2,nwk) = joff
                    enddo
                    jobs = 0
                    jrow = iaj
                  endif
                  jobs = jobs + 1
      if (jobs.gt.mobs) then
      write(6,'(a)') 'MOBS exceeded'
      call flush(6)
      endif
                  aobs(1,jobs) = rt + 3.5 + 7.0 * float(nwk-1)
                  aobs(2,jobs) = obs
                  aobs(3,jobs) = ai
                  aobs(4,jobs) = aj
                  aobs(5,jobs) = ak
                  aobs(6,jobs) = err
                endif
              enddo
              jtobs(jrow,1,nwk) = jobs
              jtobs(jrow,2,nwk) = joff
              if (jobs > 0) then
                do j=1,jobs
                  write(iotmpao,REC=joff+j)(aobs(i,j),i=1,6) 
                enddo
              endif
              joff = joff + jobs
              do j = jrow+1,jecomp
                jtobs(j,2,nwk) = joff
              enddo
            endif
!
          enddo
!
          call relunit (nu)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, error)
      enddo
#endif
#ifdef asm_sal
!----------------------------------------------------------------------
!  read and save S(z) data to direct access file
!----------------------------------------------------------------------
!
      do n=1,num_processors
        if (pn == n) then
          allocate( ktobs(jscomp:jecomp,2,nwkobs) )
          ktobs(:,:,:) = 0
!
!         write(stdout,'(/,a,i3)')
          write(6,'(a,i3)')
     & ' starting to read and sort S(z) data on pn=', pn
          call flush(6)
!
          write (sal_fle,'(a,i2.2)') 'sala.ud.', pn
          write (opt_dasm,'(a)') ' unformatted direct words = 6'
          call getunit (iosalao, sal_fle, opt_dasm)
          call getunit (nu, 'sala.mom'
     &,                    'unformatted sequential rewind ieee')
          joff = 0
          do nwk = 1,nwkobs
            read (nu) obstmp(nwk), icnt
!           write(stdout,'(a10,i3,a7,i5)') ' week no. ',nwk,' icnt= ',icnt
            if (icnt.gt.0) then
              jobs = 0
              jrow = jscomp
              do ic = 1,icnt
                read (nu) cll,rt,obs,ai,aj,ak,err
                iai=ai
                iaj=aj
                iak=ak
!
                ierr = 0
!
! Retain only data needed by each processor *********************
                if (iaj < jscomp .or. iaj > jecomp) ierr = ierr + 1
!
! Drop some data close to topography ****************************
                if (iak .gt. kass) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj-1)) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj-1)) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj-1)) ierr = ierr + 1
! End of topography check****************************************
!
                if (ierr.eq.0) then
                  if (iaj.ne.jrow) then
                    ktobs(jrow,1,nwk) = jobs
                    ktobs(jrow,2,nwk) = joff
                    if (jobs > 0) then
                      do j=1,jobs
                        write(iosalao,REC=joff+j)(aobs(i,j),i=1,6) 
                      enddo
                    endif
                    joff = joff + jobs
                    do j = jrow+1,iaj-1
                      ktobs(j,2,nwk) = joff
                    enddo
                    jobs = 0
                    jrow = iaj
                  endif
                  jobs = jobs + 1
      if (jobs.gt.mobs) then
      write(6,'(a)') 'MOBS exceeded'
      call flush(6)
      endif
                  aobs(1,jobs) = rt + 3.5 + 7.0 * float(nwk-1)
                  aobs(2,jobs) = obs
                  aobs(3,jobs) = ai
                  aobs(4,jobs) = aj
                  aobs(5,jobs) = ak + float(kass)
                  aobs(6,jobs) = err
                endif
              enddo
              ktobs(jrow,1,nwk) = jobs
              ktobs(jrow,2,nwk) = joff
              if (jobs > 0) then
                do j=1,jobs
                  write(iosalao,REC=joff+j)(aobs(i,j),i=1,6) 
                enddo
              endif
              joff = joff + jobs
              do j = jrow+1,jecomp
                ktobs(j,2,nwk) = joff
              enddo
            endif
          enddo
!
          call relunit (nu)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, error)
      enddo
#endif
#ifdef asm_ssh
!
!----------------------------------------------------------------------
!  read and save ssh data to direct access file
!----------------------------------------------------------------------
!
      do n=1,num_processors
        if (pn == n) then
          allocate( ltobs(jscomp:jecomp,2,nwkobs) )
          ltobs(:,:,:) = 0
!
!         write(stdout,'(/,a,i3)')
          write(6,'(a,i3)')
     & ' starting to read and sort SSH data on pn=', pn
          call flush(6)
!
          write (ssh_fle,'(a,i2.2)') 'ssha.ud.', pn
          write (opt_dasm,'(a)') ' unformatted direct words = 6'
          call getunit (iosshao, ssh_fle, opt_dasm)
          call getunit (nu, 'ssha.mom'
     &,                      'unformatted sequential rewind ieee')
          ak = 0.0
          iak = 1
          joff = 0
          do nwk = 1,nwkobs
            read (nu) obstmp(nwk), icnt
            write(stdout,'(a)') obstmp(nwk)
            write(stdout,'(a10,i3,a7,i5)') 
     &                    ' week no. ',nwk,' icnt= ',icnt
            call flush(6)
!
            if (icnt.gt.0) then
              jobs = 0
              jrow = jscomp
              do ic = 1,icnt
                read (nu) cll,rt,obs,ai,aj,err
! DEBUG
!     if (ic <= 10) write(stdout,*) cll,rt,obs,ai,aj,err
! DEBUG
                iai=ai
                iaj=aj
!
                ierr = 0
!
! Retain only data needed by each processor *********************
                if (iaj < jscomp .or. iaj > jecomp) ierr = ierr + 1
!
! Drop some data close to topography ****************************
                if (iak.gt.kma(iai  ,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai  ,iaj-1)) ierr = ierr + 1
                if (iak.gt.kma(iai+1,iaj-1)) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj  )) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj+1)) ierr = ierr + 1
                if (iak.gt.kma(iai-1,iaj-1)) ierr = ierr + 1
! End of topography check****************************************
!
                if (ierr == 0) then
                  if (iaj /= jrow) then
                    ltobs(jrow,1,nwk) = jobs
                    ltobs(jrow,2,nwk) = joff
                    if (jobs > 0) then
                      do j=1,jobs
                        write(iosshao,REC=joff+j)(aobs(i,j),i=1,6)
                      enddo
                    endif
                    joff = joff + jobs
                    do j = jrow+1,iaj-1
                      ltobs(j,2,nwk) = joff
                    enddo
                    jobs = 0
                    jrow = iaj
                  endif
                  jobs = jobs + 1
                  aobs(1,jobs) = rt + 3.5 + 7.0 * float(nwk-1)
                  aobs(2,jobs) = obs
                  aobs(3,jobs) = ai
                  aobs(4,jobs) = aj
                  aobs(5,jobs) = ak
                  aobs(6,jobs) = err
                endif
              enddo
              ltobs(jrow,1,nwk) = jobs
              ltobs(jrow,2,nwk) = joff
              if (jobs > 0) then
                do j=1,jobs
                  write(iosshao,REC=joff+j)(aobs(i,j),i=1,6)
                enddo
              endif
              joff = joff + jobs
              do j = jrow+1,jecomp
                ltobs(j,2,nwk) = joff
              enddo
            endif
!
          enddo
!
          call relunit (nu)
!         write(stdout,'(i3,a)')
!    &      nwkobs,' weeks of ssh obs. processed and saved'
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, error)
      enddo
#endif
!----------------------------------------------------------------------
!   Read time stamps and set relative times for observation periods
!   The orginal time stamps in the data represent the midpoint of
!   the data period (week, hence midpoint is Wednesday, noon).
!   Here the time stamps are changed to mark the begining and end
!   of each data period. There are "nwkobs" periods and "obstmp"
!   is dimensioned (0:nwkobs).
!   See the module "assim_mod" in assim_mod.F90 for descriptions
!   of the variables.
!----------------------------------------------------------------------
      call get_stamp (obstmp(1), iyr, mon, idy, ihr, imn, isc)
      t_time = set_date(iyr, mon, idy, ihr, imn, isc)
      o_time(1) = decrement_time(t_time, 12*60*60, 3)
      call get_date(o_time(1), iyr, mon, idy, ihr, imn, isc)
      call set_stamp (obstmp(0), iyr, mon, idy, ihr, imn, isc)
      o_time(2) = increment_time(t_time, 12*60*60, 3)
      call get_date(o_time(2), iyr, mon, idy, ihr, imn, isc)
      call set_stamp (obstmp(1), iyr, mon, idy, ihr, imn, isc)
      oprd  = 7.0
      toprd(1)  = 0.5*oprd
      toprd(2)  = toprd(1) + oprd
!
      do nwk = 3,nwkobs+1
        call get_stamp (obstmp(nwk-1), iyr, mon, idy, ihr, imn, isc)
        t_time = set_date(iyr, mon, idy, ihr, imn, isc)
        o_time(nwk) = increment_time(t_time, 12*60*60, 3)
        call get_date(o_time(nwk), iyr, mon, idy, ihr, imn, isc)
        call set_stamp (obstmp(nwk-1), iyr, mon, idy, ihr, imn, isc)
        if (nwk .le. nwkobs) then
          toprd(nwk) = toprd(nwk-1) + oprd
        endif
      enddo
!
      do nwk = 1,nwkobs
        write (stdout,'(a,i3,4a)') ' data period #',nwk
     &, ' extends from ',obstmp(nwk-1),' to ',obstmp(nwk)
      enddo
!
!----------------------------------------------------------------------
!   Find current location of model within the set of observations,
!   measured in days after obs_start_time
!----------------------------------------------------------------------
      obs_start_time = o_time(1)
      obs_end_time   = o_time(nwkobs+1)
      dayasm = position_within_obs (model_time, obs_start_time
     &,                              obs_end_time)
!----------------------------------------------------------------------
!   Set pointers to the weekly observation periods on disk.
!   These pointers are used by get_aobs, called from comobs.
!----------------------------------------------------------------------
      ndxsst = 0
      ndxtmp = 0
      ndxsal = 0
      ndxssh = 0
#ifdef asm_sst
      do nwk = 1,nwksst
        ipsstd(nwk) = nwk
      enddo
      ndxsst = 1
      call timeassim (dayasm, ndxsst, toprd, oprd, nwkobs, ipsstd,
     &                                                        nwksst)
#endif
#ifdef asm_tmp
      do nwk = 1,nwktmp
        iptmpd(nwk) = nwk
      enddo
      ndxtmp = ndxsst + 1
      call timeassim (dayasm, ndxtmp, toprd, oprd, nwkobs, iptmpd,
     &                                                        nwktmp)
#endif
#ifdef asm_sal
      do nwk = 1,nwksal
        ipsald(nwk) = nwk
      enddo
      ndxsal = ndxtmp + 1
      call timeassim (dayasm, ndxsal, toprd, oprd, nwkobs, ipsald,
     &                                                        nwksal)
#endif
#ifdef asm_ssh
      do nwk = 1,nwkssh
        ipsshd(nwk) = nwk
      enddo
      ndxssh = ndxsal + 1
      call timeassim (dayasm, ndxssh, toprd, oprd, nwkobs, ipsshd,
     &                                                        nwkssh)
#endif
!
!----------------------------------------------------------------------
!  create constant used for empirical correction for topography
!----------------------------------------------------------------------
!
      aeval=0.
      do n=1,npits-1
        aeval=1./float(n)+aeval
      enddo
      aeval=aeval/float(npits)
!
!----------------------------------------------------------------------
!     b2 = (0.5 * b * pi / 180)**2 , where b (deg) is from exp(-(x/b)**2)
!     if b = 4 deg then b2 = 1.2146e-3 .   hrscl = b
!----------------------------------------------------------------------
!
      pi = 4.0 * atan(1.0)
      b2 = (hrscl * pi / 360.0)**2
!
      write(stdout,'(/,a,e17.10,a,e17.10,a,e17.10)')
     & ' beta=',b2,' acoef=',acoef,' aeval=',aeval
!
!----------------------------------------------------------------------
!  Allocate temporary arrays to be used in subroutine lpwghts
!  The laplacian smoother is decomposed in the vertical, so each
!  processor needs these arrays allocated with the global N-S
!  dimension (jmt).
!----------------------------------------------------------------------
!
      allocate( wgns(jmt) )
      allocate( s1(imt,jmt) )
!
      do j=1,jmt
        wgns(j)=.5
      enddo
!
      call lpwghts
!
      do j=1,jmt
        if (wgns(j) .gt. 1. .or. wgns(j) .lt. 0.) then
          write(stdout,'(/,a,1pe12.4)') ' wgns nonconvergence',wgns(j)
          call flush(stdout)
          call abort()
        endif
      enddo
!
!----------------------------------------------------------------------
!   square root of weights taken because used before and after in lpsmthr
!----------------------------------------------------------------------
!
      do j=1,jmt
        wgta(j)=sqrt(acoef/wgns(j))
      enddo
!
      deallocate(wgns)
      deallocate(s1)
!
      contains
!
        subroutine lpwghts
!
!-----------------------------------------------------------------------
!   this subroutine calculates weights used in lpsmthr.
!   the approximation to [e] is made by a series of multiplications by
!   1+laplacian.  some of the initial computations of wgns here parallel
!   that of h in lpsmthr.F
!-----------------------------------------------------------------------
!
        dbsq   = 0.5*aeval*aeval
        radius = 6370.
        radsq  = radius*radius
        acon   = b2/float(npits)
!
        do j=1,jmt
          elipt(j) = 1.0
        enddo
        do j=2,jmtm1
          if (yt(j).lt.-10.0) then
            elipt(j) = 1.0 - 0.75*exp(-((yt(j)+10.0)**2/900.0))
          else if (yt(j).gt.10.0) then
            elipt(j) = 1.0 - 0.75*exp(-((yt(j)-10.0)**2/900.0))
          else
            elipt(j) = 0.25
          endif
        enddo
!
!       do j=2,jmtm1
!         elipt(j) = (dyt(j)/dyt(20))**1.81
!       enddo
!
        fyc      = 1.e10*radsq*acon
!
        do j=2,jmtm1
!
! As long as elipt = 1. do not divide by elipt...this wastes cp time
!
!       vale = fyc*cstr(j)*cstr(j)*csu(j)
!       vale = fyc*cstr(j)*cstr(j)
          vale = fyc*cstr(j)*cstr(j)*csu(j)/elipt(j)
          valn = fyc*csu(j)*cstr(j)*csu(j)*dytr(j)*dyur(j)
          vals = fyc*csu(j-1)*cstr(j)*csu(j)*dytr(j)*dyur(j-1)
!
!     print *,' j,vale,valn,vals,fyc,cstr,csu,dytr,dyur = ',
!    $j,vale,valn,vals,fyc,cstr(j),csu(j),dytr(j),dyur(j)
!
          do i=1,imt
            wew(i,j) = vale*dxtr(i)*dxur(i)
!         wno(i,j) = .5e10*csu(j)*(cstr(j)*dytr(j)+cstr(j+1)*dytr(j+1))*dyur(j)
            wno(i,j) = valn
            wso(i,j) = vals
          enddo
        enddo
!
#ifdef cyclic
        wew(1,:)   = wew(imtm1,:)
        wso(1,:)   = wso(imtm1,:)
        wno(1,:)   = wno(imtm1,:)
        wew(imt,:) = wew(2,:)
        wso(imt,:) = wso(2,:)
        wno(imt,:) = wno(2,:)
#else
        wew(1,:)   = 0.0
        wso(1,:)   = 0.0
        wno(1,:)   = 0.0
        wew(imt,:) = 0.0
        wso(imt,:) = 0.0
        wno(imt,:) = 0.0
#endif
!
        wno(:,jmt)   = 0.0
        wso(:,jmt)   = 0.0
        wew(:,jmt)   = 0.0
        wno(:,1)     = 0.0
        wso(:,1)     = 0.0
        wew(:,1)     = 0.0
!
        npid2 = npits/2
!
!-----------------------------------------------------------------------
!   set up coefficients for central point of 5 point star.
!-----------------------------------------------------------------------
!
        scra(:,1)   = 0.0
        scra(:,jmt) = 0.0
        temp(:,1)   = 0.0
        temp(:,jmt) = 0.0
!
        temp(1,:)   = 0.0
        temp(imt,:) = 0.0
!
        do j=2,jmtm1
          do i=1,imt
            scra(i,j) = 1.0-wso(i,j)-wno(i,j)-2.0*wew(i,j)
          enddo
        enddo
!
        do i=1,imt
          con           = wso(i,2)
          col           = con*con/((con+aeval)*con+dbsq)
          scra(i,2)     = scra(i,2)+con*col
          con           = wno(i,jmtm1)
          col           = con*con/((con+aeval)*con+dbsq)
          scra(i,jmtm1) = scra(i,jmtm1)+con*col
        enddo
!
        do jj=2,jmtm1
          s1(:,:) = 0.0
          s1(imt/2,jj) = 1.0
!
!-----------------------------------------------------------------------
!   begin series of multiplications by 1+laplacians. note within each
!   loop first multiplied by 1+laplacian then transpose of 1+laplacian
!   to ensure symmetry of [e] matrix
!-----------------------------------------------------------------------
!
          do l=1,npid2
            do j=2,jmtm1
              do i=2,imtm1
                temp(i,j) = scra(i,j) * s1(i,j) +
     &                      wso(i,j) * s1(i,j-1) +
     &                      wno(i,j) * s1(i,j+1) +
     &                      wew(i,j) * (s1(i+1,j) + s1(i-1,j))
              enddo
            enddo
#ifdef cyclic
            temp(imt,:) = temp(2,:)
            temp(1,:)   = temp(imtm1,:)
#endif
            do j=2,jmtm1
              do i=2,imtm1
                s1(i,j) = scra(i,j) * temp(i,j) +
     &                    wno(i,j-1) * temp(i,j-1) +
     &                    wso(i,j+1) * temp(i,j+1) +
     &                    wew(i,j) * (temp(i+1,j) + temp(i-1,j))
              enddo
            enddo
#ifdef cyclic
            s1(imt,:) = s1(2,:)
            s1(1,:)   = s1(imtm1,:)
#endif
          enddo
          wgns(jj) = s1(imt/2,jj)
        enddo
!
        end subroutine lpwghts
!
#else
      real :: setadum
#endif
      end subroutine setassim
