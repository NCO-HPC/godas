!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  GODAS_EDITPRF
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2003-07-23
!
! ABSTRACT:  Perform quality control on subsurface temperature data
!   for use in Global Ocean Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2003-07-23  David W. Behringer
! 2005-10-06  Diane C. Stokes - discard profiles in Gulf of Mexico.
!                 Recent deployment of frequently reporting profiling 
!                 buoys in this region is overwhelming analysis.  
!                 This is temporary until a more suitable and generic 
!                 fix is available.
! 2005-11-01  David W. Behringer - profiles in the Gulf of Mexico are once
!                 again allowed.  A modification to the program avePrfDly.f
!                 identifies platforms that are reporting too frequently
!                 and super-obs them.  The editing in editPrf.f has also
!                 been tuned.  A second method for fixing near surface
!                 spikes and a method for removing deep inversions have
!                 been added.
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - TEMPERATURE PROFILE DATA IN IEEE
!              - GODAS TIME_MEAN FILE IN NETCDF
!
!   OUTPUT FILES:
!     UNIT 51  - QC'D TEMPERATURE PROFILE DATA IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     UNIT 80  - SCRATCH FILE FOR GODAS PROFILE DATA
!     UNIT 81  - SCRATCH FILE FOR GODAS - OBS PROFILE DATA
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - indx, zintrp
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!       NETCDF   - nf_open, nf_inq_varid, nf_get_var_real, nf_inq_var,
!                  nf_get_vara_real, nf_close, nf_strerror
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - indx, zintrp
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!       NETCDF   - nf_open, nf_inq_varid, nf_get_var_real, nf_inq_var,
!                  nf_get_vara_real, nf_close, nf_strerror
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  10 - NETCDF ERROR
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
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
      program editPrf
!
!   "buoy"  profiles have dtyp = TR or TA or BU
!   "bathy" profiles have dtyp = BA
!   "tesac" profiles have dtyp = TE
!
      include 'netcdf.inc'
!
      integer imx,jmx,kmx,tmx
      parameter (imx=360,jmx=200,kmx=40,tmx=1)
      integer ii,jj,kk,ncid,status
      integer ndims, nvars, natts, idunlm
      integer dsiz(20), vtype, nvdm, vdm(10), nvatts
      integer start(4), count(4)
      character*20 dname, vname
!
      real xt(imx), yt(jmx), zt(kmx), temp(kmx)
!
      character csign*8,ukey*1,sid*2,dtyp*2,qkey*1,msign*8,dsign*8
      real, allocatable, dimension(:) :: pt, pz, tti, dt, rfA, rfB
!     real pt(10), pz(10), tti(10), dt(10)
      logical noerr, despk, spke, mooring
!
      call w3tagb('GODAS_EDITPRF',2005,0319,0319,'NP23')
!
! set some limits and some parameters
!
      nzerr = 0
!
      z1max = 25.0
      nz1err = 0
!
      gpmax = 500.0
      nlgp = 0
!
      t1smax = 0.01
      nt1spk = 0
!
      bmxl = 0.25
      dmxl = 0.11
      fmxl = 0.2
      nt2spk = 0
!
      tbsmax = 0.06
      zbsmin = 500.0
      ntbspk = 0
!
      nspkpss = 2
!     nsmax = 5
      nsmax = 500
      zsAmin = 40.0
      rfA0 = 10.0
!     sdfA = 3.0
      sdfA = 4.0
      zsBmin = 250.0
      rfB0 = 10.0
!     sdfB = 2.0
      sdfB = 3.0
      ntfspk = 0
!
      gdfm = 4.0
      mtomax = 4.0
!
      gdf1 = 4.0
      gdf2 = 6.0
      gdf3 = 4.0
      glat1 = 45.0
      glat2 = 30.0
      ngdiff = 0
      gdffn = 8.0
!
!     zsd = 250.0
      zsd = 500.0
      dtomax = 2.0
      ndoff = 0
!
      nnoise = 0
      nobnds = 0
!
      ninv = 0
      sltinv = -35.0
      nltinv = 35.0
!
! open existing netCDF ocean file
!
      status = nf_open('time_mean.nc', 0, ncid)
      if (status .ne. NF_NOERR) go to 100
!
! get xt
!
      status = nf_inq_varid(ncid, 'xt_i', nid)
      if (status .ne. NF_NOERR) go to 100
      status = nf_get_var_real(ncid, nid, xt)
      if (status .ne. NF_NOERR) go to 100
!
! get yt
!
      status = nf_inq_varid(ncid, 'yt_j', nid)
      if (status .ne. NF_NOERR) go to 100
      status = nf_get_var_real(ncid, nid, yt)
      if (status .ne. NF_NOERR) go to 100
!
! get zt
!
      status = nf_inq_varid(ncid, 'zt_k', nid)
      if (status .ne. NF_NOERR) go to 100
      status = nf_get_var_real(ncid, nid, zt)
      if (status .ne. NF_NOERR) go to 100
!
! set fixed start / count
!
      count(1) = 1
      count(2) = 1
      start(3) = 1
      count(3) = kmx
      start(4) = 1
      count(4) = 1
      msign = 'modelPrf'
!
! open profile file and file for edited profiles
!
      open (11, form='unformatted', status='old', access='sequential', &
            &  err=110)
      open (51, form='unformatted', access='sequential')
!
! open file for model profiles and file for obs - model difference
!
      open (80, status='scratch', form='unformatted', &
              & access='sequential')
      open (81, status='scratch', form='unformatted', &
              & access='sequential')
!
! begin loop on profile file
!
      nprf = 0
      npmx = 0
      do while (.true.)
        read (11, end=10, err=120) iyear,idate,csign, &
                               &  sid,dtyp,qkey,yp,xp,np
        nprf = nprf + 1
        if (np .gt. npmx) npmx = np
      end do
  10  continue
!
      allocate(pt(npmx))
      allocate(pz(npmx))
      allocate(tti(npmx))
      allocate(dt(npmx))
      allocate(rfA(npmx))
      allocate(rfB(npmx))
!
      rewind 11
!
      npout = 0
      do iprf=1,nprf
        read (11, err=120) iyear,idate,csign,sid,dtyp,qkey,yp,xp, &
                           &  np,(pz(k),pt(k),k=1,np)
        npi = np

        noerr = .true.
        if (dtyp .eq. 'TA' .or. dtyp .eq. 'TR') then
          mooring = .true.
        else
          mooring = .false.
        end if
!
! treat moorings (largely TAO and Triton) separately
!
        if (mooring) then
!
! check for non-monotonic depths
!
          if (noerr) then
            do k=2,np
              if (pz(k) .le. pz(k-1)) then
                noerr = .false.
                nzerr = nzerr + 1
                exit
              end if
            end do
          end if
!
! get 'temp' profile from ocean file
!
          ii = indx(xp,xt,imx,.true.)
          jj = indx(yp,yt,jmx,.false.)
!
          if (noerr .and. ii .gt. 0 .and. jj .gt. 0) then
            status = nf_inq_varid(ncid, 'temp', nid)
            if (status .ne. NF_NOERR) go to 100
            status = nf_inq_var(ncid,nid,vname,vtype,nvdm,vdm,nvatts)
            if (status .ne. NF_NOERR) go to 100
            start(1) = ii
            start(2) = jj
            status = nf_get_vara_real(ncid, nid, start, count, temp)
            if (status .ne. NF_NOERR) go to 100
!
            kk = 0
            do k=1,kmx
              if (temp(k) .gt. -990. .and. temp(k) .lt. 990.) kk = k
            end do
!
            if (kk .gt. 0) then
              call zintrp(temp,zt,kk,tti,pz,np,npi)
              rms = 0.0
              dtoff = 0.0
              do k=1,npi
                dt(k) = pt(k) - tti(k)
                rms = rms + dt(k)**2
                if (abs(dt(k)) .gt. dtoff) dtoff = abs(dt(k))
              end do
              rms = sqrt(rms/npi)
!
! check for large offsets
!
              if (dtoff .gt. mtomax) then
                noerr = .false.
                ngdiff = ngdiff + 1
              end if
!
! check for large rms difference
!
              if (rms .gt. gdfm) then
                noerr = .false.
                ngdiff = ngdiff + 1
              end if
!
            else
              noerr = .false.
              nobnds = nobnds + 1
            end if
!
          else
            noerr = .false.
            nobnds = nobnds + 1
          end if
!
        else
!
! check for non-monotonic depths
!
          if (noerr) then
            do k=2,np
              if (pz(k) .le. pz(k-1)) then
                noerr = .false.
                nzerr = nzerr + 1
                exit
              end if
            end do
          end if
!
! check for deep start
!
          if (noerr) then
            if (pz(1) .gt. z1max) then
              noerr = .false.
              nz1err = nz1err + 1
            end if
          end if
!
! check for large gaps
!
          if (noerr) then
            do k=2,np
              if (pz(k) - pz(k-1) .gt. gpmax) then
                np = k - 1
                nlgp = nlgp + 1
                exit
              end if
            end do
          end if
!
! fix surface spike, type 1
!
          if (noerr .and. np .gt. 1) then
            if (abs((pt(1)-pt(2))/(pz(2)-pz(1))) .gt. t1smax) then
              pt(1) = pt(2)
              nt1spk = nt1spk + 1
            end if
          end if
!
! fix bottom spike
!
          if (noerr) then
            spke = .false.
            do while (pz(np) .gt. zbsmin) 
              if (abs((pt(np-1)-pt(np))/(pz(np)-pz(np-1)))  &
                                        &     .gt. tbsmax) then
                np = np - 1
                spke = .true.
              else
                exit
              end if
            end do
            if (spke) ntbspk = ntbspk + 1
          end if
!
! get 'temp' profile from ocean file
!
          ii = indx(xp,xt,imx,.true.)
          jj = indx(yp,yt,jmx,.false.)
!
          if (noerr .and. ii .gt. 0 .and. jj .gt. 0) then
            status = nf_inq_varid(ncid, 'temp', nid)
            if (status .ne. NF_NOERR) go to 100
            status = nf_inq_var(ncid,nid,vname,vtype,nvdm,vdm,nvatts)
            if (status .ne. NF_NOERR) go to 100
            start(1) = ii
            start(2) = jj
            status = nf_get_vara_real(ncid, nid, start, count, temp)
            if (status .ne. NF_NOERR) go to 100
!
            kk = 0
            do k=1,kmx
              if (temp(k) .gt. -990. .and. temp(k) .lt. 990.) kk = k
            end do
!
            if (kk .gt. 0) then
              call zintrp(temp,zt,kk,tti,pz,np,npi)
!
! fix surface spike, type 2
!
              kmxl = npi
              do k=2,npi
                if (tti(1)-tti(k) .ge. bmxl) then
                  kmxl = k
                  exit
                end if
              end do
              if (kmxl .gt. 5 .and. kmxl .lt. npi) then
                do k=1,kmxl
                  mdt = 0
                  do ks=1,kmxl
                    if (abs(pt(k)-pt(ks)) .le. dmxl) mdt = mdt + 1
                  end do
                  if (float(mdt)/float(kmxl) .gt. fmxl) then
                    ks = k
                    exit
                  end if
                  ks = 0
                end do
!
                if (ks .gt. 1) then
                  nt2spk = nt2spk + 1
                  do k=1,ks
                    pt(k) = pt(ks)
                  end do
                end if
!
              end if
!
              rms = 0.0
              rmss = 0.0
              rmsd = 0.0
              knts = 0
              kntd = 0
              dtoff = 0.0
              do k=1,npi
                dt(k) = pt(k) - tti(k)
                rms = rms + dt(k)**2
                if (k .eq. 1) then
                  dtdz = abs((tti(2)-tti(1))/(pz(2)-pz(1)))
                else if (k .eq. npi) then
                  dtdz = abs((tti(npi)-tti(npi-1))/(pz(npi)-pz(npi-1)))
                else
                  dtdz = abs((tti(k+1)-tti(k-1))/(pz(k+1)-pz(k-1)))
                endif
                rfA(k) = rfA0 * dtdz + sdfA
                rfB(k) = rfB0 * dtdz + sdfB
                if (pz(k) .le. zsd) then
                  rmss = rmss + dt(k)**2
                  knts = knts + 1
                else
                  rmsd = rmsds + dt(k)**2
                  kntd = kntd + 1
                end if
                if (pz(k) .gt. zsd) then
                  if (abs(dt(k)) .gt. dtoff) dtoff = abs(dt(k))
                end if
              end do
!
              rms = sqrt(rms/npi)
              if (knts .gt. 0) rmss = sqrt(rmss/knts)
              if (kntd .gt. 0) rmsd = sqrt(rmsd/kntd)
!
              if (yp .gt. glat1) then
                if (rms .gt. gdf1) then
                  noerr = .false.
                  ngdiff = ngdiff + 1
                end if
              else if (yp .gt. glat2) then
                if (rms .gt. gdf2) then
                  noerr = .false.
                  ngdiff = ngdiff + 1
                end if
              else
                if (rms .gt. gdf3) then
                  noerr = .false.
                  ngdiff = ngdiff + 1
                end if
              end if
!
              if (noerr) then
                if (abs(yp) .ge. glat1) then
                  if (rmsd .gt. rmss) then
                    noerr = .false.
                    ndoff = ndoff + 1
                  end if
                else
                  if (rmsd .gt. 2.0*rmss) then
                    noerr = .false.
                    ndoff = ndoff + 1
                  end if
                end if
              end if
!
              if (noerr) then
                if (dtoff .ge. dtomax) then
                  noerr = .false.
                  ndoff = ndoff + 1
                end if
              end if
!
            else
              noerr = .false.
              nobnds = nobnds + 1
            end if
!
          else
            noerr = .false.
            nobnds = nobnds + 1
          end if
!
! fix spikes - first check
!
          if (noerr) then
            ns = 0
            do k=1,npi
              if (pz(k) .gt. zsAmin &
                      &   .and. abs(dt(k)) .ge. rfA(k)*rms) then
                ns = ns + 1
              end if
            end do
            if (ns .gt. nsmax) then
              noerr = .false.
              nnoise = nnoise + 1
            else
              despk = .false.
              do n=1,nspkpss
                do k=1,npi
                  if (pz(k) .gt. zsAmin &
                         &  .and. abs(dt(k)) .ge. rfA(k)*rms) then
                    if (k .eq. 1) then
                      kp = 2
                      do while (abs(dt(kp)) .ge. rfA(kp)*rms &
                                         &  .and. kp .lt. npi)
                        kp = kp + 1
                      end do
                      if (kp-k .gt. 3) then
                        noerr = .false.
                        nnoise = nnoise + 1
                        exit
                      end if
                      do kk=k,kp-1
                        pt(kk) = pt(kp)
                        dt(kk) = pt(kk) - tti(kk)
                      end do
                    else if (k .eq. npi) then
                      km = npi - 1
                      do while (abs(dt(km)) .ge. rfA(km)*rms &
                                            &  .and. km .gt. 1)
                        km = km - 1
                      end do
                      if (k-km .gt. 3) then
                        noerr = .false.
                        nnoise = nnoise + 1
                        exit
                      end if
                      npi = km
                      np = km
                    else
                      km = k - 1
                      do while (abs(dt(km)) .ge. rfA(km)*rms &
                                             &  .and. km .gt. 1)
                        km = km - 1
                      end do
                      kp = k + 1
                      do while (abs(dt(kp)) .ge. rfA(kp)*rms &
                                             &  .and. kp .lt. npi)
                        kp = kp + 1
                        end do
                      if (kp-km .gt. 4) then
                        noerr = .false.
                        nnoise = nnoise + 1
                        exit
                      end if
                      if (kp .eq. npi &
                        & .and. abs(dt(kp)) .ge. rfA(kp)*rms) then
                        npi = km
                        ntbspk = ntbspk + 1
                        exit
                      endif
                      do kk=km+1,kp-1
                        dzm = (pz(kk) - pz(km)) / (pz(kp) - pz(km))
                        dzp = (pz(kp) - pz(kk)) / (pz(kp) - pz(km))
                        pt(kk) = pt(km)*dzp + pt(kp)*dzm
                        dt(kk) = pt(kk) - tti(kk)
                      end do
                    end if
                    despk = .true.
                  end if
                end do
!
                rms = 0.0
                do k=1,npi
                  dt(k) = pt(k) - tti(k)
                  rms = rms + dt(k)**2
                end do
                rms = sqrt(rms/npi)
!
              end do
              if (despk) ntfspk = ntfspk + 1
            end if
          end if
!
! fix spikes - second check
!
          if (noerr .and. pz(npi) .gt. zsBmin) then
            np0 = 1
            do k=1,npi
              if (pz(k) .ge. zsBmin) then
                np0 = k
                exit
              end if
            end do
!
            rms = 0.0
            do k=np0,npi
              dt(k) = pt(k) - tti(k)
              rms = rms + dt(k)**2
            end do
            rms = sqrt(rms/(npi-np0+1))
!
            ns = 0
            do k=np0,npi
              if (abs(dt(k)) .ge. rfB(k)*rms) then
                ns = ns + 1
              end if
            end do
            if (ns .gt. nsmax) then
              noerr = .false.
              nnoise = nnoise + 1
            else
              despk = .false.
              do n=1,nspkpss
                do k=np0,npi
                  if (abs(dt(k)) .ge. rfB(k)*rms) then
                    if (k .eq. npi) then
                      km = npi - 1
                      do while (abs(dt(km)) .ge. rfB(km)*rms &
                                            &  .and. km .gt. 1)
                        km = km - 1
                      end do
                      if (k-km .gt. 3) then
                        noerr = .false.
                        nnoise = nnoise + 1
                        exit
                      end if
                      npi = km
                      np = km
                    else
                      km = k - 1
                      do while (abs(dt(km)) .ge. rfB(km)*rms &
                                             &  .and. km .gt. 1)
                        km = km - 1
                      end do
                      kp = k + 1
                      do while (abs(dt(kp)) .ge. rfB(kp)*rms &
                                             &  .and. kp .lt. npi)
                        kp = kp + 1
                      end do
                      if (kp-km .gt. 4) then
                        noerr = .false.
                        nnoise = nnoise + 1
                        exit
                      end if
                      if (kp .eq. npi &
                        & .and. abs(dt(kp)) .ge. rfB(kp)*rms) then
                        npi = km
                        ntbspk = ntbspk + 1
                        exit
                      endif
                      do kk=km+1,kp-1
                        dzm = (pz(kk) - pz(km)) / (pz(kp) - pz(km))
                        dzp = (pz(kp) - pz(kk)) / (pz(kp) - pz(km))
                        pt(kk) = pt(km)*dzp + pt(kp)*dzm
                        dt(kk) = pt(kk) - tti(kk)
                      end do
                    end if
                    despk = .true.
                  end if
                end do
!
                rms = 0.0
                do k=np0,npi
                  dt(k) = pt(k) - tti(k)
                  rms = rms + dt(k)**2
                end do
                rms = sqrt(rms/(npi-np0+1))
!
              end do
              if (despk) ntfspk = ntfspk + 1
            end if
          end if
!
! do final gross difference check
!
          if (noerr) then
            rms = 0.0
            do k=1,npi
              dt(k) = abs(pt(k) - tti(k))
              rms = rms + dt(k)**2
            end do
            rms = sqrt(rms/npi)
!
            do k=1,npi
              if (dt(k) .ge. gdffn) then
                noerr = .false.
                ngdiff = ngdiff + 1
                exit
              endif
            end do
          end if
!
! do final check for inversions in mid-latitudes
!
          if (noerr) then
            if (yp .ge. sltinv .and. yp .le. nltinv) then
              komxl = -1
              do k=kmxl,npi
                if (pt(k) .lt. pt(kmxl)) then
                  komxl = k
                  exit
                endif
              end do
!
              if (komxl .gt. 0) then
                npt = np
                do k=komxl+1,np
                  if (pt(k) .gt. pt(k-1)) then
                    npt = k-1
                    ninv = ninv + 1
                    exit
                  end if
                end do
                np = npt
              end if
            end if
          end if
!
        end if
!
! write edited, model and difference profiles
!
        if (noerr) then
          npout = npout + 1
!
          write (51) iyear,idate,csign,sid,dtyp,qkey,yp,xp, &
                           &  np,(pz(k),pt(k),k=1,np)
!
          write (80) iyear,idate,msign,sid,dtyp,qkey,yt(jj), &
                      & xt(ii),npi,(pz(k),tti(k),k=1,npi)
!
          write(dsign,'(f6.1,2x)') rms
          write (81) iyear,idate,dsign,sid,dtyp,qkey,yt(jj), &
                       & xt(ii),npi,(pz(k),dt(k),k=1,npi)
!
        end if
!
      end do
!
      status = nf_close(ncid)
      close (11)
      close (51)
      close (80)
      close (81)
!
      write(6,'(a)') 'Results of editing profiles'
      write(6,'(i5,a)') nprf, ' profiles read'
      write(6,'(i5,a)') npout, ' profiles retained'
      write(6,'(i5,a)') nzerr, ' deleted for non-mono depths'
      write(6,'(i5,a,f3.0)') nz1err, ' deleted for deep z1 > ', z1max
      write(6,'(i5,a)') nobnds, ' deleted as out of bounds'
      write(6,'(i5,a)') ngdiff, ' deleted for gross differences'
      write(6,'(i5,a)') ndoff, ' deleted for deep offset'
      write(6,'(i5,a)') nnoise, ' deleted as noisy'
      write(6,'(i5,a,f4.0)') nlgp, ' truncated for large gap > ', gpmax
      write(6,'(i5,a)') nt1spk, ' type 1 surface spikes fixed'
      write(6,'(i5,a)') nt2spk, ' type 2 surface spikes fixed'
      write(6,'(i5,a)') ntbspk, ' bottom spikes fixed'
      write(6,'(i5,a)') ntfspk, ' profiles had other spikes fixed'
      write(6,'(i5,a)') ninv, ' profiles truncated below inversion'
      write(6,*)
!
      call w3tage('GODAS_EDITPRF')
      call errexit(0)
      call exit(0)
!
  100 continue
      write(6,'(a)') 'Error reading time-mean netCDF file'
      write(6,'(a)') nf_strerror(status)
      call w3tage('GODAS_EDITPRF')
      call errexit(10)
!     call exit(10)
!
  110 write(6,'(a)') 'Error opening profile file on unit 11'
      call w3tage('GODAS_EDITPRF')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading profile file on unit 11'
      call w3tage('GODAS_EDITPRF')
      call errexit(12)
!     call exit(12)
!
      end program editPrf

! -------------------------------------------------------------------

      integer function indx(p, pt, nmx, cyclic)
      real p, pt(*)
      integer nmx
      logical cyclic
!
      if (p .ge. pt(1) .and. p .le. pt(nmx)) then
        do n=2,nmx
          if (p .ge. pt(n-1) .and. p .le. pt(n)) then
            if (p-pt(n-1) .lt.  pt(n)-p) then
              indx = n-1
            else
              indx = n
            end if
          end if
        end do
      else
        if (cyclic) then
          indx = 1
        else
          indx = -1
        end if
      end if
!
      end function indx
!
! -------------------------------------------------------------------
!
      subroutine zintrp(ta,za,na,tb,zb,nb,nbi)
      integer na, nb
      real ta(*), za(*), tb(*), zb(*)
!
      nbi = -1
      do k=1,nb
        if (zb(k) .le. za(1)) then
          tb(k) = ta(1)
          nbi = k
        else if (zb(k) .gt. za(na)) then
          exit
        else
          do ka=2,na
            if (zb(k) .gt. za(ka-1) .and. zb(k) .le. za(ka)) then
              dzm = (zb(k) - za(ka-1)) / (za(ka) - za(ka-1))
              dzp = (za(ka) - zb(k)) / (za(ka) - za(ka-1))
              tb(k) = ta(ka-1)*dzp + ta(ka)*dzm
              nbi = k
              exit
            end if
          end do
        end if
      end do
!
      end subroutine zintrp
