  program mkEvNc
!
  implicit none
!
  include 'netcdf.inc'

! void getTz(), writeFile(), usage();
! void nrmFile(), smthFile(), fltrFile(), cielFile(), sfcMnMx();

! mode_t mode = 0644;

  integer :: imx, jmx, kmx
  integer :: kass, iav, jav, kav
  real(kind=4), allocatable, dimension(:) :: xt, yt, zt, dz
  real(kind=4), allocatable, dimension(:) :: zm, dzm, wtk, wtkp
  real(kind=4), allocatable, dimension(:) :: ab, ap, az
  real(kind=4), allocatable, dimension(:,:) :: buf
  real(kind=4), allocatable, dimension(:,:,:) :: a, b, ev
  real(kind=4) :: bmx
  integer :: nSqrt, Sfc
  integer :: Nrml, ltS, ltN, lgW, lgE
  integer :: year, month, day
  integer, allocatable, dimension(:,:) :: mask
  real(kind=4) :: VsMn, VsMx, gscl, LsTpr, LnTpr, period, gciel
  character(len=32) :: stamp
  character(len=121) :: str
  character(len=91) :: fileIn, fileOut, ulTitle, urDate
  character(len=3) :: bsn
  integer :: i, j, k, ka, n, narg
  logical :: SdZ, zTpr, prdc, flg

  real(kind=4), parameter :: spv = 999.999
! rpd = PI / 180.0
  real(kind=4), parameter :: rpd = 0.017453293

! netCDF

  integer :: ncid, xtVid, ytVid, ztVid, tVid
  integer :: ndims, nvars, ngatts, xdimid, nspdms
  integer, dimension(4) :: vstart, vcount
  real(kind=4) :: msV
  character(len=121) :: attnm, vname, dname
  integer :: nfstat, idunlm
  integer :: dimid, varid, pVid, len

! set some defaults
  kass = 30
  iav = 11
  jav = 5
  kav = 5
  SdZ = .true.
  nSqrt = 2
  zTpr = .false.
  prdc = .true.
  Sfc = 1
  Nrml = 1
  ltS = -2
  ltN = 2
  lgW = 130
  lgE = 250
!
  VsMn = 0.30
  VsMx = 0.80
  gscl = 1.0
  LsTpr = -90.0
  LnTpr = 90.0
  gciel = 1.0
!
  fileIn = 'EMPTY'
  fileOut = 'EMPTY'
  narg = iargc()
  n = 1
  do while (n .le. narg)
    call getarg(n,str)
    if (str(1:2) .eq. '-f') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      fileIn = str
    else if (str(1:2) .eq. '-o') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      fileOut = str
    else if (str(1:2) .eq. '-k') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i2)') kass
    else if (str(1:4) .eq. '-LtN') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i3)') ltS
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i3)') ltN
      Nrml = 1
    else if (str(1:4) .eq. '-LgN') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i4)') lgW
      if (lgW .lt. 0) lgW = lgW + 360;
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i4)') lgE
      if (lgE .lt. 0) lgE = lgE + 360;
      Nrml = 1
    else if (str(1:3) .eq. '-gS') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(f4.1)') gscl
    else if (str(1:3) .eq. '-gC') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(f4.1)') gciel
    else if (str(1:3) .eq. '-Av') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i3)') iav
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i3)') jav
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i3)') kav
    else if (str(1:5) .eq. '-Vsfc') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(f5.2)') VsMn
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(f5.2)') VsMx
      Sfc = 1
    else if (str(1:4) .eq. '-SdZ') then
      SdZ = .true.
    else if (str(1:5) .eq. '-SqRt') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(i1)') nSqrt
    else if (str(1:5) .eq. '-zTpr') then
      zTpr = .true.
    else if (str(1:5) .eq. '-yTpr') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(f4.1)') LsTpr
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(f4.1)') LnTpr
    else if (str(1:2) .eq. '-h') then
      call usage
      call exit(99)
    else
      call usage
      call exit(99)
    endif
    n = n + 1
  enddo

  if (fileIn .eq. 'EMPTY') then
    write(6,'(a)') 'Error: An input file must be given.'
    call usage
    call exit(99)
  endif
  if (fileOut .eq. 'EMPTY') then
    write(6,'(a)') 'Error: An output file must be given.'
    call usage
    call exit(99)
  endif
  if (.not.Nrml) then
    gciel = -99.0
  endif

  call getTz

 call fltrFile
!
 call nrmFile
!
 call cielFile
!
 call smthFile
!
 call nrmFile
!
 call sfcMnMx
!
  call writeFile
!
  contains
!
! ===================================================================

  subroutine getTz

!  Computing Tz from fileIn

  nfstat = nf_open(fileIn, NF_NOWRITE, ncid)
  if (nfstat .ne. NF_NOERR) then
    write(6,'(a,a)') 'Could not open ',trim(fileIn)
    call exit(99)
  endif

  nfstat = nf_inq(ncid, ndims, nvars, ngatts, idunlm)
  if (nfstat .ne. NF_NOERR) then
    write(6,'(a,a)') 'Could not inquire about ',trim(fileIn)
    call exit(99)
  endif

  do i=1,ngatts
    nfstat = nf_inq_attname(ncid, NF_GLOBAL,i,attnm)
    if (nfstat .ne. NF_NOERR) then
      write(6,'(a,a)') 'Could not read globsl attname in ',trim(fileIn)
      call exit(99)
    endif
    if (attnm .eq. 'filename') then
      nfstat = nf_get_att_text(ncid, NF_GLOBAL, 'filename', str)
      if (nfstat .ne. NF_NOERR) then
        write(6,'(a,a)') 'Could not read global attname ','filename'
        call exit(99)
      endif
      read(str(11:22),'(i6,1x,i2,1x,i2)') year, month, day
    endif
  enddo

  do varid=1,nvars
    nfstat = nf_inq_varname(ncid, varid, vname)
    if (nfstat .ne. NF_NOERR) then
      write(6,'(a,a)') 'Could not read varname in ',trim(fileIn)
      call exit(99)
    endif
    if (vname .eq. 'period') then
      vstart(1) = 1
      vcount(1) = 1
      nfstat = nf_get_vara_real(ncid, varid, vstart, vcount, period)
      if (nfstat .ne. NF_NOERR) then
        write(6,'(a,a)') 'Could not read variable ',trim(vname)
        call exit(99)
      endif
    endif
  enddo

  do dimid=1,ndims
    nfstat = nf_inq_dim(ncid, dimid, dname, len)
    if (nfstat .ne. NF_NOERR) then
      write(6,'(a,a)') 'Could not read dimensions in ',trim(fileIn)
      call exit(99)
    endif
    if (dname .eq. 'xt_i') then
      imx = len
    else if (dname .eq. 'yt_j') then
      jmx = len
    else if (dname .eq. 'zt_k') then
      kmx = len
    endif
  enddo

  allocate (xt(imx))
  allocate (yt(jmx))
  allocate (zt(kmx))
  allocate (mask(imx,jmx))

  do varid=1,nvars
    nfstat = nf_inq_varname(ncid, varid, vname)
    if (nfstat .ne. NF_NOERR) then
      write(6,'(a,a)') 'Could not read varname in ',trim(fileIn)
      call exit(99)
    endif
    if (vname .eq. 'xt_i') then
      xtVid = varid
    else if (vname .eq. 'yt_j') then
      ytVid = varid
    else if (vname .eq. 'zt_k') then
      ztVid = varid
    else if (vname .eq. 'temp') then
      tVid = varid
    endif
  enddo

  vstart(1) = 1
  vcount(1) = imx
  nfstat = nf_get_vara_real(ncid, xtVid, vstart, vcount, xt)
  if (nfstat .ne. NF_NOERR) then
    write(6,'(a,a)') 'Could not read variable ','xt_i'
    call exit(99)
  endif
  vcount(1) = jmx
  nfstat = nf_get_vara_real(ncid, ytVid, vstart, vcount, yt)
  if (nfstat .ne. NF_NOERR) then
    write(6,'(a,a)') 'Could not read variable ','yt_j'
    call exit(99)
  endif
  vcount(1) = kmx
  nfstat = nf_get_vara_real(ncid, ztVid, vstart, vcount, zt)
  if (nfstat .ne. NF_NOERR) then
    write(6,'(a,a)') 'Could not read variable ','zt_k'
    call exit(99)
  endif

  allocate (zm(kmx))
  allocate (dzm(kmx))
  allocate (wtk(kmx))
  allocate (wtkp(kmx))
  zm(1) = 0.0
  do k=2,kmx
    zm(k) = 0.5*(zt(k-1) + zt(k))
    dzm(k) = zt(k) - zt(k-1)
  enddo
  dzm(1) = dzm(2)
  do k=1,kmx-1
    wtkp(k) =  (zm(k+1) - zt(k)) / (zm(k+1) - zm(k))
    wtk(k) =  (zt(k) - zm(k)) / (zm(k+1) - zm(k))
  enddo
  wtkp(kmx) = wtkp(kmx-1)
  wtk(kmx) = wtk(kmx-1)


  nfstat = nf_get_att_real(ncid, tVid, 'missing_value', msV)
  if (nfstat .ne. NF_NOERR) then
    write(6,'(a,a)') 'Could not variable ','temp missing_value'
    call exit(99)
  endif

  allocate(b(imx,jmx,kmx))
  allocate(ab(kmx))
  allocate(ap(kmx))
  allocate(az(kmx))
  allocate(buf(imx,kmx))

  vstart(1) = 1
  vcount(1) = imx
  vstart(2) = 1
  vcount(2) = 1
  vstart(3) = 1
  vcount(3) = kmx
  vstart(4) = 1
  vcount(4) = 1

  do j=1,jmx
    vstart(2) = j
    nfstat = nf_get_vara_real(ncid, tVid, vstart, vcount, buf)
    if (nfstat .ne. NF_NOERR) then
      write(6,'(a,a,i5)') 'Could not read variable ','temp',j
      call exit(99)
    endif
    do i=1,imx
      if (buf(i,1) .ne. msV) then
        ka = 1
        do while (ka+1 .lt. kmx .and. buf(i,ka+1) .ne. msV)
          ka = ka + 1
        enddo
        mask(i,j) = ka
        do k=1,ka
          ab(k) = buf(i,k)
        enddo
        do k=2,ka
          ap(k) = (ab(k-1) - ab(k)) / dzm(k)
        enddo
        ap(1) = ap(2)
        do k=2,ka-1
          az(k) = ap(k) * wtkp(k)  +  ap(k+1) * wtk(k)
        enddo
        az(1) = az(2);
        az(ka) = az(ka-1)
        do k=1,kmx
          if (k .le. ka) then
            b(i,j,k) = az(k)
          else
            b(i,j,k) = spv
          endif
        enddo
      else
        do k=1,kmx
          b(i,j,k) = spv
        enddo
        mask(i,j) = 0
      endif
    enddo
  enddo
  nfstat = nf_close(ncid)
  if (nfstat .ne. NF_NOERR) then
    write(6,'(a,a)') 'Error closing ',trim(fileIn)
    call exit(99)
  endif

  bmx = 0.0
  do j=1,jmx
    do i=1,imx
      do k=1,mask(i,j)
        if (b(i,j,k) .gt. bmx) bmx =  b(i,j,k)
      enddo
    enddo
  enddo
  if (bmx .gt. 0.0) then
    do j=1,jmx
      do i=1,imx
        do k=1,mask(i,j)
          if (b(i,j,k) < 0.0) then
            b(i,j,k) = 0.0
          else
            b(i,j,k) = b(i,j,k) / bmx
          endif
        enddo
      enddo
    enddo
  endif

  allocate(dz(kmx))
  dz(1) = 2.0* zt(1)
  do k=2,kmx
    dz(k) = 2.0*(zt(k) - zt(k-1)) - dz(k-1)
  enddo

  do j=1,jmx
    do i=1,imx
      if (mask(i,j) .gt. 1) then
        b(i,j,1) = b(i,j,2)
      endif
    enddo
  enddo
!
  end subroutine getTz

! ===================================================================
!
  subroutine nrmFile
!
  real(kind=4) lats, latn, lngw, lnge, fctr
!
  bmx = 0.0
  lats = real(ltS,4)
  latn = real(ltN,4)
  lngw = real(lgW,4)
  lnge = real(lgE,4)
  do j=1,jmx
    if (yt(j) .gt. latn) exit
    if (yt(j) .ge. lats) then
      do i=1,imx
        if (xt(i) .gt. lnge) exit
        if (xt(i) .ge. lngw) then
          do k=1,mask(i,j)
            if (b(i,j,k) .gt. bmx) bmx = b(i,j,k)
          enddo
        endif
      enddo
    endif
  enddo
  fctr = gscl / bmx
  do j=1,jmx
    do i=1,imx
      do k=1,mask(i,j)
        b(i,j,k) = b(i,j,k) * fctr
      enddo
    enddo
  enddo
!
  end subroutine nrmFile
!
! ===================================================================
!
  subroutine cielFile
!
  if (gciel .lt. 0.0) return
!
  do j=1,jmx
    do i=1,imx
      do k=1,mask(i,j)
        if (b(i,j,k) .gt. gciel) b(i,j,k) = gciel
      enddo
    enddo
  enddo
!
  end subroutine cielFile
!
! ===================================================================
!
  subroutine smthFile
!
  integer :: ii, i3, im, ip, jj, jm, jp, kk, km, kp, cnt;

  allocate(a(imx,jmx,kmx))
  if (iav .ne. 0 .and. jav .ne. 0) then
    iav = iav/2
    jav = jav/2
    do k=1,kmx
      do j=1,jmx
        jm = j-jav
        if (jm .lt. 1) jm = 1
        jp = j+jav
        if (jp .gt. jmx) jp = jmx
        do i=1,imx
          if (k .le. mask(i,j)) then
            if (prdc) then
              im = i-iav
              ip = i+iav
            else
              im = i-iav
              if (im .lt. 1) im = 1
              ip = i+iav
              if (ip .gt. imx) ip = imx
            endif
            cnt = 0
            a(i,j,k) = 0.0
            do jj=jm,jp
              do i3=im,ip
                ii = i3;
                if (ii .lt. 1) ii = ii + imx;
                if (ii .gt. imx) ii = ii - imx;
                if (k .le. mask(ii,jj)) then
                  a(i,j,k) = a(i,j,k) + b(ii,jj,k)
                  cnt = cnt + 1
                endif
              enddo
            enddo
            if (cnt .gt. 0) a(i,j,k) = a(i,j,k) / real(cnt,4)
          else
            a(i,j,k) = spv
          endif
        enddo
      enddo
    enddo
    iav = 2*iav + 1
    jav = 2*jav + 1
  else
    do k=1,kmx
      do j=1,jmx
        do i=1,imx
          a(i,j,k) = b(i,j,k)
        enddo
      enddo
    enddo
  endif

  if (kav .gt. 1) then
    kav = kav/2;
    do j=1,jmx
      do i=1,imx
        do k=1,kmx
          if (k .le. mask(i,j)) then
            km = k-kav
            if (km .lt. 1) km = 1
            kp = k+kav
            if (kp .gt. kmx) kp = kmx
            cnt = 0
            b(i,j,k) = 0.0
            do kk=km,kp
              if (kk .le. mask(i,j)) then
                 b(i,j,k) = b(i,j,k) + a(i,j,kk)
                cnt = cnt + 1
              endif
            enddo
            if (cnt .gt. 0) b(i,j,k) =  b(i,j,k) / real(cnt,4)
          else
            b(i,j,k) = spv
          endif
        enddo
      enddo
    enddo
    kav = 2*kav + 1
  else
    do k=1,kmx
      do j=1,jmx
        do i=1,imx
           b(i,j,k) = a(i,j,k)
        enddo
      enddo
    enddo
  endif
!
  end subroutine smthFile

! ===================================================================
!
  subroutine fltrFile
!
  integer :: kt
  real(kind=4) :: zmx

  if (SdZ) then
    do j=1,jmx
      do i=1,imx
        do k=1,mask(i,j)
          b(i,j,k) = b(i,j,k) / dz(k)
        enddo
      enddo
    enddo
  endif
  if (nSqrt > 0) then
    do j=1,jmx
      do i=1,imx
        do k=1,mask(i,j)
          if (b(i,j,k) .gt. 0.0) then
            b(i,j,k) = sqrt(b(i,j,k))
            if (nSqrt .ge. 2) then
              b(i,j,k) = sqrt(b(i,j,k))
            endif
          else
            b(i,j,k) = 0.0;
          endif
        enddo
      enddo
    enddo
  endif
  if (zTpr) then
    zmx = zt(kmx)
    do j=1,jmx
      do i=1,imx
        kt = 0
        bmx = 0.0
        do k=1,mask(i,j)
          if (b(i,j,k) .gt. bmx) then
            bmx = b(i,j,k)
            kt = k
          endif
        enddo
        do k=kt,mask(i,j)
          b(i,j,k) = b(i,j,k)*(zmx-zt(k) + 0.25*(zt(k)-zt(kt))/(zmx-zt(kt)))
        enddo
      enddo
    enddo
  endif
!
  end subroutine fltrFile

! ===================================================================
!
  subroutine sfcMnMx
!
  integer :: klv
!
  do j=1,jmx
    do i=1,imx
      if (b(i,j,1) .lt. VsMn) then
        klv = -1
        do k=1,mask(i,j)
          if (b(i,j,k) .ge. VsMn) then
            klv = k
            exit
          endif
        enddo
        do k=1,klv-1
          b(i,j,k) = VsMn
        enddo
      else if (b(i,j,1) .gt. VsMx) then
        klv = -1
        do k=1,mask(i,j)
          if (b(i,j,k) .le. VsMx) then
            klv = k;
            exit
          endif
        enddo
        do k=1,klv-1
          b(i,j,k) = VsMx
        enddo
      endif
    enddo
  enddo
!
  end subroutine sfcMnMx
!
! ===================================================================
!
  subroutine writeFile
!
  integer :: ime, jme, kme

  ime = imx + 2
  jme = jmx + 2
  kme = kass
  allocate(ev(ime,jme,kme))
!
  do k=1,kme
    do j=1,jmx
      do i=1,imx
        ev(i+1,j+1,k) = b(i,j,k)
      enddo
    enddo
  enddo
  do k=1,kme
    do j=1,jme-1
      ev(1,j,k) = ev(ime-1,j,k)
      ev(ime,j,k) = ev(2,j,k)
    enddo
  enddo
  do k=1,kme
    do i=1,ime
      ev(i,1,k) = ev(i,2,k)
      ev(i,jme,k) = ev(i,jme-1,k)
    enddo
  enddo
!
  write(stamp,'(a,i2,a,i2,a,i4,a)') 'm/d/y=',month,'/',day,'/',year,', h:m:s= 0: 0: 0'
  open(11,file=fileOut,form='UNFORMATTED',access='SEQUENTIAL',status='REPLACE')
  write(11) stamp, period, ime, jme, kme, ev
  close(11)
!
  end subroutine writeFile
!
! ===================================================================

subroutine usage
!
  write(6,'(a)') 'Usage:'
  write(6,'(a)') 'mkEvNc -f fileIn -o fileOut ...'
  write(6,'(a)') '    -f fileIn    - input archive netCDF file'
  write(6,'(a)') '    -o fileOut   - output MOM3 file'
  write(6,'(a)') '    -k kass      - number of levels in assimilation'
  write(6,'(a)') ' Filter options:'
  write(6,'(a)') '    -SdZ         - scale by layer thickness'
  write(6,'(a)') '    -SqRt n      - take square root of Tz (n= 1 or 2)'
  write(6,'(a)') '    -Vsfc Mn Mx  - minimum/maximum values of variance at surface'
  write(6,'(a)') '    -zTpr        - taper portion of Tz below maximum'
  write(6,'(a)') '    -yTpr Ls Ln  - taper Tz south (north) of Ls (Ln)'
  write(6,'(a)') ' Smoothing options:'
  write(6,'(a)') '    -Av i j k    - # grid points in averaging window'
  write(6,'(a)') ' Normalization options:'
  write(6,'(a)') '    -LtN ltS ltN - latitude limits for normalizing'
  write(6,'(a)') '    -LgN lgW lgE - longitude limits for normalizing'
  write(6,'(a)') '    -gS real     - global scale factor'
  write(6,'(a)') '    -gC real     - global cieling (applied after gS)'
  write(6,*)
  write(6,'(a)') 'Defaults:'
  write(6,'(a,i2)') '    -k ', kass
  write(6,'(a,l4)') '    -SdZ', SdZ
  write(6,'(a,l4)') '    -zTpr', zTpr
  write(6,'(a,i2)') '    -SqRt ', nSqrt
  write(6,'(a,2f5.1)') '    -Vsfc ', VsMn, VsMx
  write(6,'(a,2i5)') '    -LtN ', ltS, ltN
  write(6,'(a,2i5)') '    -LgN ', lgW, lgE
  write(6,'(a,f5.1)') '    -gS ', gscl
  write(6,'(a,f5.1)') '    -gC ', gciel
  write(6,'(a,3i5)') '    -Av ', iav, jav, kav
  write(6,*)
!
  end subroutine usage
!
  end program mkEvNc
