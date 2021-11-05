  program mkSEv
!
  implicit none
!
! readFile, writeFile, usage
!

  integer :: imx, jmx, kmx
  integer :: iav, jav, kav
  real(kind=4), allocatable, dimension(:) :: xt, yt, zt, dz
  real(kind=4), allocatable, dimension(:,:,:) :: a, b
  integer, allocatable, dimension(:,:) :: mask
  character(len=32) :: stamp
  character(len=121) :: str
  character(len=91) :: fileIn, fileOut
  real(kind=4) :: Tz, fctr, period
  integer :: i, j, k, km, n, narg
  logical :: prdc

  real(kind=4), parameter :: spv = 999.999

! set some defaults
  Tz = 0.4
  fctr = 3.162e-5
  prdc = .true.
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
    else if (str(1:2) .eq. '-t') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(f8.2)') Tz
    else if (str(1:2) .eq. '-c') then
      n = n + 1
      call getarg(n,str)
      if (n .gt. narg .or. str(1:1) .eq. '-') then
        call usage
        call exit(99)
      endif
      read(str,'(f12.0)') fctr
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

  call readFile

  do j=1,jmx
    do i=1,imx
      if (mask(i,j) .gt. 0) then
        km = 0
        do k=1,mask(i,j)
          if (b(i,j,k) .gt. Tz) km = k
        enddo
        do k=1,km
          b(i,j,k) = 1.0
        enddo
        km = km + 1
        do k=km,mask(i,j)
          b(i,j,k) = b(i,j,k) / Tz
        enddo
      endif
    enddo
  enddo

  do j=1,jmx
    do i=1,imx
      do k=1,mask(i,j)
         b(i,j,k) = b(i,j,k) * fctr
      enddo
    enddo
  enddo

  call writeFile

  contains
!
! ===================================================================
!
  subroutine readFile
!
  open(41,file=fileIn,form='UNFORMATTED',access='SEQUENTIAL',status='OLD')
  read(41) stamp, period, imx, jmx, kmx
  allocate(b(imx,jmx,kmx))
  rewind(41)
! close(41)
! open(41,file=fileIn,form='UNFORMATTED',access='SEQUENTIAL',status='OLD')

  read(41) stamp, period, imx, jmx, kmx, b
  close(41)

  allocate(mask(imx,jmx))
  do j=1,jmx
    do i=1,imx
      mask(i,j) = 0
      do k=1,kmx
        if (b(i,j,k) .ne. spv) then
          mask(i,j) = k
        endif
      enddo
    enddo
  enddo
!
  end subroutine readFile
!
! ===================================================================
!
  subroutine writeFile
!
  open(11,file=fileOut,form='UNFORMATTED',access='SEQUENTIAL',status='REPLACE')
  write(11) stamp, period, imx, jmx, kmx, b
  close(11)
!
  end subroutine writeFile

! ===================================================================

subroutine usage
!
  write(6,'(a)') 'Usage:'
  write(6,'(a)') 'mkEvNc -f fileIn -o fileOut ...'
  write(6,'(a)') '    -f fileIn    - input temperature variance file'
  write(6,'(a)') '    -o fileOut   - output salinity variance file'
  write(6,'(a)') '    -t Tz        -  lower limit on input Tz field'
  write(6,'(a)') '    -c fctr      -  multiply field by constant factor'
  write(6,*)
  write(6,'(a)') 'Defaults:'
  write(6,'(a,f6.2)') '    -t ', Tz
  write(6,'(a,1pe10.3)') '    -c ', fctr
!
  end subroutine usage
!
  end program mkSEv
