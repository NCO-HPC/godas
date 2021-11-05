!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  MRGPRF
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2003-07-25
!
! ABSTRACT:  Perform quality control on subsurface temperature data
!   for use in Global Ocean Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2003-07-25  David W. Behringer
! 2004-01-30  David W. Behringer - bug fix in routine merge_files
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - TEMPERATURE PROFILE DATA IN IEEE
!     UNIT 12  - TEMPERATURE PROFILE DATA IN IEEE
!
!   OUTPUT FILES:
!     UNIT 51  - MERGED TEMPERATURE PROFILE DATA IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  (INCLUDING SCRATCH FILES)
!     None
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - merge_files
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - merge_files
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  21 - ERROR OPENING UNIT 12
!     COND =  22 - ERROR READING UNIT 12
!     COND =  51 - ERROR OPENING UNIT 51
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
      program mrgPrf
!
!  mrgPrf merges profile files. It assumes the files are sorted by time.
!
      character csign1*8,sid1*2,dtyp1*2,qkey1*1
      real, allocatable, dimension(:) :: t1, z1
      character csign2*8,sid2*2,dtyp2*2,qkey2*1
      real, allocatable, dimension(:) :: t2, z2
!
      call w3tagb('GODAS_MRGPRF',2003,0164,0164,'NP23')
!
! open profile file
!
      open (11, form='unformatted', status='old', &
                            & access='sequential', err=110)
!
      npmx = 0
      do while (.true.)
        read (11, end=100, err=120) iyear,idate,csign1,sid1, &
                                       &  dtyp1,qkey1,yp,xp,np
        if (np .gt. npmx) npmx = np
      end do
  100 continue
      rewind (11)
!
! open 2nd profile file
!
      open (12, form='unformatted', status='old', &
                            & access='sequential', err=210)
!
      do while (.true.)
        read (12, end=200, err=220) iyear,idate,csign2,sid2, &
                                       &  dtyp2,qkey2,yp,xp,np
        if (np .gt. npmx) npmx = np
      end do
  200 continue
      rewind (12)
!
      allocate(t1(npmx))
      allocate(z1(npmx))
      allocate(t2(npmx))
      allocate(z2(npmx))
!
! open profile file
!
      open (51, form='unformatted', access='sequential', err=510)
      call merge_files(11,12,51)
!
      close (11)
      close (12)
      close (51)
!
      call w3tage('GODAS_MRGPRF')
      call errexit(0)
!     call exit(0)
!
  110 write(6,'(a)') 'Error opening profile file on unit 11'
      call w3tage('GODAS_MRGPRF')
      call errexit(11)
!     call exit(11)
!
  120 write(6,'(a)') 'Error reading profile file on unit 11'
      call w3tage('GODAS_MRGPRF')
      call errexit(12)
!     call exit(12)
!
  210 write(6,'(a)') 'Error opening profile file on unit 12'
      call w3tage('GODAS_MRGPRF')
      call errexit(21)
!     call exit(21)
!
  220 write(6,'(a)') 'Error reading profile file on unit 12'
      call w3tage('GODAS_MRGPRF')
      call errexit(22)
!     call exit(22)
!
  510 write(6,'(a)') 'Error opening profile file on unit 51'
      call w3tage('GODAS_MRGPRF')
      call errexit(51)
!     call exit(51)
!
      contains
!
! -------------------------------------------------------------- 
!
      subroutine merge_files(nu1,nu2,nu3)
!
      integer nu1, nu2, nu3
      integer(kind=8) :: iyd1, iyd2
      integer(kind=8) :: scl=100000000
      logical flg1, flg2
!
      flg1 = .true.
      flg2 = .true.
!
      read (nu1) iyear1,idate1,csign1,sid1,dtyp1,qkey1,y1,x1, &
                           &  np1,(z1(k),t1(k),k=1,np1)
      iyd1 = iyear1*scl + idate1
      read (nu2) iyear2,idate2,csign2,sid2,dtyp2,qkey2,y2,x2, &
                           &  np2,(z2(k),t2(k),k=1,np2)
      iyd2 = iyear2*scl + idate2
!
      do while (flg1 .or. flg2)
        if (flg1 .and. flg2) then
          if (iyd1 .le. iyd2) then
            write (nu3) iyear1,idate1,csign1,sid1,dtyp1,qkey1,y1,x1, &
                           &  np1,(z1(k),t1(k),k=1,np1)
            read (nu1, iostat=ios) iyear1,idate1,csign1,sid1,dtyp1, &
                           &  qkey1,y1,x1,np1,(z1(k),t1(k),k=1,np1)
            iyd1 = iyear1*scl + idate1
            if (ios .ne. 0) flg1 = .false.
          else
            write (nu3) iyear2,idate2,csign2,sid2,dtyp2,qkey2,y2,x2, &
                           &  np2,(z2(k),t2(k),k=1,np2)
            read (nu2, iostat=ios) iyear2,idate2,csign2,sid2,dtyp2, &
                           &  qkey2,y2,x2,np2,(z2(k),t2(k),k=1,np2)
            iyd2 = iyear2*scl + idate2
            if (ios .ne. 0) flg2 = .false.
          end if
        else if (flg1) then
          write (nu3) iyear1,idate1,csign1,sid1,dtyp1,qkey1,y1,x1, &
                           &  np1,(z1(k),t1(k),k=1,np1)
          read (nu1, iostat=ios) iyear1,idate1,csign1,sid1,dtyp1, &
                           &  qkey1,y1,x1,np1,(z1(k),t1(k),k=1,np1)
          iyd1 = iyear1*scl + idate1
          if (ios .ne. 0) flg1 = .false.
        else if (flg2) then
          write (nu3) iyear2,idate2,csign2,sid2,dtyp2,qkey2,y2,x2, &
                           &  np2,(z2(k),t2(k),k=1,np2)
          read (nu2, iostat=ios) iyear2,idate2,csign2,sid2,dtyp2, &
                           &  qkey2,y2,x2,np2,(z2(k),t2(k),k=1,np2)
          iyd2 = iyear2*scl + idate2
          if (ios .ne. 0) flg2 = .false.
        end if
      end do
!
      rewind (nu1)
      rewind (nu2)
      rewind (nu3)
!
      end subroutine merge_files
!
      end program mrgPrf
