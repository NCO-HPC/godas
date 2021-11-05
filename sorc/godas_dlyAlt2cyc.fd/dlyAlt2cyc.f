!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  DLYALT2CYC
!   PRGMMR: David W. Behringer  ORG: NP23        DATE: 2006-08-08
!
! ABSTRACT:  Convert daily files of altimetry data to cycle files
!   as first step in preparing altimetry for use in Global Ocean
!   Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2008-08-08  David W. Behringer
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - LIST OF DAILY ALTIMETRY FILES IN ASCII
!     UNIT 21  - DAILY ALTIMETRY FILE IN ASCII. SEE REMARKS.
!
!   OUTPUT FILES:
!     UNIT 51  - ALTIMETRY CYCLE FILE IN ASCII. SEE REMARKS.
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   WORK FILES:  NONE
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - none
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - none
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage, errexit
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =  11 - ERROR OPENING UNIT 11
!     COND =  12 - ERROR READING UNIT 11
!     COND =  21 - ERROR OPENING UNIT 21
!     COND =  24 - ERROR BACKSPACING UNIT 21
!     COND =  51 - ERROR OPENING UNIT 51
!     COND =  52 - ERROR WRITING UNIT 51
!     COND =  53 - ERROR READING UNIT 51
!     COND =  54 - ERROR BACKSPACING UNIT 51
!
! REMARKS AND IMPORTANT LOCAL VARIABLES:
!     UNIT 21 is not assigned externally, file names are read from UNIT 11
!     UNIT 51 is not assigned externally, file names computed internally
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
!
      program dlyAlt2cyc
!
!  dlyAlt2cyc converts a list of daily ASCII files to cycle files. The daily
!  files may contain data from more than one cycle.
!
      implicit none
      integer :: dte, hr, mn, cycle, orbit
      integer :: odte, ohr, omn, ocycle, oorbit
      integer :: istatus
      integer :: nbad_order=0
      real :: lat, lon, ssh
      logical :: first, okay, exCF
      character(len=50) :: listFile, dayFile, cycleFile
      character(len=30) :: str
!
      call w3tagb('GODAS_DLYALT2CYC',2006,0164,0164,'NP23')
!
      first = .TRUE.
!
!  open an ASCII file listing the names of the daily ASCII files
!
      open(unit=11,form='FORMATTED', status='OLD', &
     &                             access='SEQUENTIAL',err=110)
      do while (.TRUE.)
!
!  read a daily file name from unit 11 and open the daily file. continue
!  until the list is exhausted, then jump to 200.
!
        read(11,'(a)',end=200,err=120) dayFile
        open(unit=21,file=trim(dayFile),form='FORMATTED', &
     &                      status='OLD', access='SEQUENTIAL',err=210)
        do while (.TRUE.)
!
!  read a record of daily file on unit 21.  continue until the file is
!  exhausted, then jump to 100
!
          read(21,'(i9,2i3,2i5,2f12.5,f8.3)',end=100,err=80)  &
     &                 dte, hr, mn, cycle, orbit, lat, lon, ssh
          if (first) then
!
!  if first record of first daily file, open first cycle file.  cycle files are
!  named with the convention: jb_xxx.txt where "xxx" is the cycle number.
!
            write(cycleFile,'(a4,i3.3,a4)',iostat=istatus) 'jb_c',cycle,'.txt'
            if(istatus.ne.0)then
              write(6,*) 'Error writing internal file cycleFile for cycle ',cycle
              cycle ! skip to next data record
            endif
            write(6,'(a,a)') 'Writing ', cycleFile
            inquire(file=trim(cycleFile),exist=exCF)
            if (exCF) then
              open(unit=51,file=trim(cycleFile),form='FORMATTED', &
     &          access='SEQUENTIAL',status='OLD',position='APPEND',err=510)
              backspace(51,err=540)
              read(51,'(i9,2i3,2i5,2f12.5,f8.3)',err=530) dte, hr, mn, &
     &                                cycle, orbit, lat, lon, ssh
              odte = dte
              ohr = hr
              omn = mn
              ocycle = cycle
              oorbit = orbit
              backspace(21,err=240)
            else
              odte = dte
              ohr = hr
              omn = mn
              ocycle = cycle
              oorbit = orbit
              open(unit=51,file=trim(cycleFile),form='FORMATTED', &
     &          access='SEQUENTIAL',status='NEW',err=510)
              write(51,'(i9,2i3,2i5,2f12.5,f8.3)',err=520) dte, hr, mn, &
     &                                cycle, orbit, lat, lon, ssh
            endif
            first = .FALSE.
          else if (cycle .eq. ocycle) then
!
!  if the cycle number is unchanged, continue writing to the same cycle
!  file.  only write a record if the orbit number, date and time remain
!  the same or increase.
!
            okay = .TRUE.
            if (dte .lt. odte) then
              okay = .FALSE.
            else if (dte .eq. odte .and. hr .lt. ohr) then
              okay = .FALSE.
            else if (dte .eq. odte .and. hr .eq. ohr .and. mn .lt. omn) then 
              okay = .FALSE.
            else if (orbit .lt. oorbit) then
              okay = .FALSE.
            endif
            if (okay) then
              write(51,'(i9,2i3,2i5,2f12.5,f8.3)',err=520) dte, hr, mn, &
     &                                  cycle, orbit, lat, lon, ssh
              odte = dte
              ohr = hr
              omn = mn
              oorbit = orbit
            else
              nbad_order=nbad_order+1
            endif
          else
!
!  if the cycle number changes, close unit 51 and open the next cycle file,
!  again on unit 51.
!
            close(51)
            write(cycleFile,'(a4,i3.3,a4)') 'jb_c',cycle,'.txt'
            write(6,'(a,a)') 'Writing ', cycleFile
            inquire(file=trim(cycleFile),exist=exCF)
            if (exCF) then
              open(unit=51,file=trim(cycleFile),form='FORMATTED', &
     &          access='SEQUENTIAL',status='OLD',position='APPEND',err=510)
              backspace(51,err=540)
              read(51,'(i9,2i3,2i5,2f12.5,f8.3)',err=530) dte, hr, mn, &
     &                                cycle, orbit, lat, lon, ssh
              odte = dte
              ohr = hr
              omn = mn
              ocycle = cycle
              oorbit = orbit
              backspace(21,err=240)
            else
              odte = dte
              ohr = hr
              omn = mn
              ocycle = cycle
              oorbit = orbit
              open(unit=51,file=trim(cycleFile),form='FORMATTED', &
     &          access='SEQUENTIAL',status='NEW',err=510)
              write(51,'(i9,2i3,2i5,2f12.5,f8.3)',err=520) dte, hr, mn, &
     &                                cycle, orbit, lat, lon, ssh
            endif
          endif
        enddo
        close(21)  ! this line & cycle shouldn't be reached, but included in case of future changes to above
        cycle
 80     continue  !  bad record found.  move on to the next file
        write(6,'(a,a)') 'Error reading a daily altimetry file on unit 21:',trim(dayFile)
        write(6,'(a)') 'Rest of file ignored.  Continue to process any remaining files.'
100     continue
        close(21)
      enddo
200   continue
      close(11)
      close(51)
      write (6,*) 
      write (6,*) nbad_order, 'OBS DISCARDED FOR BAD TIME SEQUENCE'
!
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(0)
!     call exit(0)
!
110   write(6,'(a)') 'Error opening list of daily altimetry files on unit 11'
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(11)
!     call exit(11)
!
120   write(6,'(a)') 'Error reading list of daily altimetry files on unit 11'
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(12)
!     call exit(12)
!
210   write(6,'(a)') 'Error opening a daily altimetry file on unit 21'
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(21)
!     call exit(21)
!
240   write(6,'(a)') 'Error backspacing daily altimetry file on unit 21'
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(24)
!     call exit(24)
!
510   write(6,'(a)') 'Error opening an altimetry cycle file on unit 51'
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(51)
!     call exit(51)
!
520   write(6,'(a)') 'Error writing an altimetry cycle file on unit 51'
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(52)
!     call exit(52)
!
530   write(6,'(a)') 'Error writing an altimetry cycle file on unit 51'
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(53)
!     call exit(53)
!
540   write(6,'(a)') 'Error backspacing an altimetry cycle file on unit 51'
      call w3tage('GODAS_DLYALT2CYC')
      call errexit(54)
!     call exit(54)
!
!
      end program dlyAlt2cyc
