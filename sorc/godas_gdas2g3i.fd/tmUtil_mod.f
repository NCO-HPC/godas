  module tmutil_mod
!
! The module tmutil_mod (time utility) contains several functions
! and subroutines for public use.
!
  public jlnDayR, jlnDayC, calDayR, calDayC, weekDay

  contains
!
! ---------------------------------------------------------------------
!
      function jlnDayR (year, month, day, year0)
!
      integer :: jlnDayR
      integer :: year, month, day, year0
      integer :: n, cn, cn0, yr, jdy
!
      cn0 = year0 / 100
      cn = year / 100
!
      if (cn0 .lt. cn) then
        yr = 100*cn0 + 99
        jdy = jlnDayC(yr,12,31) - jlnDayC(year0,1,1) + 1
!
        do n = cn0+1,cn-1
          yr = 100*n + 99
          jdy = jdy + jlnDayC(yr,12,31)
        enddo
!
        jdy = jdy + jlnDayC(year,month,day)
      else if (cn0 .eq. cn) then
        jdy = jlnDayC(year,month,day) - jlnDayC(year0,1,1) + 1
      else
        jdy = -1
      endif
!
      jlnDayR = jdy
!
      end function jlnDayR
!
! ---------------------------------------------------------------------
!
      function jlnDayC (year, month, day)
!
      integer :: jlnDayC
      integer :: year, month, day
      integer :: yr0, dyr, jdy
      integer, dimension(13), save :: smMn(13)
      data smMn /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/
!
      yr0 = (year / 100) * 100
      dyr = year - yr0
!
      jdy = 365 * dyr + (dyr-1)/4 + smMn(month) + day
!
      if (month .gt. 2 .and. mod(year,4) .eq. 0 &
                & .and. mod(year,100) .ne. 0) jdy = jdy + 1
!
      if (mod(yr0,400) .eq. 0 &
           & .and. (month .gt. 2 .or. year .ne. yr0)) jdy = jdy + 1
!
      jlnDayC = jdy
!
      end function jlnDayC
!
! ---------------------------------------------------------------------
!
      subroutine calDayR (year, month, day, year0, jday)
!
!     This routine returns the year, month, day, given the number of
!     days, jday, counting from the begining of year0
!
      integer :: year, month, day, year0, jday
      integer :: yr, mo, dy, yd, ytyp
      integer, save, dimension(13,2) :: smMn
      data smMn /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365, &
               & 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366/
!
      yr = year0
      dy = jday
!
      if (mod(yr,400) .eq. 0) then
        yd = 366
      else if (mod(yr,4) .eq. 0 .and. mod(yr,100) .ne. 0) then
        yd = 366
      else
        yd = 365
      end if
!
      do while (yd .lt. dy)
        dy = dy - yd
        yr = yr + 1
        if (mod(yr,400) .eq. 0) then
          yd = 366
        else if (mod(yr,4) .eq. 0 .and. mod(yr,100) .ne. 0) then
          yd = 366
        else
          yd = 365
        end if
      end do
!
      if (mod(yr,400) .eq. 0) then
        ytyp = 2
      else if (mod(yr,4) .eq. 0 .and. mod(yr,100) .ne. 0) then
        ytyp = 2
      else
        ytyp = 1
      endif
!
      mo = 1
      do while (dy .gt. smMn(mo+1,ytyp))
        mo = mo + 1
      enddo
!
      year = yr
      month = mo
      day = dy - smMn(month,ytyp)
!
      end subroutine calDayR
!
! ---------------------------------------------------------------------
!
      subroutine calDayC (year, month, day, jday)
!
!     This routine returns the month, day, given the number of days, 
!     jday, counting from the begining of the century in which year lies.
!
      integer :: year, month, day, jday
      integer :: yr0, yr, mo, dy, jjdy
!
      yr0 = (year / 100) * 100
      jdy = jday
!
      call calDayR(yr, mo, dy, yr0, jdy)
!
      year = yr
      month = mo
      day = dy
!
      end subroutine calDayC
!
! ---------------------------------------------------------------------
!
      function weekDay (year, month, day)
!
      integer :: weekDay
      integer :: year, month, day, rdow
      integer, parameter :: r19dow=2, r20dow=7
!
      if (year .ge. 2000) then
        rdow = r20dow
      else
        rdow = r19dow
      end if
!
      weekDay = mod((jlnDayC(year,month,day)+rdow-2),7) + 1
!
      end function weekDay
!
  end module tmutil_mod
