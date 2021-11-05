  module gnutil_mod
!
! The module gnutil_mod (general utility) contains several functions
! and subroutines for public use.
!
  public expndLndMsk, expndOcnFld
  public intrp2d

  contains

! ----------------------------------------------------------------
!
      subroutine expndLndMsk (mask, im, jm, maxit, msk0)
!
!  This routine expands the land area represented in the array mask by 0's
!
      integer :: im, jm, maxit
      integer, dimension(im,jm) :: mask, msk0
      integer :: i, j, mxiter, iter, jsum, itot
!
      do j=1,jm
        do i=1,im
          msk0(i,j) = mask(i,j)
        enddo
      enddo
!
      mxiter = iabs(maxit)
      iter   = 0
      jsum = 1
      do while (iter .lt. mxiter .and. jsum .gt. 0)
        jsum = 0
        do j=2,jm-1
          do i=2,im-1
            if (mask(i,j) .ne. 0) then               ! check only ocean points
              itot = mask(i,j-1) + mask(i,j+1) &
                     & + mask(i-1,j) + mask(i+1,j)
              if (itot .lt. 4) then                  ! coastal ocean point
                msk0(i,j) = 0                        ! convert to land point
                jsum = jsum + 1
              endif
            endif
          enddo
        enddo
!
        do j = 2,jm-1
          do i = 2,im-1
            mask(i,j) = msk0(i,j)
          enddo
        enddo
        iter = iter + 1
      enddo
!
      end subroutine expndLndMsk
!
! ----------------------------------------------------------------
!
      subroutine expndOcnFld (fld, im, jm, x, y, mask, maxit, wrk, msk0, msk1)
!
!  This routine expands a data field in the array fld from ocean points
!  onto land points.  In the array mask, ocean points are 1's and land
!  points are 0's.
!
!  maxit.gt.0: no. of iterations to expand onto the land.
!  maxit.lt.0: don't change mask on return.
!  maxit.eq.0: fill all land points with extrapolated ocean values.
!
      integer :: im, jm
      integer, dimension(im,jm) :: mask, msk0, msk1
      real, dimension(im,jm) :: fld, wrk
      real, dimension(im) :: x
      real, dimension(jm) :: y
      integer :: i, j, mxiter, iter, jsum, itot
      real :: wtxm, wtxp, wtym, wtyp, dd
!
      do j=1,jm
        do i=1,im
          msk1(i,j) = mask(i,j)
          msk0(i,j) = mask(i,j)
          wrk(i,j) = fld(i,j)
        enddo
      enddo
!
      mxiter = iabs(maxit)
      iter = 0
      jsum = 1
      do while (jsum.gt.0 .and. (iter.lt.mxiter .or. mxiter.eq.0))
        jsum=0
        do j=2,jm-1
          do i=2,im-1
            if (msk1(i,j) .eq. 0) then                ! check only land points
              itot = msk1(i,j-1) + msk1(i,j+1) &
                       & + msk1(i-1,j) + msk1(i+1,j)
              if (itot .gt. 0) then                   ! coastal land point
                msk0(i,j) = 1                         ! convert to ocean point
                wtxm = (x(  i) - x(i-1))*msk1(i+1,j)
                wtxp = (x(i+1) - x(  i))*msk1(i-1,j)
                wtym = (y(  j) - y(j-1))*msk1(i,j+1)
                wtyp = (y(j+1) - y(  j))*msk1(i,j-1)
                dd   = wtxm+wtxp+wtym+wtyp
                wrk(i,j) = (wtxm*fld(i+1,j) + wtxp*fld(i-1,j) &
                            & + wtym*fld(i,j+1) + wtyp*fld(i,j-1))/dd
                jsum = jsum + 1
              endif
            endif
          enddo
        enddo
!
        do j=2,jm-1
          do i=2,im-1
            msk1(i,j) = msk0(i,j)
            fld (i,j) = wrk (i,j)
          enddo
        enddo
        iter = iter + 1
      enddo
!
      if (maxit .ge. 0) then
        do j = 1,jm
          do i = 1,im
            mask(i,j) = msk1(i,j)
          enddo
        enddo
      endif
!
      end subroutine expndOcnFld
!
! ----------------------------------------------------------------
!
      subroutine intrp2d (cn, cx, cy, ic, jc, gn, gx, gy, ig, jg, nO)
!
!  This routine is a simple 2D linear interpolation
!
      integer :: ic, jc, ig, jg, nO
      real, dimension(ic,jc) :: cn
      real, dimension(ic) :: cx
      real, dimension(jc) :: cy
      real, dimension(ig,jg) :: gn
      real, dimension(ig) :: gx
      real, dimension(jg) :: gy
      integer :: i, j, iw, ie, js, jn, igi
      real :: dxw, dxe, dx, dys, dyn, dy

      igi = ig - nO
      do j=1,jg
        js = 1
        do while (cy(js+1) .lt. gy(j))
          js = js + 1
        enddo
        jn = js + 1
        dyn = cy(jn) - gy(j)
        dys = gy(j) - cy(js)
        dy = cy(jn) - cy(js)
!
        do i = 1,igi
          if (cx(ic) .lt. gx(i)) then
            iw = ic - 1
            ie = 1
            dxe = cx(ie) + 360.0 - gx(i)
            dxw = gx(i) - cx(iw)
            dx = cx(ie) + 360.0 - cx(iw)
          else
            iw = 1
            do while (cx(iw+1) .lt. gx(i))
              iw = iw + 1
            enddo
            ie = iw + 1
            dxe = cx(ie) - gx(i)
            dxw = gx(i) - cx(iw)
            dx = cx(ie) - cx(iw)
          endif
          gn(i,j) = (dyn*(dxe*cn(iw,js) + dxw*cn(ie,js)) &
                     & + dys*(dxe*cn(iw,jn) + dxw*cn(ie,jn))) / (dx*dy)
        enddo
      enddo
!
      if (nO .gt. 0) then
        do j=1,jg
          do i = 1,nO
            iw = ig - nO + i
            gn(iw,j) = gn(i,j)
          enddo
        enddo
      endif
!
      end subroutine intrp2d
!
  end module gnutil_mod
