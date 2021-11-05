      subroutine combine_obs
#ifdef assm_dta
!
!-----------------------------------------------------------------------
!   This subroutine combines all observation types in preparation for
!   the initial comparison with the first guess.
!
!   Each processor handles only the data in the region assigned to it
!   by the decomposition.
!
!   The decomposition is by latitude, using the decomposition computed
!   for the model. The relevant information is passed through taskrows.h
!-----------------------------------------------------------------------
!
      use assim_mod
      use asm_x_mod
!
      include 'mpif.h'
#include "param.h"
#include "coord.h"
#include "grdvar.h"
#include "taskrows.h"
!
      nobs(:) = 0
      nsobs(:) = 0
      nsobs(jscomp) = 1
      nrsiz = (jecomp - jscomp + 1)*mobs
!
      do j=jscomp,jecomp
!
        if (j < jmtm1) then
!
!-----------------------------------------------------------------------
!   The call to get_all_obs merges all the data types for the row j
!     rtm   contains the observation times relative to the current model time
!     val   contains the observation values
!     aip   contains the i positions
!     ajp   contains the j positions
!     akp   contains the k positions
!     aerr  contains the estimated error values
!-----------------------------------------------------------------------
!
          call get_all_obs(j)
!
!-----------------------------------------------------------------------
!   modify expected error (aerr) based on time difference
!      T(z) has a +/- 14 day window  (akp >= 1)
!      S(z) has a +/- 14 day window  (akp >= 1+kass)
!      SST has a +/- 7 day window    (akp = 1)
!      SSH has a +/- 7 day window    (akp = 0)
!-----------------------------------------------------------------------
!
          do n=1,nobs(j)
            artm = abs(rtm(n))
!
            if (artm .gt. 7.0 .and. akp(n) .le. 1.5) aerr(n) = -99.0
            if (artm .gt. 14.0) aerr(n) = -99.0
!
            if (akp(n) .le. 1.5) then
              aerr(n) = aerr(n)*(1.0-artm*r7)
            else
              aerr(n) = aerr(n)*(1.0-artm*r14)
            endif
          enddo
          lmn = 0
          do n=1,nobs(j)
            if (aerr(n) .gt. 0.0) then
              lmn = lmn + 1
              rtm(lmn) = rtm(n)
              val(lmn) = val(n)
              aip(lmn) = aip(n)
              ajp(lmn) = ajp(n)
              akp(lmn) = akp(n)
              aerr(lmn) = aerr(n)
            endif
          enddo
          nobs(j) = lmn
          nsobs(j+1) = nsobs(j) + lmn
!
!-----------------------------------------------------------------------
!   Save the data
!-----------------------------------------------------------------------
!
          do n=1,nobs(j)
!
            if (nsobs(j)-1+n > nrsiz) then
              write(6,'(a,i3,a,2i8)')
     &          'pn=',pn,'   SIP SIZE EXCEEDED',nsobs(j)-1+n,nrsiz
              call flush(6)
            endif
!
            srtm(nsobs(j)-1+n) = rtm(n)
            sval(nsobs(j)-1+n) = val(n)
            sip(nsobs(j)-1+n)  = aip(n)
            sjp(nsobs(j)-1+n)  = ajp(n)
            skp(nsobs(j)-1+n)  = akp(n)
            serr(nsobs(j)-1+n) = aerr(n)
!
          enddo
!
        endif
!
      enddo
#else
      real :: cmbs_dum
#endif
      end subroutine combine_obs
