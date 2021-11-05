      subroutine add_correction (joff, js, je)
#ifdef assm_dta
!
!=======================================================================
!     add the correction inside the mw
!
!     input:
!
!      joff   = offset between row j in the MW and latitude jrow on disk
!      js     = starting row for calculations
!      je     = ending row for calculations
!=======================================================================
!
      use memory_window
      use assim_mod
      use asm_x_mod
!
#include "param.h"
#include "mw.h"
!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------
!
      if (js .le. je) then
!
        do j=js,je
          do k=1,kass
            do i=1,imt
              t(i,k,j,1,tau) = t(i,k,j,1,tau) + ares(i,k,j+joff)
#if defined cor_sal || defined asm_sal
              t(i,k,j,2,tau) = t(i,k,j,2,tau) + ares(i,k+kass,j+joff)
#endif
            enddo
          enddo
        enddo
!
      endif
#else
      real :: ac_dum
#endif
      end subroutine add_correction
