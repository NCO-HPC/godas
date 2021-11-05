      subroutine extract_guess (joff, js, je)
#ifdef assm_dta
!
!=======================================================================
!     extract 1st guess field from inside the mw
!
!     input:
!
!      joff   = offset between row j in the MW and latitude jrow on disk
!      js     = starting row for calculations
!      je     = ending row for calculations
!=======================================================================
!
      use arrays_2d
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
              arex(i,k,j+joff) = t(i,k,j,1,tau)
#if defined cor_sal || defined asm_sal
              arex(i,k+kass,j+joff) = t(i,k,j,2,tau)
#endif
            enddo
          enddo
        enddo
!
#ifdef asm_ssh
        do j=js,je
          do i=1,imt
            etax(i,j+joff) = etat(i,j+joff,1)
          enddo
        enddo
#endif
!
      endif
#else
      real :: eg_dum
#endif
      end subroutine extract_guess
