      subroutine comp_f
#ifdef assm_dta
!
!-----------------------------------------------------------------------
!   This routine calculates the f vector
!     f^n = e^n + trnsD invF D d^n
!-----------------------------------------------------------------------
!
      use assim_mod
!
      include 'mpif.h'
#include "param.h"
#include "taskrows.h"
!
      integer :: tag = 1, req = 0, error
      integer :: stat(MPI_STATUS_SIZE)
!
!-----------------------------------------------------------------------
!   Begin with the second half of the expression.
!   Set f to zero.
!-----------------------------------------------------------------------
!
      f_cg = 0.0
!
      do j=jscomp,jecomp
        if (j < jmtm1) then
          if (nobs(j) .ne. 0) then
            do n=1,nobs(j)
              aip(n)  = sip(nsobs(j)-1+n)
              ajp(n)  = sjp(nsobs(j)-1+n)
              akp(n)  = skp(nsobs(j)-1+n)
              aerr(n) = serr(nsobs(j)-1+n)
              il      = aip(n)
              ai      = il
              dx      = aip(n)-ai
              dy      = ajp(n)
              dxm1    = 1.0-dx
              dym1    = 1.0-dy
              ils(n)  = il
              a111(n) = dxm1*dym1
              a211(n) = dx*dym1
              a221(n) = dx*dy
              a121(n) = dxm1*dy
            enddo
!
!-----------------------------------------------------------------------
!   interpolate d to observation locations
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              il  = ils(n)
              ilp = il+1
#ifdef asm_ssh
              if (akp(n).gt.0.5) then
                kl = akp(n)
                ov(n) = a221(n) * d_cg(ilp,kl,j+1)
     &                  + a211(n) * d_cg(ilp,kl,j)
     &                  + a121(n) * d_cg(il,kl,j+1)
     &                  + a111(n) * d_cg(il,kl,j)
              else
                ov(n) = 0.0
                do k=1,kass
                  ov(n) = ov(n) + cdnz(k) *
     &                          ( a221(n) * d_cg(ilp,k,j+1)
     &                          + a211(n) * d_cg(ilp,k,j)
     &                          + a121(n) * d_cg(il,k,j+1)
     &                          + a111(n) * d_cg(il,k,j) )
# if defined cor_sal || defined asm_sal
     &                          + cdnzs(k) *
     &                          ( a221(n) * d_cg(ilp,k+kass,j+1)
     &                          + a211(n) * d_cg(ilp,k+kass,j)
     &                          + a121(n) * d_cg(il,k+kass,j+1)
     &                          + a111(n) * d_cg(il,k+kass,j) )
# endif
                enddo
              endif
#else
              kl = akp(n)
              ov(n) = a221(n) * d_cg(ilp,kl,j+1)
     &                + a211(n) * d_cg(ilp,kl,j)
     &                + a121(n) * d_cg(il,kl,j+1)
     &                + a111(n) * d_cg(il,kl,j)
#endif
            enddo
!
!-----------------------------------------------------------------------
!   multiply by inverse of F
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              ov(n) = ov(n) * aerr(n)
            enddo
!
!-----------------------------------------------------------------------
!   multiply by transpose of [d]
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              a221(n) = a221(n) * ov(n)
              a121(n) = a121(n) * ov(n)
              a211(n) = a211(n) * ov(n)
              a111(n) = a111(n) * ov(n)
            enddo
!
!-----------------------------------------------------------------------
!   accumulate result in f
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              il  = ils(n)
              ilp = il + 1
#ifdef cyclic
              if (il .eq. 1) il = imtm1
              if (ilp .eq. imt) ilp = 2
#endif
#ifdef asm_ssh
              if (akp(n) .gt. 0.5) then
                kl = akp(n)
                f_cg(ilp,kl,j+1) = f_cg(ilp,kl,j+1)+a221(n)
                f_cg(il,kl,j+1)  = f_cg(il,kl,j+1)+a121(n)
                f_cg(ilp,kl,j)   = f_cg(ilp,kl,j)+a211(n)
                f_cg(il,kl,j)    = f_cg(il,kl,j)+a111(n)
              else
                do k=1,kass
                  f_cg(ilp,k,j+1) = f_cg(ilp,k,j+1)+a221(n)*cdnz(k)
                  f_cg(il,k,j+1)  = f_cg(il,k,j+1)+a121(n)*cdnz(k)
                  f_cg(ilp,k,j)   = f_cg(ilp,k,j)+a211(n)*cdnz(k)
                  f_cg(il,k,j)    = f_cg(il,k,j)+a111(n)*cdnz(k)
# if defined cor_sal || defined asm_sal
                  kl = k+kass
                  f_cg(ilp,kl,j+1) = f_cg(ilp,kl,j+1)+a221(n)*cdnzs(k)
                  f_cg(il,kl,j+1)  = f_cg(il,kl,j+1)+a121(n)*cdnzs(k)
                  f_cg(ilp,kl,j)   = f_cg(ilp,kl,j)+a211(n)*cdnzs(k)
                  f_cg(il,kl,j)    = f_cg(il,kl,j)+a111(n)*cdnzs(k)
# endif
                enddo
              endif
#else
              kl  = akp(n)
              f_cg(ilp,kl,j+1) = f_cg(ilp,kl,j+1)+a221(n)
              f_cg(il,kl,j+1)  = f_cg(il,kl,j+1)+a121(n)
              f_cg(ilp,kl,j)   = f_cg(ilp,kl,j)+a211(n)
              f_cg(il,kl,j)    = f_cg(il,kl,j)+a111(n)
#endif
            enddo
!
          endif
        endif
!
      enddo
# ifdef cyclic
      f_cg(1,:,:)   = f_cg(imtm1,:,:)
      f_cg(imt,:,:) = f_cg(2,:,:)
# endif
!
!-----------------------------------------------------------------------
!   For processes pn=2,nproc transmit f_cg(:,:,jecomp+1) from pn-1 to a 
!   buffer on pn.  Add the buffer to f_cg(:,:,jscomp) on pn.
!-----------------------------------------------------------------------
!
      jecp1 = jecomp + 1
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
      if (pn > 1) then
        call MPI_IRECV(r_cg,imtka,MPI_REAL8,pn-2,tag,
     &                                    MPI_COMM_WORLD,req,error)
      endif
      if (pn < num_processors) then
        call MPI_SEND(f_cg(:,:,jecp1),imtka,MPI_REAL8,pn,tag,
     &                                    MPI_COMM_WORLD,error)
      endif
      if (pn > 1) then
        call MPI_WAIT(req, stat, error)
        f_cg(:,:,jscomp) = f_cg(:,:,jscomp) + r_cg(:,:)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
!-----------------------------------------------------------------------
!   add e
!-----------------------------------------------------------------------
!
      f_cg = f_cg + e_cg
!
#else
      real :: cf_dum
#endif
      end subroutine comp_f
