      subroutine init_grad
#ifdef assm_dta
!
!-----------------------------------------------------------------------
!   This subroutine computes the initial estimate of the gradient of the
!   functional (g).
!-----------------------------------------------------------------------
!
      use assim_mod
      use asm_x_mod
!
      include 'mpif.h'
#include "param.h"
#include "tmngr.h"
#include "coord.h"
#include "grdvar.h"
#include "taskrows.h"
!
      integer :: tag = 1, req = 0, error
      integer :: stat(MPI_STATUS_SIZE)
!
!-----------------------------------------------------------------------
!   data types are encoded in akp 
!      T(z)        1 <  akp <= kass
!      S(z)   1+kass <= akp <= 2*kass
!      SST              akp  = 1
!      SSH              akp  = 0
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!   Set g to zero.
!-----------------------------------------------------------------------
!
      g_cg = 0.0
!
      do j=jscomp,jecomp
        if (j < jmtm1) then
          if (nobs(j) .ne. 0) then
#ifdef exRes
            nrtmp = 0
# ifdef asm_sal
            nrsal = 0
# endif
# ifdef asm_ssh                                                                 
            nreta = 0
# endif
#endif
            aj = j
            do n=1,nobs(j)
              rtm(n)  = srtm(nsobs(j)-1+n)
              val(n)  = sval(nsobs(j)-1+n)
              aip(n)  = sip(nsobs(j)-1+n)
              ajp(n)  = sjp(nsobs(j)-1+n)
              akp(n)  = skp(nsobs(j)-1+n)
              aerr(n) = serr(nsobs(j)-1+n)
              il      = aip(n)
              ai      = il
              dx      = aip(n)-ai
              dy      = ajp(n)-aj
              ajp(n)  = dy
              dxm1    = 1.0-dx
              dym1    = 1.0-dy
              ils(n)  = il
              a111(n) = dxm1*dym1
              a211(n) = dx*dym1
              a221(n) = dx*dy
              a121(n) = dxm1*dy
#ifdef exRes
              if (abs(rtm(n)) .le. rsdys) then
# ifdef asm_ssh
                if (akp(n).gt.0.5) then
                  kl    = akp(n)
                  if (kl .le. kass) then
                    nrtmp = nrtmp + 1
                  else
                    nrsal = nrsal + 1
                  endif
                else
                  nreta = nreta + 1
                endif
# else
#  ifdef asm_sal
                kl    = akp(n)
                if (kl .le. kass) then
                  nrtmp = nrtmp + 1
                else
                  nrsal = nrsal + 1
                endif
#  else
                nrtmp = nrtmp + 1
#  endif
# endif
              endif
#endif
            enddo
!
!-----------------------------------------------------------------------
!   interpolate model solution to observation pations
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              il    = ils(n)
              ilp   = il+1
#ifdef asm_ssh
              if (akp(n).gt.0.5) then
                kl    = akp(n)
                ov(n) = a221(n)*arex(ilp,kl,j+1)
     &                  + a211(n)*arex(ilp,kl,j)
     &                  + a121(n)*arex(il,kl,j+1)
     &                  + a111(n)*arex(il,kl,j)
              else
                ov(n) = a221(n)*dssh(ilp,j+1)
     &                  + a211(n)*dssh(ilp,j)
     &                  + a121(n)*dssh(il,j+1)
     &                  + a111(n)*dssh(il,j)
              endif
#else
              kl    = akp(n)
              ov(n) = a221(n)*arex(ilp,kl,j+1)
     &                + a211(n)*arex(ilp,kl,j)
     &                + a121(n)*arex(il,kl,j+1)
     &                + a111(n)*arex(il,kl,j)
#endif
            enddo
!
#ifdef exRes
!
!-----------------------------------------------------------------------
!   save model values, observations and residuals to files
!-----------------------------------------------------------------------
!
            if (nassim .eq. 1) then
              if (nrtmp .gt. 0)
     &          write(ioextr,'(1x,a3,i7,1x,a32,1x,a2,i7)')
     &                            'ts=',itt,stamp,'n=',nrtmp
# ifdef asm_sal
              if (nrsal .gt. 0)
     &          write(ioexsr,'(1x,a3,i7,1x,a32,1x,a2,i7)')
     &                            'ts=',itt,stamp,'n=',nrsal
# endif
# ifdef asm_ssh                                                                 
              if (nreta .gt. 0)
     &          write(ioexer,'(1x,a3,i7,1x,a32,1x,a2,i7)')
     &                            'ts=',itt,stamp,'n=',nreta                    
# endif
              do n=1,nobs(j)
                if (abs(rtm(n)) .le. rsdys) then
                  kl      = akp(n)
                  il      = aip(n)
                  ai      = il
                  dx      = aip(n)-ai
                  dy      = ajp(n)
                  xo = xt(il) + dx*( xt(il+1) - xt(il) )
                  yo = yt(j) + dy*( yt(j+1) - yt(j) )
# ifdef asm_ssh
                  if (kl .eq. 0) then
                    zo = 0.0                                          
! DEBUG
!     write(ioexer,'(8i8)') 
!    & 1, il, il+1, imt, jstask, j, j+1, jetask
!     write(ioexer,'(5f8.2,4f8.3)') 
!    & rtm(n),xt(il),xt(il+1),yt(j),yt(j+1),
!    & dssh(il,j),dssh(il,j+1),dssh(il+1,j),dssh(il+1,j+1)
! DEBUG
                    write(ioexer,'(4f8.2,3f8.3)') rtm(n),
     &                     xo, yo, zo, ov(n), val(n), ov(n) - val(n)
                  else if (kl .le. kass) then
                    zo = zt(kl)*1.0e-2
                    write(ioextr,'(4f8.2,3f8.3)') rtm(n),
     &                     xo, yo, zo, ov(n), val(n), ov(n) - val(n)
                  else
                    zo = zt(kl-kass)*1.0e-2
                    write(ioexsr,'(4f8.2,3f8.3)') rtm(n),
     &                     xo, yo, zo, 35.0+1.0e3*ov(n),
     &                     35.0+1.0e3*val(n), 1.0e3*(ov(n) - val(n))
                  endif
# else
#  ifdef asm_sal
                  if (kl .le. kass) then
                    zo = zt(kl)*1.0e-2
                    write(ioextr,'(4f8.2,3f8.3)') rtm(n),
     &                     xo, yo, zo, ov(n), val(n), ov(n) - val(n)
                  else
                    zo = zt(kl-kass)*1.0e-2
                    write(ioexsr,'(4f8.2,3f8.3)') rtm(n),
     &                     xo, yo, zo, 35.0+1.0e3*ov(n), 
     &                     35.0+1.0e3*val(n), 1.0e3*(ov(n) - val(n))
                  endif
#  else
                  zo = zt(kl)*1.0e-2
                  write(ioextr,'(4f8.2,3f8.3)') rtm(n),
     &                     xo, yo, zo, ov(n), val(n), ov(n) - val(n)
#  endif
# endif
                endif
              enddo
            endif
#endif
!
!-----------------------------------------------------------------------
!   calculate difference between model solution and observation
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              ov(n) = ov(n) - val(n)
            enddo
            do n=1,nobs(j)
              val(n) = abs(ov(n))
            enddo
!
!-----------------------------------------------------------------------
!   Adjust estimated error based on difference between model and obs.
!-----------------------------------------------------------------------
!   for temperature:
!     if (dT > 10) 1./error = 0.
!     if (dtlim < dT < 10)  1/error = 1./(error*(1.+abs(dT)-dtlim)**2)
!-----------------------------------------------------------------------
# ifdef Nasm_sal
!   for salinity:
!     if (dS > 5) 1./error = 0.
!     if (dslim < dS < 5)  1/error = 1./(error*(1.+abs(dS)-dslim)**2)
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              kl = akp(n)
              if (kl .le. kass) then
                if (val(n) .gt. 10.0) then
                  aerr(n) = -99.0
                else if (val(n) .gt. dtlim) then
                  aerr(n) = aerr(n)/(1.+val(n)-dtlim)**2
                endif
              else
                if (val(n) .gt. 5.0) then
                  aerr(n) = -99.0
                else if (val(n) .gt. dslim) then
                  aerr(n) = aerr(n)/(1.+val(n)-dslim)**2
                endif
              endif
            enddo
# else
!
            do n=1,nobs(j)
              if (val(n) .gt. 10.0) then
                aerr(n) = -99.0
              else if (val(n) .gt. dtlim) then
                aerr(n) = aerr(n)/(1.+val(n)-dtlim)**2
              endif
            enddo
# endif
!
!-----------------------------------------------------------------------
!   Strip out observations that are too far from the model.
!-----------------------------------------------------------------------
!
            lmn = 0
            do n=1,nobs(j)
              if (aerr(n) .gt. 0.0) then
                lmn = lmn + 1
                ov(lmn) = ov(n)
                a111(lmn) = a111(n)
                a211(lmn) = a211(n)
                a221(lmn) = a221(n)
                a121(lmn) = a121(n)
                ils(lmn) = ils(n)
                aip(lmn)  = aip(n)
                ajp(lmn)  = ajp(n)
                akp(lmn)  = akp(n)
                aerr(lmn) = aerr(n)
              endif
            enddo
            nobs(j) = lmn
!
!-----------------------------------------------------------------------
!   Multiply difference by inverse of estimated error and project back
!   on grid by using transpose of D
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              ov(n) = ov(n)*aerr(n)
            enddo
            do n=1,nobs(j)
              a221(n) = a221(n)*ov(n)
              a211(n) = a211(n)*ov(n)
              a111(n) = a111(n)*ov(n)
              a121(n) = a121(n)*ov(n)
            enddo
!
!-----------------------------------------------------------------------
!   Accumulate the result in g
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              il  = ils(n)
              ilp = il+1
#ifdef cyclic
              if (il .eq. 1) il = imtm1
              if (ilp .eq. imt) ilp = 2
#endif
#ifdef asm_ssh
              if (akp(n).gt.0.5) then
                kl = akp(n)
                g_cg(ilp,kl,j+1) = g_cg(ilp,kl,j+1)+a221(n)
                g_cg(il,kl,j+1)  = g_cg(il,kl,j+1)+a121(n)
                g_cg(ilp,kl,j)   = g_cg(ilp,kl,j)+a211(n)
                g_cg(il,kl,j)    = g_cg(il,kl,j)+a111(n)
              else
                do k=1,kass
                  g_cg(ilp,k,j+1)  = g_cg(ilp,k,j+1)+a221(n)*cdnz(k)
                  g_cg(il,k,j+1)   = g_cg(il,k,j+1)+a121(n)*cdnz(k)
                  g_cg(ilp,k,j)    = g_cg(ilp,k,j)+a211(n)*cdnz(k)
                  g_cg(il,k,j)     = g_cg(il,k,j)+a111(n)*cdnz(k)
# if defined cor_sal || defined asm_sal
                  kl = k+kass
                  g_cg(ilp,kl,j+1) = g_cg(ilp,kl,j+1)+a221(n)*cdnzs(k)
                  g_cg(il,kl,j+1)  = g_cg(il,kl,j+1)+a121(n)*cdnzs(k)
                  g_cg(ilp,kl,j)   = g_cg(ilp,kl,j)+a211(n)*cdnzs(k)
                  g_cg(il,kl,j)    = g_cg(il,kl,j)+a111(n)*cdnzs(k)
# endif
                enddo
              endif
#else
              kl = akp(n)
              g_cg(ilp,kl,j+1) = g_cg(ilp,kl,j+1)+a221(n)
              g_cg(il,kl,j+1)  = g_cg(il,kl,j+1)+a121(n)
              g_cg(ilp,kl,j)   = g_cg(ilp,kl,j)+a211(n)
              g_cg(il,kl,j)    = g_cg(il,kl,j)+a111(n)
#endif
            enddo
!
!-----------------------------------------------------------------------
!   Save data positions and errors again to account for data that may 
!   have been stripped out.
!-----------------------------------------------------------------------
!
            do n=1,nobs(j)
              sip(nsobs(j)-1+n)  = aip(n)
              sjp(nsobs(j)-1+n)  = ajp(n)
              skp(nsobs(j)-1+n)  = akp(n)
              serr(nsobs(j)-1+n) = aerr(n)
            enddo
!
          endif
        endif
!
      enddo
# ifdef cyclic
      g_cg(1,:,:)   = g_cg(imtm1,:,:)
      g_cg(imt,:,:) = g_cg(2,:,:)
# endif
!
!-----------------------------------------------------------------------
!   For processes pn=2,nproc transmit g_cg(:,:,jecomp+1) from pn-1 to a 
!   buffer on pn.  Add the buffer to g_cg(:,:,jscomp) on pn.
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
        call MPI_SEND(g_cg(:,:,jecp1),imtka,MPI_REAL8,pn,tag,
     &                                    MPI_COMM_WORLD,error)
      endif
      if (pn > 1) then
        call MPI_WAIT(req, stat, error)
        g_cg(:,:,jscomp) = g_cg(:,:,jscomp) + r_cg(:,:)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
#else
      real :: fg_dum
#endif
      end subroutine init_grad
