      subroutine Eg_lpsmthr
#ifdef assm_dta
!
!-----------------------------------------------------------------------
!   This subroutine multiplies g by an approximation to the first guess
!   error covariance matrix [e] to get the vector h. The approximation to
!   [e] is made by a series of multiplications by 1+laplacian.

!   At the start and end the decomposition is by latitude, using the 
!   decomposition computed for the model (the relevant information is
!   passed through taskrows.h). After some initial computation the array
!   s2v is transposed to the array s2h and the 2D horizontal laplacian 
!   smoother is applied under a vertical decomposition. Afterward, the
!   array s2h is transposed back to the array s2v.
!-----------------------------------------------------------------------
!
      use assim_mod
      use asm_x_mod
!
      include 'mpif.h'
# include "param.h"
# include "taskrows.h"
!
      integer :: tag = 1, req = 0, error
      integer :: stat(MPI_STATUS_SIZE)
      real :: tmpk(kass2)
!
      dpth = 1.0
!
      dbsq  = 0.5*aeval*aeval
      npid2 = npits/2
!
!-----------------------------------------------------------------------
!   load g into s2v
!-----------------------------------------------------------------------
!
      do j=jscomp,jecomp
        s2v(:,:,j) = g_cg(:,:,j)
      enddo
!
!-----------------------------------------------------------------------
!   multiply g by the square root of the local vertical background
!   error covariance
!-----------------------------------------------------------------------
!
      do j=jscomp,jecomp
        do i=2,imtm1
          call getcovsr
          do k=1,kass
            tmpk(k) = 0.
            do kk=1,kass
              tmpk(k) = tmpk(k) + s2v(i,kk,j) * covsr(k,kk)
            enddo
          enddo
# if defined cor_sal || defined asm_sal
          do k=kass+1,kass2
            tmpk(k) = 0.
            do kk=kass+1,kass2
              tmpk(k) = tmpk(k) + s2v(i,kk,j) * covsr(k,kk)
            enddo
          enddo
# endif
          do k=1,kass2
            s2v(i,k,j) = tmpk(k)
          enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!   TRANSPOSE vertical slabs (s2v) to horizontal slabs (s2h)
!-----------------------------------------------------------------------
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      do m=1,num_processors
        if (pn /= m) then
          call MPI_IRECV(s2buf,ns2b,MPI_REAL8,m-1,tag,
     &                                  MPI_COMM_WORLD,req,error)
          call MPI_WAIT(req, stat, error)
          do k=kscomp,kecomp
            do j=jscomp3(m),jecomp3(m)
              do i=1,imt
                s2h(i,j,k) = s2buf(i,j-jscomp3(m)+1,k-kscomp+1)
              enddo
            enddo
          enddo
        else
          do n=1,num_processors
            if (pn /= n) then
              do k=kscomp3(n),kecomp3(n)
                do j=jscomp,jecomp
                  do i=1,imt
                    s2buf(i,j-jscomp+1,k-kscomp3(n)+1) = s2v(i,k,j)
                  enddo
                enddo
              enddo
              call MPI_SEND(s2buf,ns2b,MPI_REAL8,n-1,tag,
     &                                     MPI_COMM_WORLD,error)
            else
              do k=kscomp,kecomp
                do j=jscomp,jecomp
                  do i=1,imt
                    s2h(i,j,k) = s2v(i,k,j)
                  enddo
                enddo
              enddo
            endif
          enddo
        endif
      enddo
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
!-----------------------------------------------------------------------
!   begin series of multiplications by 1+laplacians. note within each 
!   loop first multiplied by 1+laplacian then transpose of 1+laplacian 
!   to ensure symmetry of [e] matrix
!-----------------------------------------------------------------------
!
      do k=kscomp,kecomp
!
        temp(:,1)    = 0.0
        temp(:,jmt)  = 0.0
        s2h(:,1,k)   = 0.0
        s2h(:,jmt,k) = 0.0
!
        temp(1  ,:)  = 0.0
        temp(imt,:)  = 0.0
        s2h(1  ,:,k) = 0.0
        s2h(imt,:,k) = 0.0
!
!-----------------------------------------------------------------------
!   set up coefficients for central point of 5 point star.  attempt to
!    partially account for orography
!-----------------------------------------------------------------------
!
        do j=2,jmtm1
          do i=1,imt
            scra(i,j) = 1.0-wso(i,j)-wno(i,j)-2.0*wew(i,j)
          enddo
        enddo
        do j=2,jmtm1
          do i=1,imt
            con = wso(i,j)
            if (kma(i,j+1) .lt. k)      con = wno(i,j)
            col = con*con/((con+aeval)*con+dbsq)
            if (kma(i,j-1) .lt. k) scra(i,j) = scra(i,j)+con*col
            if (kma(i,j+1) .lt. k) scra(i,j) = scra(i,j)+con*col
          enddo
        enddo
        do j=2,jmtm1
          do i=2,imtm1
            con = wew(i,j)
            col = con*con*con/((con+aeval)*con+dbsq)
            if (kma(i-1,j) .lt. k) scra(i,j) = scra(i,j)+col
            if (kma(i+1,j) .lt. k) scra(i,j) = scra(i,j)+col
          enddo
        enddo
!
!-----------------------------------------------------------------------
!   set up dpth array
!-----------------------------------------------------------------------
!
        do j = 1,jmt
          do i = 1,imt
            if (kma(i,j).lt.k) dpth(i,j) = 0.0
          enddo
        enddo
!
!-----------------------------------------------------------------------
!   begin series of multiplications by 1+laplacians. note within each
!   loop first multiplied by 1+laplacian then transpose of 1+laplacian
!   to ensure symmetry of [e] matrix
!-----------------------------------------------------------------------
!
        do j=1,jmt
          do i=1,imt                                    
            s2h(i,j,k) = wgta(j)* s2h(i,j,k)
          enddo
        enddo
!
        do l=1,npid2
          do j=2,jmtm1                                                 
            do i=2,imtm1                                                
              temp(i,j) = ( scra(i,j)* s2h(i,j,k)
     &                    + wso(i,j)*  s2h(i,j-1,k)
     &                    + wno(i,j)*  s2h(i,j+1,k)
     &                    + wew(i,j)*( s2h(i+1,j,k)
     &                              + s2h(i-1,j,k) ) )
     &                    * dpth(i,j)
            enddo
          enddo
#ifdef cyclic
          do j=2,jmtm1
            temp(imt,j) = temp(2,j)
            temp(1,j)   = temp(imtm1,j)
          enddo
#endif
          do j=2,jmtm1                                  
            do i=2,imtm1                                         
              s2h(i,j,k) = ( scra(i,j)*temp(i,j)
     &                     + wno(i,j-1)*temp(i,j-1)
     &                     + wso(i,j+1)*temp(i,j+1)
     &                     + wew(i,j)*( temp(i+1,j)
     &                               + temp(i-1,j) ) )
     &                     * dpth(i,j)
            enddo
          enddo
#ifdef cyclic
          do j=2,jmtm1
            s2h(imt,j,k) = s2h(2,j,k)
            s2h(1,j,k)   = s2h(imtm1,j,k)
          enddo
#endif
        enddo
!
        do j=1,jmt
          do i=1,imt                                    
            s2h(i,j,k) = wgta(j)* s2h(i,j,k)
          enddo
        enddo
!
      enddo                        !        end of k-loop
!
!-----------------------------------------------------------------------
!   TRANSPOSE horizontal slabs (s2h) to vertical slabs (s2v)
!-----------------------------------------------------------------------
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      do m=1,num_processors
        if (pn /= m) then
          call MPI_IRECV(s2buf,ns2b,MPI_REAL8,m-1,tag,
     &                                  MPI_COMM_WORLD,req,error)
          call MPI_WAIT(req, stat, error)
          do k=kscomp3(m),kecomp3(m)
            do j=jscomp,jecomp
              do i=1,imt
                s2v(i,k,j) = s2buf(i,j-jscomp+1,k-kscomp3(m)+1)
              enddo
            enddo
          enddo
        else
          do n=1,num_processors
            if (pn /= n) then
              do k=kscomp,kecomp
                do j=jscomp3(n),jecomp3(n)
                  do i=1,imt
                    s2buf(i,j-jscomp3(n)+1,k-kscomp+1) = s2h(i,j,k)
                  enddo
                enddo
              enddo
              call MPI_SEND(s2buf,ns2b,MPI_REAL8,n-1,tag,
     &                                    MPI_COMM_WORLD,error)
            else
              do k=kscomp,kecomp
                do j=jscomp,jecomp
                  do i=1,imt
                    s2v(i,k,j) = s2h(i,j,k)
                  enddo
                enddo
              enddo
            endif
          enddo
        endif
      enddo
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
# if defined cor_sal || defined asm_sal
!
!-----------------------------------------------------------------------
!   Repeat process for salinity
!-----------------------------------------------------------------------
!   TRANSPOSE vertical slabs (s2v) to horizontal slabs (s2h)
!-----------------------------------------------------------------------
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      do m=1,num_processors
        if (pn /= m) then
          call MPI_IRECV(s2buf,ns2b,MPI_REAL8,m-1,tag,
     &                                  MPI_COMM_WORLD,req,error)
          call MPI_WAIT(req, stat, error)
          do k=kscomp,kecomp
            do j=jscomp3(m),jecomp3(m)
              do i=1,imt
                s2h(i,j,k) = s2buf(i,j-jscomp3(m)+1,k-kscomp+1)
              enddo
            enddo
          enddo
        else
          do n=1,num_processors
            if (pn /= n) then
              do k=kscomp3(n),kecomp3(n)
                do j=jscomp,jecomp
                  do i=1,imt
                    s2buf(i,j-jscomp+1,k-kscomp3(n)+1) = s2v(i,k+kass,j)
                  enddo
                enddo
              enddo
              call MPI_SEND(s2buf,ns2b,MPI_REAL8,n-1,tag,
     &                                     MPI_COMM_WORLD,error)
            else
              do k=kscomp,kecomp
                do j=jscomp,jecomp
                  do i=1,imt
                    s2h(i,j,k) = s2v(i,k+kass,j)
                  enddo
                enddo
              enddo
            endif
          enddo
        endif
      enddo
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
!-----------------------------------------------------------------------
!   begin series of multiplications by 1+laplacians. note within each
!   loop first multiplied by 1+laplacian then transpose of 1+laplacian
!   to ensure symmetry of [e] matrix
!-----------------------------------------------------------------------
!
      do k=kscomp,kecomp
!
        temp(:,1)    = 0.0
        temp(:,jmt)  = 0.0
        s2h(:,1,k)   = 0.0
        s2h(:,jmt,k) = 0.0
!
        temp(1  ,:)  = 0.0
        temp(imt,:)  = 0.0
        s2h(1  ,:,k) = 0.0
        s2h(imt,:,k) = 0.0
!
!-----------------------------------------------------------------------
!   set up coefficients for central point of 5 point star.  attempt to
!    partially account for orography
!-----------------------------------------------------------------------
!
        do j=2,jmtm1
          do i=1,imt
            scra(i,j) = 1.0-wso(i,j)-wno(i,j)-2.0*wew(i,j)
          enddo
        enddo
        do j=2,jmtm1
          do i=1,imt
            con = wso(i,j)
            if (kma(i,j+1) .lt. k)      con = wno(i,j)
            col = con*con/((con+aeval)*con+dbsq)
            if (kma(i,j-1) .lt. k) scra(i,j) = scra(i,j)+con*col
            if (kma(i,j+1) .lt. k) scra(i,j) = scra(i,j)+con*col
          enddo
        enddo
        do j=2,jmtm1
          do i=2,imtm1
            con = wew(i,j)
            col = con*con*con/((con+aeval)*con+dbsq)
            if (kma(i-1,j) .lt. k) scra(i,j) = scra(i,j)+col
            if (kma(i+1,j) .lt. k) scra(i,j) = scra(i,j)+col
          enddo
        enddo
!
!-----------------------------------------------------------------------
!   begin series of multiplications by 1+laplacians. note within each
!   loop first multiplied by 1+laplacian then transpose of 1+laplacian
!   to ensure symmetry of [e] matrix
!-----------------------------------------------------------------------
!
        do j=1,jmt
          do i=1,imt
            s2h(i,j,k) = wgta(j)* s2h(i,j,k)
          enddo
        enddo
!
#ifdef cyclic
        do j=2,jmtm1
          s2h(imt,j,k) = s2h(2,j,k)
          s2h(1,j,k)   = s2h(imtm1,j,k)
        enddo
#endif
!
        do l=1,npid2
          do j=2,jmtm1
            do i=2,imtm1
              temp(i,j) = ( scra(i,j)* s2h(i,j,k)
     &                    + wso(i,j)*  s2h(i,j-1,k)
     &                    + wno(i,j)*  s2h(i,j+1,k)
     &                    + wew(i,j)*( s2h(i+1,j,k)
     &                              + s2h(i-1,j,k) ) )
     &                    * dpth(i,j)
            enddo
          enddo
#ifdef cyclic
          do j=2,jmtm1
            temp(imt,j) = temp(2,j)
            temp(1,j)   = temp(imtm1,j)
          enddo
#endif
          do j=2,jmtm1
            do i=2,imtm1
              s2h(i,j,k) = ( scra(i,j)*temp(i,j)
     &                     + wno(i,j-1)*temp(i,j-1)
     &                     + wso(i,j+1)*temp(i,j+1)
     &                     + wew(i,j)*( temp(i+1,j)
     &                               + temp(i-1,j) ) )
     &                     * dpth(i,j)
            enddo
          enddo
#ifdef cyclic
          do j=2,jmtm1
            s2h(imt,j,k) = s2h(2,j,k)
            s2h(1,j,k)   = s2h(imtm1,j,k)
          enddo
#endif
        enddo
!
        do j=1,jmt
          do i=1,imt
            s2h(i,j,k) = wgta(j)* s2h(i,j,k)
          enddo
        enddo
!
      enddo                        !        end of k-loop
!
!-----------------------------------------------------------------------
!   TRANSPOSE horizontal slabs (s2h) to vertical slabs (s2v)
!-----------------------------------------------------------------------
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      do m=1,num_processors
        if (pn /= m) then
          call MPI_IRECV(s2buf,ns2b,MPI_REAL8,m-1,tag,
     &                                  MPI_COMM_WORLD,req,error)
          call MPI_WAIT(req, stat, error)
          do k=kscomp3(m),kecomp3(m)
            do j=jscomp,jecomp
              do i=1,imt
                s2v(i,k+kass,j) = s2buf(i,j-jscomp+1,k-kscomp3(m)+1)
              enddo
            enddo
          enddo
        else
          do n=1,num_processors
            if (pn /= n) then
              do k=kscomp,kecomp
                do j=jscomp3(n),jecomp3(n)
                  do i=1,imt
                    s2buf(i,j-jscomp3(n)+1,k-kscomp+1) = s2h(i,j,k)
                  enddo
                enddo
              enddo
              call MPI_SEND(s2buf,ns2b,MPI_REAL8,n-1,tag,
     &                                    MPI_COMM_WORLD,error)
            else
              do k=kscomp,kecomp
                do j=jscomp,jecomp
                  do i=1,imt
                    s2v(i,k+kass,j) = s2h(i,j,k)
                  enddo
                enddo
              enddo
            endif
          enddo
        endif
      enddo
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
#endif
!
!-----------------------------------------------------------------------
!   multiply by the square root of the local vertical background
!   error covariance
!-----------------------------------------------------------------------
!
      do j=jscomp,jecomp
        do i=2,imtm1
          call getcovsr
          do k=1,kass
            tmpk(k) = 0.
            do kk=1,kass
              tmpk(k) = tmpk(k) + s2v(i,kk,j) * covsr(kk,k)
            enddo
          enddo
# if defined cor_sal || defined asm_sal
          do k=kass+1,kass2
            tmpk(k) = 0.
            do kk=kass+1,kass2
              tmpk(k) = tmpk(k) + s2v(i,kk,j) * covsr(kk,k)
            enddo
          enddo
# endif
          do k=1,kass2
            s2v(i,k,j) = tmpk(k)
          enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!   save result back in h_cg
!-----------------------------------------------------------------------
!
      do j=jscomp,jecomp
        h_cg(:,:,j) = s2v(:,:,j)
# ifdef cyclic
        h_cg(imt,:,j) = h_cg(2,:,j)
        h_cg(1,:,j) = h_cg(imtm1,:,j)
# endif
      enddo
!
!---------------------------------------------------------------------
!    UPDATE HALO of h_cg
!    Only jecomp+1 of pn needs to be updated by jscomp of pn+1
!---------------------------------------------------------------------
!
      jecp1 = jecomp + 1
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      if (pn < num_processors) call MPI_IRECV(h_cg(:,1,jecp1),imtka,
     &                      MPI_REAL8,pn,tag,MPI_COMM_WORLD,req,error)
      if (pn > 1) call MPI_SEND(h_cg(:,1,jscomp),imtka,MPI_REAL8,
     &                                pn-2,tag,MPI_COMM_WORLD,error)
      if (pn < num_processors) call MPI_WAIT(req, stat, error)
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
!
!-----------------------------------------------------------------------
!   end of main part of Eg_lpsmthr
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!  vertical covariance
!-----------------------------------------------------------------------
!
      subroutine getcovsr
!
      dimension ev(kass)
!
      do l = 1,kass2
        do k = 1,kass2
          covsr(k,l) = 0.0
        enddo
      enddo
!
      if (kma(i,j).eq.0) return
!
      ka = min(kass, kma(i,j))
      do k = 1,ka
# ifdef fix_vv
        ev(k) = sqrt(vtmp(k,i,j))
# else
        ev(k) = sqrt(wprvv*vtmp(k,i,j,iprvvm) +
     &                         (1.0-wprvv)*vtmp(k,i,j,inxtvm))
# endif
      enddo
!
      if (kma(i,j).ge.kass) then
        do l = 1,kass
          do k = 1,kass
            covsr(k,l) = cvn(k,l) * ev(k) * ev(l)
          enddo
        enddo
      else
        do k = 1,ka
          covsr(k,k) = ev(k) * ev(k)
        enddo
      endif
!
# if defined cor_sal || defined asm_sal
!
! ------------------------------------------------------------------------
!  repeat for salinity
! ------------------------------------------------------------------------
      do k = 1,ka
#  ifdef fix_vv
        ev(k) = sqrt(vsal(k,i,j))
#  else
        ev(k) = sqrt(wprvv*vsal(k,i,j,iprvvm) +
     &                         (1.0-wprvv)*vsal(k,i,j,inxtvm))
#  endif
      enddo
!
      if (kma(i,j).ge.kass) then
        do l = kass+1,kass+ka
          do k = kass+1,kass+ka
            covsr(k,l) = 
     &               cvnsalt(k-kass,l-kass) * ev(k-kass) * ev(l-kass)
          enddo
        enddo
      else
        do k = kass+1,kass+ka
          covsr(k,k) = ev(k-kass) * ev(k-kass)
        enddo
      endif
# endif
        end subroutine getcovsr
#else
      real lps_dum
#endif
      end subroutine Eg_lpsmthr
