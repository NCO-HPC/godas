      subroutine assim
#ifdef assm_dta
!
!   A VERTICAL SLAB ROUTINE
!
      use time_manager_mod
#ifndef fix_vv
      use timeinterp_mod
#endif
      use assim_mod
      use asm_x_mod
      use timeassim_mod
!
!-----------------------------------------------------------------------
!   Subroutine assim is the main driving routine for data assimilation
!
!   This version uses MPI message passing only.
!
!   All comments concerning data assimilation refer to paper
!   Derber and Rosati (1989) A global oceanic data assimilation system.
!   in J. Physical Oceanography.
!-----------------------------------------------------------------------
!
      include 'mpif.h'
#include "param.h"
#include "tmngr.h"
#include "switch.h"
#include "taskrows.h"
!
      integer :: tag = 1, req = 0, error
      integer :: stat(MPI_STATUS_SIZE)
      character*9 chstr
!
!-----------------------------------------------------------------------
!   beta   = beta in paper, factor controlling amount of previous search
!            direction to be included in current direction
!   gh_old = <g,h at previous iteration (< ,  represents the inner product)
!   gh_new = <g,h at current iteration
!   df_val = <d,f at current iteration
!   alpha  = alpha in paper, stepsize
!-----------------------------------------------------------------------
!
      real :: beta, gh_new, gh_old, df_val, alpha
!
      if (eorun) return
!
#ifndef fix_vv
!-----------------------------------------------------------------------
!     Update background error covriance data
!-----------------------------------------------------------------------
!
      call update_cov
!
#endif
      if (rassim) nassim = 0
!
      if (nassim .ge. 1 .and. nassim .le. nassdur) then
!
!-----------------------------------------------------------------------
!  Update data pointers
!-----------------------------------------------------------------------
!
        call update_obs_ptrs
!
        if (nassim.eq.1) write(stdout,999) itt,stamp
!
!-----------------------------------------------------------------------
!  Gather data falling within specified time windows
!-----------------------------------------------------------------------
!
        call combine_obs
!
# ifdef asm_ssh
!-----------------------------------------------------------------------
!  Subtract mean climatology from the current model sea surface height
!-----------------------------------------------------------------------
!
        dssh = etax - sshc
! DEBUG
!     write(stdout,'(a)') 'DEBUG DSSH ETAX SSHC'
!     do j=jstask,jetask
!     do i=1,imt
!     if (abs(dssh(i,j)) .ge. 100.0) then
!     write(stdout,'(2i4,a,1pe12.3,a,e12.3,a,e12.3)') 
!    & i,j,'  ETAX',etax(i,j),'  SSHC',sshc(i,j),'  DSSH',dssh(i,j)
!     write(stdout,'(a,1p5e12.3)') 'ETAX',(etax(i,j),i=178,182)
!     write(stdout,'(a,1p5e12.3)') 'SSHC',(sshc(i,j),i=178,182)
!     write(stdout,'(a,1p5e12.3)') 'DSSH',(dssh(i,j),i=178,182)
!     endif
!     enddo
!     enddo
!     write(stdout,'(a)') 'DEBUG END'
! DEBUG
# endif
!
!-----------------------------------------------------------------------
!  Find the first interation of the gradient of the functional (g^1)
!  by setting the intial guess for the correction field to zero (T^1 = 0),
!  comparing the model with the observations, weighting their difference
!  with the inverse of the observation error covariance (F) and projecting
!  this onto the model grid.
!       T^1 = 0
!       g^1 = -trnsD invF To
!-----------------------------------------------------------------------
!
        t_cg = 0.0
!
        call init_grad
!
!-----------------------------------------------------------------------
!  Do the first multiplication of the gradient by the background
!  error covariance matrix (E).
!       h^1 = E g^1
!  In this version a laplace smoother is used.
!-----------------------------------------------------------------------
!
        call Eg_lpsmthr
!
!-----------------------------------------------------------------------
!  Set the initial search directions to zero.
!       d^0 = 0
!       e^0 = 0
!-----------------------------------------------------------------------
!
        d_cg = 0.0
        e_cg = 0.0
!
!-----------------------------------------------------------------------
!  Set the initial value of beta to zero
!-----------------------------------------------------------------------
!
        beta = 0.0
!
!-----------------------------------------------------------------------
!  Begin the iteration loop
!-----------------------------------------------------------------------
!
        miter = maxits
        do iter=1,miter
!
!-----------------------------------------------------------------------
!  Update the search directions
!-----------------------------------------------------------------------
!
          d_cg = beta * d_cg - h_cg
          e_cg = beta * e_cg - g_cg
!
!-----------------------------------------------------------------------
!  Compute f
!      f^n = e^n + trnsD invF D d^n
!-----------------------------------------------------------------------
!
          call comp_f
!
!-----------------------------------------------------------------------
!  Compute the inner products <g,h  and <d,f and update alpha
!  (only over the computational part of the processor domain)
!-----------------------------------------------------------------------
!
          gh_new = 0.0
          df_val = 0.0
          do j=jscomp,jecomp
            do k=1,kass2
              do i=2,imtm1
                gh_new = gh_new + g_cg(i,k,j) * h_cg(i,k,j)
                df_val = df_val + d_cg(i,k,j) * f_cg(i,k,j)
              enddo
            enddo
          enddo
!
!---------------------------------------------------------------------
!    gh_new and df_val are REDUCED and DISPERSED to all processes
!---------------------------------------------------------------------
!
          call MPI_BARRIER(MPI_COMM_WORLD, error)
          call MPI_ALLREDUCE(gh_new,ssum,1,MPI_REAL8,MPI_SUM,
     &                                  MPI_COMM_WORLD,error)
          gh_new = ssum
!
          call MPI_BARRIER(MPI_COMM_WORLD, error)
          call MPI_ALLREDUCE(df_val,ssum,1,MPI_REAL8,MPI_SUM,
     &                                  MPI_COMM_WORLD,error)
          df_val = ssum
!
          alpha = gh_new / df_val
!
!-----------------------------------------------------------------------
!  Update the field correction (T) and the gradient (g)
!      T^(n+1) = T^n + alpha d^n
!      g^(n+1) = g^n + alpha f^n
!-----------------------------------------------------------------------
!
          t_cg = t_cg + alpha * d_cg
!
          if (iter .lt. miter) then
            g_cg = g_cg + alpha * f_cg
!
!-----------------------------------------------------------------------
!  Update h by multiplying the new gradient ( g^(n+1) ) by the
!  background error covariance E.
!       h^(n+1) = E g^(n+1)
!  In this version a laplace smoother is used.
!-----------------------------------------------------------------------
!
            call Eg_lpsmthr
!
!-----------------------------------------------------------------------
!  Compute a new inner product <g,h and update beta
!  (only over the computational part of the processor domain)
!-----------------------------------------------------------------------
!
            gh_old = gh_new
            gh_new = 0.0
            do j=jscomp,jecomp
              do k=1,kass2
                do i=2,imtm1
                  gh_new = gh_new + g_cg(i,k,j) * h_cg(i,k,j)
                enddo
              enddo
            enddo
!
!---------------------------------------------------------------------
!    gh_new is REDUCED and DISPERSED to all processes
!---------------------------------------------------------------------
!
            call MPI_BARRIER(MPI_COMM_WORLD, error)
            call MPI_ALLREDUCE(gh_new,ssum,1,MPI_REAL8,MPI_SUM,
     &                                  MPI_COMM_WORLD,error)
            gh_new = ssum
!
            beta = gh_new / gh_old
!
          endif
!
        enddo
!-----------------------------------------------------------------------
!  End of iteration loop
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Put field correction back in ares
!-----------------------------------------------------------------------
!
        ares = t_cg
!
!---------------------------------------------------------------------
!    UPDATE HALO of ares
!---------------------------------------------------------------------
!
        itk2 = 2*imtka
        jrecv = jecomp + 1
        jsend = jscomp
!
        call MPI_BARRIER(MPI_COMM_WORLD, error)
!
        if (pn < num_processors) call MPI_IRECV(ares(:,1,jrecv),itk2,
     &                      MPI_REAL8,pn,tag,MPI_COMM_WORLD,req,error)
        if (pn > 1) call MPI_SEND(ares(:,1,jsend),itk2,MPI_REAL8,
     &                                pn-2,tag,MPI_COMM_WORLD,error)
        if (pn < num_processors) call MPI_WAIT(req, stat, error)
!
        call MPI_BARRIER(MPI_COMM_WORLD, error)
!
        jrecv = jstask
        jsend = jecomp - 1
!
        if (pn > 1) call MPI_IRECV(ares(:,1,jrecv),itk2, MPI_REAL8,
     &                             pn-2,tag,MPI_COMM_WORLD,req,error)
        if (pn < num_processors) call MPI_SEND(ares(:,1,jsend),itk2,
     &                        MPI_REAL8, pn,tag,MPI_COMM_WORLD,error)
        if (pn > 1) call MPI_WAIT(req, stat, error)
!
        call MPI_BARRIER(MPI_COMM_WORLD, error)
!
      endif
      nassim = nassim + 1
!
  999 format(1x,'ts=',i7,1x,a32,', data assimilation')
!
!-----------------------------------------------------------------------
!   end of main part of assim
!-----------------------------------------------------------------------
!
      contains
!
#ifndef fix_vv
!-----------------------------------------------------------------------
!  updating the background error covariance data
!-----------------------------------------------------------------------
!
      subroutine update_cov
!
!-----------------------------------------------------------------------
!     determine the disk pointers, time weight interpolation factor,
!     and whether or not it is time to bring in new variance data from
!     disk based on the time (days)
!-----------------------------------------------------------------------
      dayvv = position_within_data (model_time, vv_start_time
     &,                              vv_end_time, vvprd)
      call timeinterp (dayvv, indxvv, vvrec, vprec, nvvrec, vvprd
     &,         mthdvv, inxtvd, iprvvd, wprvv, rdvv, inxtvm, iprvvm)

      if (rdvv) then
!
!    read in the next data record from disk
!
        do j = jscomp,jecomp
          nrec = j + (inxtvd-1)*jmt
          read (iotvv,rec=nrec) bufik
          do i = 1,imt
            do k = 1,kass
              vtmp(k,i,j,inxtvm) = bufik(i,k)
            enddo
          enddo
# if defined cor_sal || defined asm_sal
          read (iosvv,rec=nrec) bufik
          do i = 1,imt
            do k = 1,kass
              vsal(k,i,j,inxtvm) = bufik(i,k)
            enddo
          enddo
# endif
        enddo
!
        write (stdout,9000) inxtvd, inxtvm, stamp
9000    format (/1x,'==> Reading temperature variance from rec ',i3
     &, ' to index ',i3,' at MOM time ',a,/)
# if defined cor_sal || defined asm_sal
        write (stdout,9001) inxtvd, inxtvm, stamp
9001    format (/1x,'==> Reading salinity variance from rec ',i3
     &, ' to index ',i3,' at MOM time ',a,/)
# endif
      endif
!
      end subroutine update_cov
#endif
!
!-----------------------------------------------------------------------
!  updating the pointers to observations
!-----------------------------------------------------------------------
!
      subroutine update_obs_ptrs
!
!-----------------------------------------------------------------------
!  All the sets of observations are organised in weekly periods
!  Therefore, all will be due for updating at the same time
!  No explicit updating is necessary, because at each assimilation
!  observations are read from the word addressable disk files. The
!  weeks that are read are given by the arrays ip***d that are kept
!  up to date by the calls to timeassim.
!-----------------------------------------------------------------------
!
      dayasm = position_within_obs (model_time, obs_start_time
     &,                              obs_end_time)
!
#ifdef asm_sst
      call timeassim (dayasm, ndxsst, toprd, oprd, nwkobs
     &,                ipsstd, nwksst)
#endif
#ifdef asm_tmp
      call timeassim (dayasm, ndxtmp, toprd, oprd, nwkobs
     &,                iptmpd, nwktmp)
#endif
#ifdef asm_sal
      call timeassim (dayasm, ndxsal, toprd, oprd, nwkobs
     &,                ipsald, nwksal)
#endif
#ifdef asm_ssh
      call timeassim (dayasm, ndxssh, toprd, oprd, nwkobs
     &,                ipsshd, nwkssh)
#endif
!
      end subroutine update_obs_ptrs
!
#else
      real :: a_dum
#endif
      end subroutine assim
