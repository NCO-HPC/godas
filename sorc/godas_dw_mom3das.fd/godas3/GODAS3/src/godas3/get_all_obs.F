      subroutine get_all_obs(jrow)
#ifdef assm_dta
!
!----------------------------------------------------------------------
!   get_all_obs merges the specified number of weeks for each data type
!   into one dataset.  It is called from combine_obs.
!----------------------------------------------------------------------
!
      use assim_mod
!
# include "param.h"
# include "taskrows.h"
!
      j = jrow
      jobs = 1
      aobs(:,:) = 0.0
!
#ifdef asm_sst
      do m = 1,nwksst
        nwk = ipsstd(m)
        joff = itobs(jrow,2,nwk)
        jnum = itobs(jrow,1,nwk)
!
        if (jobs + jnum - 1 > mobs) then
          write(6,'(a)') 'Error reading SST obs in get_aobs'
          write(6,'(a)') 'Maximum number of obs per row exceeded'
          jbc = jobs + jnum - 1
          write(6,'(a,i4,i6,a,i6)') 'Row=',j,jbc,' >',mobs
          call flush(6)
          call abort()
        endif
        if (jnum > 0) then
          do n=1,jnum
            read(iosstao,REC=joff+n)(aobs(i,jobs),i=1,6)
            jobs = jobs + 1
          enddo
        endif
      enddo
#endif
#ifdef asm_tmp
      do m = 1,nwktmp
        nwk = iptmpd(m)
        joff = jtobs(jrow,2,nwk)
        jnum = jtobs(jrow,1,nwk)
!
        if (jobs + jnum - 1 > mobs) then
          write(6,'(a)') 'Error reading temperature obs in get_aobs'
          write(6,'(a)') 'Maximum number of obs per row exceeded'
          jbc = jobs + jnum - 1
          write(6,'(a,i4,i6,a,i6)') 'Row=',jrow,jbc,' >',mobs
          call flush(6)
          call abort()
        endif
        if (jnum > 0) then
          do n=1,jnum
            read(iotmpao,REC=joff+n)(aobs(i,jobs),i=1,6)
            jobs = jobs + 1
          enddo
        endif
      enddo
#endif
#ifdef asm_sal
      do m = 1,nwksal
        nwk = ipsald(m)
        joff = ktobs(jrow,2,nwk)
        jnum = ktobs(jrow,1,nwk)
!
        if (jobs + jnum - 1 > mobs) then
          write(6,'(a)') 'Error reading salinity obs in get_aobs'
          write(6,'(a)') 'Maximum number of obs per row exceeded'
          jbc = jobs + jnum - 1
          write(6,'(a,i4,i6,a,i6)') 'Row=',j,jbc,' >',mobs
          call flush(6)
          call abort()
        endif
        if (jnum > 0) then
          do n=1,jnum
            read(iosalao,REC=joff+n)(aobs(i,jobs),i=1,6)
            jobs = jobs + 1
          enddo
        endif
      enddo
#endif
#ifdef asm_ssh
      do m = 1,nwkssh
        nwk = ipsshd(m)
        joff = ltobs(jrow,2,nwk)
        jnum = ltobs(jrow,1,nwk)
!
        if (jobs + jnum - 1 > mobs) then
          write(6,'(a)') 'Error reading SSH obs in get_aobs'
          write(6,'(a)') 'Maximum number of obs per row exceeded'
          jbc = jobs + jnum - 1
          write(6,'(a,i4,i6,a,i6)') 'Row=',j,jbc,' >',mobs
          call flush(6)
          call abort()
        endif
        if (jnum > 0) then
          do n=1,jnum
            read(iosshao,REC=joff+n)(aobs(i,jobs),i=1,6)
            jobs = jobs + 1
          enddo
        endif
      enddo
#endif
      jobs = jobs - 1
      do n = 1,jobs
        rtm(n) = aobs(1,n) - dayasm
        val(n) = aobs(2,n)
        aip(n) = aobs(3,n)
        ajp(n) = aobs(4,n)
        akp(n) = aobs(5,n)
        aerr(n) = aobs(6,n)
      enddo
      nobs(j) = jobs
#else
      real :: gadum
#endif
      end subroutine get_all_obs
