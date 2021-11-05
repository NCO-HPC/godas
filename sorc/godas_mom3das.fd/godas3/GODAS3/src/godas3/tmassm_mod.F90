  module timeassim_mod
!
!  The module timeassim_mod contains subroutine timeassim. 
!
!  A call of timeassim must have a unique index to refer to a particular
!  data set. 
!
!  Then in
!
!      call timeassim(time, index(i), ....
!
!  index(i) is the index of the data set and must be unique to that
!  data set.
!
!  Please note: In FORTRAN90 the type of the procedure parameters is
!  checked in more detail than in FORTRAN77. toprd and oprd must be
!  fields of type real, statements like
!      call timeassim(time, index, toprd(1,n) ..
!  are not accepted by some compilers since toprd(1,n) is a simple variable
!  of type real. Use
!      call timeassim(time, index, toprd(:,n) ..
!

  public position_within_obs, timeassim

  contains

      function position_within_obs (mom_time, obs_start_time, obs_end_time)
!
      use time_manager_mod
!
!=======================================================================
!     compute where the model time is in relation to starting and ending
!     time for the observations
!=======================================================================
!
      type(time_type) ::  mom_time, obs_start_time, obs_end_time
      type(time_type) ::  temp_time
      real :: position_within_obs
!
      if (mom_time .ge. obs_start_time) then
        call get_time (mom_time - obs_start_time, isec, iday)
      else
        call get_time (obs_start_time - mom_time, isec,iday)
        iday = -iday
        isec = -isec
      endif
!
      position_within_obs = iday + isec/86400.0
      return
      end function position_within_obs
!
!=======================================================================
!
      subroutine timeassim (tm, n, toprd, oprd, ndpd, iptod, ndpm)
!
!=======================================================================
!
!     time assimilation manager... constructs indices needed for
!     managing assimilation data defined at constant time intervals
!     (midpoints of weeks, typically) relative to the time of the 
!     current model time step.
!
!     inputs:
!
!     tm     = the time at which the data is desired (units of "toprd")
!
!     toprd  = the times at which the data periods in the dataset are
!              defined. times must be monotonically increasing and are
!              assumed to be at the centers of the data periods.
!
!     oprd  = constant length of the data periods
!              (eg: 7 days for weekly periods)
!
!     ndpd   = number of data periods in the dataset, i.e. on disk
!
!     n      = an index denoting which dataset is being assimilated
!              (each dataset should be referenced by a unique index
!               starting with 1 for the 1st, 2 for the 2nd, ...etc)
!
!     ndpm   = number of data periods in use at any time in memory
!
!     outputs:
!
!     iptod  = indices pointing to data periods on disk
!
!     author: d. w. behringer       e-mail=> david.behringer@noaa.gov
!=======================================================================
!
      dimension iptod(ndpm)
      parameter (maxsets=25, iflag=-99999)
      dimension ibold(maxsets), toprd(ndpd)
      data ibold /maxsets*iflag/
      save ibold
!
      if (n .gt. maxsets) then
        write (*,'(a,i10,a,i10)') 'Error: n=', n, ' maxsets=',maxsets
        stop '=>timeassim'
      endif
!
!       define the position of the dataset in time.
!      
      dstart = toprd(1) - 0.5*oprd
      dend   = toprd(ndpd) + 0.5*oprd
!
!       map the model time into the dataset. assume data is constant
!       before the beginning and after the end of the dataset
!
      if (tm .lt. dstart) then
        time = dstart
      elseif (tm .gt. dend) then
        time = dend 
      else
        time = tm
      endif
!
!     calculate record pointers for locating the dataset records
!     relative to the model time step.
!
!     find index of record whose center is closest to model time
!
      ndpm2 = ndpm/2
      ic = indp (time, toprd, ndpd)
      ib = ic - ndpm2
      if (ib .lt. 1) then
        ib = 1
      else if (ib .gt. ndpd-ndpm+1) then
        ib = ndpd-ndpm+1
      endif
      ibd = ib - ibold(n)
!
!     refresh pointers to memory buffers when reading disk data
!
      if (ibold(n) .ne. ib) then
        if (ibold(n) .eq. iflag) then
          do i = 1,ndpm
            iptod(i) = i+ib-1
          enddo
        else
!         do 
!           if (ibd .lt. ndpm) exit
!           ibd = ibd - ndpm
!         enddo
          do i = 1,ndpm
            iptod(i) = i+ib-1
          enddo
        endif
      endif
      ibold(n) = ib
!
      return
      end subroutine timeassim

      end module timeassim_mod

