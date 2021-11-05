      subroutine checks (errorc, vmixset, htr_mixset, hvel_mixset)
      use mom3_mpp_mod
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
      logical errorc, vmixset, htr_mixset, hvel_mixset
#include "param.h"
#include "accel.h"
#include "coord.h"
#include "derived_options.h"
#include "grdvar.h"
#include "hmixc.h"
#include "index.h"
#include "iounit.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "taskrows.h"
#include "vmixc.h"
c
      common / checks_common / icg, jcg, cgmax, gridmin, sml
c
c-----------------------------------------------------------------------
c     do consistency checks before allowing MOM to continue
c         
c     author:   R. C. Pacanowski     e-mail=> rcp@gfdl.gov
c    (mpp items)   v. balaji         e-mail=> vb@gfdl.gov
c-----------------------------------------------------------------------
c
      write (stdout,'(/,20x,a,/,8x,a,a,/)')
     & 'G E N E R A L    C O N S I S T E N C Y    C H E C K S'
     &,'(A "rudimentary intelligence" analysis of all specifications'
     &,' and their interactions)'
c
c
c     explicit free surface option checks
c
      write(stdout,9000)
     & '==> Warning: Option explicit_free_surface_nonlinear is enabled.'
     &,'             dzt(k=1)+eta must be positive in order to        '
     &,'             preclude vanishingly small top cell thicknesses. '
     &,'             If dzt(k=1)+eta < 0, the model will not conserve '
     &,'             first moments.                                   '
      write(stdout,*) 
     & '             A nominal minimum thickness of ', thickmin,'(cm) '
      write(stdout,9000)     
     & '             is used to keep the model from blowing up in case'
     &,'             of extremely small top cell thicknesses.         '
     &,'             Yet if the model is producing such thin top cells,'
     &,'             you should investigate the reason for such.       '
      if(dtuv .ne. dtts) then
      write(stdout,9000)
     & '==> Warning: Option explicit_free_surface_nonlinear is enabled'
     &,'             yet dtts is not equal to dtuv.  In order to      '
     &,'             have a conservative model, set dtts=dtuv.        '
      endif

	  write (stdout,9000)
     & '==> Warning: because options explicit_eb or explicit_efb have '
     &,'             not been enabled, the explicit free surface '
     &,'             solution will be time averaged over two baroclinic'
     &,'             timesteps to keep the barotropic solution defined'
     &,'             at the same time as the baroclinic solution. This'             
     &,'             will allow stretching of tracer timesteps relative'
     &,'             to baroclinic timesteps.'
c
c     fresh water flux option checks
c
c
c
c
c       options which have not been parallelized  
c
      if (imt .lt. 3) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  parameter "imt" less than 3 is not allowed'        
        errorc = .true.
      endif
c
      if (jmt .lt. 4) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  parameter "jmt" less than 4 is not allowed'        
        errorc = .true.
      endif
      if (num_processors .gt. 1 .and. num_processors .le. jmt-2) then
        write (stdout,'(/,1x,a,i4,a)')
     & '==> Warning:  "num_processors" =',num_processors
     &,' implies multi-tasking (parallel processing)'
      elseif (num_processors .eq. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning:  "num_processors" = 1 implies unitasking          '
      elseif (num_processors .gt. jmt-2 .or. num_processors .lt. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  "num_processors" cannot be > jmt-2 or < 1          '
        write (stdout,*)' num_processors = ',num_processors
        write (stdout,*) 
     &  'num_processors is set via namelist "procs" in run script '
	errorc = .true.
      endif
c
c
c     more options not checked out in parallel
c
c
c     some misc option checks
c
      if(smooth .le. 0.0) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: The smooth parameter for robert_time_filter is      '
     &,'           currently set <= 0.0. Model will eventually go      '
     &,'           unstable with this setting. A value of smooth=0.01  '
     &,'           should be sufficient for numerical stability.       '
        errorc = .true.
      endif
c
c
c       Check for obsolete options
c
c
c Check initial condition options
c
      n = 0
      n = n + 1
c mjh
      if (n .ne. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  Either "levitus_ic" or "idealized_ic" or "read_ic"'
     &,'            or "equatorial_thermocline" must be enabled.      '
     &,'             choose only one for initial conditions.          '
        errorc = .true.
      endif
c
c     Check boundary condition options
c
c
      n = 0
      n = n + 1
      if (n .ne. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  Either "simple_sbc" or "time_mean_sbc_data" or    '
     &,'            "time_varying_sbc_data" or "coupled" must be      '
     &,'             enabled. choose one surface boundary condition   '
        errorc = .true.
      endif
c
c     Check lateral mixing options
c 
      if (htr_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tracer_horz_mix_const"  cannot be enabled        '
     &,'    because another horizontal mixing scheme has been enabled.'
        errorc     = .true.
      else
        htr_mixset = .true.
      endif
      if (.not.htr_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No horizontal tracer diffusivity scheme enabled.  '
        errorc = .true.
      endif
      if (hvel_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "velocity_horz_mix_var"  cannot be enabled        '
     &,'    because another horizontal mixing scheme has been enabled.'
        errorc     = .true.
      else
        hvel_mixset = .true.
      endif
      if (.not.hvel_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No horizontal velocity viscosity scheme enabled.  '
        errorc = .true.
      endif
c
c     Check vertical mixing options
c
      if (.not. vmixset) then
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a/1x,a/1x,a))')
     & '==> Error: "bryan_lewis_vertical" tracer diffusion requires   '
     &,'           "constvmix", "ppvmix", "tcvmix", or "kppvmix".     '
     &,'           If "constvmix" is enabled, "bryan_lewis_vertical"  '
     &,'           coefficients replace kappa_h. If "ppvmix","tcvmix",'
     &,'           or "kppvmix" is enabled, "bryan_lewis_vertical"    '
     &,'           coefficients are added to the predicted coeffs.    '
        errorc = .true.
      endif
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No vertical mixing scheme has been enabled  '
        errorc = .true.
      endif
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: there is no vertical mixing scheme enabled         '
        errorc = .true.
      endif

	write (stdout,9000)
     & '==> Warning: option partial_cell uses the UNESCO equation of '
     &,'             state for hydrostatic pressure gradients.       '
c
c       Options for memory management of diagnostics 
c
          write(stdout,9000)
     & '==> Note: "diagnostics_himem" is enabled.                    '
     &,'          Diagnostics of 3d fields will be buffered over the '
     &,'          full model.   This approach requires much memory,  ' 
     &,'          but is considerably faster than diagnostics_lomem. '
c
c       Check for incompatible options
c        
c
c
        write (stdout,9000)
     & '==> Warning:  "parallel_1d" is enabled.                        '
        write (stdout,*) '              Note:', num_processors
     &, ' processors were requested.'
c
c     look for potential load imbalance 
c
      ncr = jecomp3(1)-jscomp3(1)+1
      i   = 0
      do n=2,num_processors
        if (ncr .ne. (jecomp3(n)-jscomp3(n)+1)) i = 1
      enddo
      if (i .eq. 1) then
        write (stdout,*)
     & '==> Warning: There is a potential load imbalance if the number '
     &,'             of computed rows is not the same on all processors'
        do n=1,num_processors
          ncr = jecomp3(n)-jscomp3(n)+1
          write (stdout,*) 'Processor #',n,' computes ',ncr,' rows'
        enddo
      endif
      if (jmw .ne. 5 .and. jmw .ne. window_size) then
        write (stdout,9000)
     & '==> Error: jmw is set wrong for a fourth_order_window'
     &,'           only the following sizes are allowed:'
        write (stdout,*)'          jmw=5 or jmw=',window_size
        write (stdout,*)'          you have set jmw=',jmw
	write (stdout,*)'          reset jmw to ',window_size
	write (stdout,*)'          in size.h'
        errorc = .true.  
      else
        write (stdout,9000)
     & '==> Note: a fourth_order_window is being used.'
      endif
c
      if (nkflds .lt. 2) then
        write (stdout,9000)
     & '==> Error: "nkflds" must be at least 2                        '
        write (stdout,*)'          nkflds is set = ',nkflds
        errorc = .true.  
      endif
c
        write (stdout,9000)
     & '==> Warning: The topography is being imported. Are you sure   '
     &,'             that the "kmt" field matches the vertical grid?  ' 
     &,'             i.e. Are the grid points within levels given by  '
     &,'             k = kmt(i,jrow) really defined at depth "zt(k)"? ' 
c
      if (damp_factor(1) .eq. c0 .or. damp_factor(2) .eq. c0) then
        write (stdout,9000)
     & '==> Error: the damping time scale "dampts" must be set > 0    '
     &,'           when enabling option "restorst"                    '
        errorc = .true.
      endif
      if (damp_factor(1) .ne. damp_factor(2)) then
        write (stdout,9000)
     & '==> Warning: are you aware that the damping time scale "dampts"'
     &,'           is set differently for T and S in option "restorst"?'
     &,'           (they can be set the same in the namelist)          '
      endif
c
c
c
        write (stdout,9000)
     & '==> Note: since redi_diffusion was specified for tracers,     '
     &,'    kappa_h and Ah will be used as background mixing coeffs   '
c
c
      if(eta_mix.le.100) then
        write (stdout,9000)
     &  '==> Warning:  Laplacian filtering of eta is on.              '  
     & ,'              eta_mix = ',eta_mix,' may not be a good choice!'
     & ,'              Check the output for the free surface part!    '
     & ,'              Check vertical velocities for checkerboards!   '
     & ,'              Redefine eta_mix in namelist.mixing.           '
      endif
c
c
c
c
c
        if ((ah+ahisop) .gt. 1.e11) then
          write (stdout,9000)
     & '==> Error: "ahisop"+"ah" is too large for the                 '
     &,'             "redi_diffusion" mixing option                   '
          errorc = .true.
        endif
c
        if (dtsf .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the external mode time step "dtsf"     '
          errorc = .true.
          dtsf = epsln
	endif
c
        if (dtuv .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the internal mode time step "dtuv"     '
          errorc = .true.
          dtuv = epsln
	endif
c
        if (dtts .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the density time step "dtts"           '
          errorc = .true.
          dtts = epsln
	endif
c
        critv = 1.e-6
        if (mod(rundays,dtts*secday) .gt. critv) then
          t1 = nint(rundays/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: there must be an integral number of density time '
     &,'    steps within "rundays" (the integration time).            '
          write (stdout,*) '               (changed "rundays" from '
     &,   rundays,' days to ', t1,' days to insure this condition)'
          rundays = t1
        endif
c
c-----------------------------------------------------------------------
c       time step checks
c-----------------------------------------------------------------------
c
      if (acor .ne. 0) then
        write (stdout,9000)
     & '==> Error: "acor" must=0 when option damp_inertial_oscillation'
     &,'           is not enabled.  '
          errorc = .true.
      else
c
c       check for marginally resolved inertial oscillation
c
        fmax = epsln
        do jrow=jscomp,jecomp
	  do i=2,imtm1
	    if (kmu(i,jrow) .ne. 0) then
              fmax = max(fmax,abs(c2*omega*sine(jrow)))
	    endif
	  enddo
        enddo
        call mpp_max(fmax)
        call mpp_sync()
        dtuvmax = (1.0/6.0)*(c2*pi)/fmax
        write (stdout,'(/a,/,11x,a,f8.0,a)')
     &' ==> Note: 6 timesteps/(min inertial period) implies a maximum'
     &,'dtuv of ',dtuvmax,' sec.'
        if (dtuv .gt. dtuvmax) then
          write (stdout,9000)
     & '==> Error: the inertial oscillation is not resolved. reduce    '
     &,'           "dtuv" or use option "damp_inertial_oscillation"    '
     &,'            or remove stop condition to continue.'
          errorc = .true.
        endif
      endif
c
c       CFL time step limit for internal gravity waves
c
        cgint = 300.0
	dtcg  = 1.e10
	icg   = 1
	jcg   = jmt/2
        dxumin = 1.0e20
	do jrow=jscomp,jecomp
	  do i=2,imt-1
	    if (kmu(i,jrow) .ne. 0) then
	      dxumin = min(dxumin,dyu(jrow),csu(jrow)*dxu(i))
	      dti = 0.5*dxumin/cgint
	      if (dti .lt. dtcg) then
	        dtcg = dti
	        icg  = i
	        jcg  = jrow
	      endif
	    endif
	  enddo
	enddo
	dtcg = nint(dtcg)
	write (stdout,*)
     & '==> Note: based on a maximum internal gravity wavespeed of'
     &, nint(cgint/100),' m/sec,'
        write (stdout,*)
     & '          the CFL restriction implies a maximum dtuv =',dtcg
     &,' sec'
        write (stdout,*)
     &'          at ','(i,jrow) = (',icg,',',jcg,') on the U grid where'
        write (stdout,*)
     &'          the minimum grid spacing = ~',nint(dxumin/1.e5),' km.'
        write (stdout,*)
     &'          The specified dtuv = ',nint(dtuv),'sec.'
c
c       CFL horizontal advective restriction for tracers
c
        itmin  = imt/2
        jtmin  = jmt/2
        dxtmin = min(dyt(jmt/2),cst(jmt/2)*dxt(imt/2))
	do jrow=jscomp,jecomp
	  do i=2,imt-1
	    if (kmt(i,jrow) .ne. 0) then
	      dxtmin2 = min(dyt(jrow),cst(jrow)*dxt(i))
              if(dxtmin2 .lt. dxtmin) then
                dxtmin = dxtmin2
                itmin  = i
                jtmin  = jrow
              endif
            endif
	  enddo
	enddo
c
        dt_adv = dtts
	c_adv = nint(0.5*dxtmin/dt_adv)
	write (stdout,*)
     & '==> Note: based on a minimum T-cell grid length of approx'
     &, nint(dxtmin/1.e5),' km'
        write (stdout,*)
     & '          at the (i,j) grid point (',itmin,',',jtmin
     &,  ') and a tracer timestep dtts =',nint(dtts),' sec,'
        write (stdout,*)
     & '          the CFL restriction implies a maximum horizontal '
        write (stdout,*)
     & '          advection velocity of ',nint(c_adv),' cm/sec.'
        write (stdout,*)
     & '          advection velocity greater than this is unstable'
c
c       CFL vertical advective restriction for tracers
c
        dt_adv = dtts
	dztmin = dzt(1)
	do k=2,km
	  if (dztmin .gt. dzt(k)) then dztmin = dzt(k)
	enddo
	c_adv = 0.5*dztmin/dt_adv
	write (stdout,*)
     & '==> Note: based on a minimum T-cell grid thickness of approx'
     &, nint(dztmin/1.e2),' m'
        write (stdout,*)
     & '          and a tracer timestep dtts =',nint(dtts),' sec,'
        write (stdout,*)
     & '          the CFL restriction implies a maximum vertical '
        write (stdout,'(a,es10.3,a)')
     & '          advection velocity of ',c_adv,' cm/sec.'
        write (stdout,*)
     & '          advection velocity greater than this is unstable'
	if (dtuv .gt. dtcg) then
	  write (stdout,9000)
     &'==> Error: reduce dtuv or remove this stop condition.'
          errorc = .true.
        endif
	if (abs(nint(dtuv/dtsf)*dtsf-dtuv) .gt. critv) then
	  write (stdout,9000)
     & '==> Warning: for explicit free surface mod(dtuv,dtsf)=0 must  '
     &,'             be satisfied. It is not. Original setting: '
	  write (stdout,*) 'dtuv =  ',dtuv,' sec'
	  write (stdout,*) 'dtsf =  ',dtsf,' sec'
	  write (stdout,9000) '             New setting:               '
          dtsf = dtuv/(int(dtuv/dtsf)+1)
	  write (stdout,*) 'dtsf =  ',dtsf,' sec'
	endif
        if( abs(dtuv/dtsf/2 - nint(dtuv/dtsf/2)) .gt. critv ) then        
          write(*,9000)
     & '==> Error: for explicit free surface dtuv/dtsf must be         '
     & ,'          an even integer to maintain proper centering        '
     & ,'          of barotropic time on the baroclinic time.          '
     & ,'          Reset dtsf and/or dtuv to satisfy this constraint.  '
	  errorc = .true.
        endif 
c
        icg   = 1
	jcg   = 1
	cgmax = 1.0
	dtcgmin = 1.e10
	gridmin = dxu(2)
        do jrow=jscomp,jecomp
	  do i=2,imt-1
	    if (kmu(i,jrow) .gt. 0) then
	      depth = zw(kmu(i,jrow))
	      cgrav = sqrt(grav*depth)
	      gridsp = min(dyu(jrow), csu(jrow)*dxu(i))
	      dtcg = 0.5*(gridsp)/cgrav
	      if (dtcg .lt. dtcgmin) then
	        gridmin = gridsp
                dtcgmin = dtcg
	        cgmax = cgrav
	        icg   = i
	        jcg   = jrow
	      endif
	    endif
	  enddo
	enddo
	dtcgmin = nint(dtcgmin)
	cgmax   = nint(cgmax)
! tp: this code does not work correcly on T90 but does on T3E
! tp: Balaji should look at this
!!        sml = dtcgmin
!!        call mpp_max(sml)
!!        if( pe.EQ.0 )then
!!            if( sml.NE.dtcgmin )
!!     &           call mpp_transmit( icg, 4, NULL_PE, icg, 4, ANY_PE  )
!!        else
!!            if( sml.EQ.dtcgmin )
!!     &           call mpp_transmit( icg, 4,       0, icg, 4, NULL_PE )
!!        endif
        call mpp_max(cgmax)
        call mpp_sync()
        call mpp_min(dtcgmin)
        call mpp_sync()
	write (stdout,9000) '==> Note: for the explicit free surface.'
        write (stdout,*)
     & '          The fastest gravity wave speed =  ',cgmax,'cm/sec'
!!     &,' at i=',icg,' and jrow=',jcg,' on the U grid.'
        write (stdout,*)
     & '          where the grid length = ~',nint(gridmin/1.e5),' km.'
	write (stdout,*)
     & '          For stability, dtsf must be less than ',dtcgmin,' sec'
	write (stdout,*)
     & '          The specified dtuv =  ',dtuv,' sec and dtsf = '
     &, dtsf,' sec. '
	write (stdout,*)
     & '          This implies ',2*nint(dtuv/dtsf),' (due to averaging)'
     &,' barotropic time steps'
        write (stdout,*) 
     & '          per baroclinic time step.'
	if (dtsf .gt. dtcgmin) then
	  write (stdout,9000)'==> Error: for explicit free surface.  '
     &,'          dtsf > time step allowed by fastest wave speed'
	  errorc = .true.
	endif
c
        if (itmb) then
          write (stdout,9000)
     & '==> Warning: "itmb" is set to "true". set it to "false" in    '
     &,'            subsequent runs to prevent the time independent   '
     &,'            basin mask from being written more than once. This'
     &,'            reduces the size of the diagnostic file.          '
	endif
c
        if (itrmb) then
          write (stdout,9000)
     & '==> Warning: "itrmb" is set to "true". set it to "false" in   '
     &,'            subsequent runs to prevent the time independent   '
     &,'            region masks from being written more than once.   '
     &,'            This reduces the size of the diagnostic file.     '
	endif
c
        if (itavg) then
          write (stdout,9000)
     & '==> Warning: "itavg" is set to "true". set it to "false" in   '
     &,'            subsequent runs to prevent the time independent   '
     &,'            region masks from being written more than once.   '
     &,'            This reduces the size of the diagnostic file.     '
	endif
      if (tmbint .gt. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "tmbint" is > zero but the'
     &,'             the "meridional_tracer_budget" option is not on. '
      endif
        if (mod(timavgint,dtts*secday).gt.critv .and. timavgint .gt. c0)
     &    then
          t1 = nint(timavgint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "timavgint" does not contain an integral number of'
     &,'    density time steps "dtts".                                '
          write (stdout,*) '               (changed "timavgint" from '
     &,   timavgint,' days to ', t1,' days to insure this condition)'
          timavgint = t1
        endif
        if (timavgint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: averaging interval "timavgint" = 0.0 implies '
     &,'   no averaging when option "time_averages" is enabled    '
        endif
        if (timavgint .gt. timavgper) then
          write (stdout,9000)
     & '==> Warning: the interval "timavgint" exceeds the averaging'
     &,'             period "timavgper" for option "time_averages" '
        endif
        if (timavgint .lt. timavgper) then
          write (stdout,9000)
     & '==> Warning: averaging period "timavgper" exceeds interval'
     &,'            "timavgint". Setting timavgper = timavgint    '
          timavgper = timavgint
        endif
        if (timavgper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "timavgper" is zero. The '
     &,'            average will be over only one time step!       '
        endif
      if (xbtint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "xbtint"  is > zero but   '
     &,'             the "save_xbts" option is not enabled            '
      endif
      if (dspint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "dspint"  is > zero but'
     &,'             option "diagnostic_surf_height" is not enabled'
      endif
c
      if ((dtuv .ne. dtsf) .or. (dtuv .ne. dtts)) then
        write (stdout,9000)
     & '==> Warning: use of unequal time steps implies the transient  '
     &,'             response may be inaccurate.                      '
      endif
c
c     check for mixing coefficients larger than stability permits
c
      dymin  = dyt(2)
      dxmin  = dxt(2)
      dzmin  = dzt(1)
      xlmax  = dtxcel(1)
      fzimax = fzisop(1)
      do jrow=2,jmtm1
        dymin  = min(dymin,dyt(jrow))
      enddo
      do i=2,imtm1
        dxmin  = min(dxmin,dxt(i))
      enddo
      do k=2,km
        xlmax  = max(xlmax,dtxcel(k))
        dzmin  = min(dzmin,dzt(k))
        fzimax = max(fzimax,fzisop(k))
      enddo
c
      if (xlmax .gt. c1) then
        write (stdout,9000)
     & '==> Warning: use of accelerated time steps implies the        '
     &,'             transient response is unimportant and multiple   '
     &,'             equilibria do not exist. stability tests will    '
     &,'             use "dtts" multiplied by the maximum "dtxcel"    '
      endif
c
c Check that size of viscosity is not too large
c
      write (stdout,9000)
     & '==> Warning: ifdef "shortwave" is enabled but "coupled"       '
     &,'             is not. Is the shortwave component of radiation  '
     &,'             defined ?                                        '
      write (stdout,9000)
     & '==> Warning: "implicitvmix" turns convective adjustment off.  '
     &,'              It will be handled by the vertical mixing       '
c
c     check range of implicit factors
c
c
c     verify that the domain boundary is valid
c
      in = 0
      is = 0
      do i=1,imt
        if (kmt(i,jstask) .ne. 0 .and. jstask .eq. 1) is = i
	if (kmt(i,jetask) .ne. 0 .and. jetask .eq. jmt) in = i
      enddo
      if (is .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along   '
     &,'           the southern boundary.                             '
        write (stdout,*) '           at j=1 and i=',is
      endif 
      if (in .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along   '
     &,'           the northern boundary.                             '
        write (stdout,*) '           at j=jmt and i=',in
      endif 
c
c     verify that each ocean point is at least 2 levels deep
c
      numk = 0
      do jrow=jscomp,jecomp
        do i=1,imt
	  if (kmt(i,jrow) .eq. 1) then
            numk = numk + 1
	    errorc = .true.
	    write (stdout,*)
     &	            ' Error: kmt(',i,',',jrow,') = 1 is not allowed'
	  endif
        enddo
      enddo
      if (numk .ne. 0) then	
        write (stdout,9000)
     & '==> Error: "kmt" must be at least 2 levels deep at all ocean  '
     &,'           points.                                            '
      endif
c

c
c

        if (ALL(kmt(1:imt,jscomp:jecomp) .eq. 0)) then
            write( stdout,'(2a,3i6)' )
     &    '==> Warning: There are no ocean KMT points in'
     &,   '            this parallel domain. Problem domain pn,js,je='
     &,           pn, jscomp, jecomp
        endif
        if (ALL(kmu(1:imt,jscomp:jecomp) .eq. 0)) then
            write( stdout,'(2a,3i6)' )
     &    '==> Warning: There are no ocean KMU points in'
     &,   '            this parallel domain. Problem domain pn,js,je='
     &,           pn, jscomp, jecomp
        endif
        
      write (stdout,'(/,20x,a,/)')
     &         ' E N D    C O N S I S T E N C Y    C H E C K S'
      call mpp_sync()
      if (errorc) then
        write (stdout,'(/a/a)')
     &  'An error flag has been detected in "checks"'
     &,' Use "grep -i error " on this printout file to find ALL errors.'
        call abort()
      endif
c
      return
9000  format (/,(1x,a))
      end
