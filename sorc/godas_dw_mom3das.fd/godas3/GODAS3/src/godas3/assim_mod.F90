module assim_mod
#ifdef assm_dta
use time_manager_mod
!
!   Assimilation data for MOM
!
!   Basic paramters
!
!   kass   = number of vertical levels in assimilation
!   kass2  = kass if only temperature is corrected; = 2*kass if salinity is also
!   mobs   = maximum obs in a single latitude row
!   maxits = maximum iterations in oi analysis scheme
!   npits  = number of iterations used in lsmth to define [e]
!   nlnmx  = max size of nolen
!
# if defined cor_sal || defined asm_sal
integer, parameter :: kass = 35, kass2 = 2*kass
# else
integer, parameter :: kass = 35, kass2 = kass
# endif
integer, parameter :: mobs = 100000, maxits = 3, npits = 110, nlnmx = 20000
integer, save :: imtka
!
!   namelist parameters (setassim)
!
!   rassint = frequency in days at which the assimilation is triggered
!             (see "casimsw.h")
!   nassdur = number of repeat assimilations at each triggering 
!             (see "asm_x_mod")
!   acoef = global scaling for model error variance
!   vcvs  = model vertical covariance scale
!   hrscl = model horizontal covariance scale
!   b2    = (hrscl * pi / 360.0)**2
!   aeval = constant used in correction for topography
!
real, save :: acoef, vcvs, hrscl, aeval, b2
integer, save :: mitrDB, pnvw
!
!   1D vertical decomposition.  The 1D horizontal decomposition of MOM3 is
!   used, but a transpose to a vertical decomposition is done in "lpsmthr"
!
!   max_tasks_v      = maximum number of processors
!   num_processors_v = requested number of processors
!   pnv              = this processor number (from 1 to num_processors_v)
!
!   kscomp3(n)  = uppermost computed level for processor "n"
!   kecomp3(n)  = lowermost computed level for processor "n"
!
!   kscomp  = uppermost computed level for this processor (pnv)
!   kecomp  = lowermost computed level for this processor (pnv)
!
integer, parameter :: max_tasks_v = 2048
integer, save :: pnv, num_processors_v
integer, save, dimension(max_tasks_v) :: kscomp3, kecomp3
integer, save :: kscomp, kecomp
!
!   temporary dummy file for output
!
integer, save :: stddum
!
!   kma(imt,jmt) = a global copy of kmt allocated in setassim
!
integer, save, dimension(:,:), allocatable :: kma
!
!   The observations
!
!   For now there is no flexibility.  All data sets are organized in
!   weekly periods, all have the same number of periods and all have
!   the same starting and ending dates.
!
!     nwkobs = number of weekly observation periods on disk.
!     obstmp = time stamps marking the beginning/end of each weekly period
!     o_time = array of time_type marking the beginning/end of each period
!     obs_start_time = start time of observations
!     obs_end_time   = end time of observations
!     dayasm = position of current model time within observations
!     toprd   = relative times of centers of observation periods (Wed, noon)
!     oprd    = length observation period. (7 days, constant)
!     aobs   = array for facillitating disk IO of observations
!
#if defined pentad || defined opr_day
integer, parameter :: nwkobs = 5
#else
integer, parameter :: nwkobs = 10
#endif
character(len=32), save, dimension(0:nwkobs) :: obstmp
type(time_type), save, dimension(nwkobs+1) :: o_time
type(time_type), save :: obs_start_time, obs_end_time
real, save, dimension(nwkobs) :: toprd
real, save :: oprd, dayasm
!  real(kind=4), dimension(6,mobs) :: aobs
real, save, dimension(6,mobs) :: aobs
character(len=128), save :: opt_dasm
!
#ifdef asm_sst
!
!     nwksst  = number of weeks in assimilation window for SST
!     ndxsst  = index for SST for assimilation time manager
!     ipsstd  = array of pointers to weeks that are in assimilation window
!     iosstao = unit number for word addressable disk
!     itobs   = array holding number of observations per row in itobs(j,1,*)
!               and sum of observations in all prior rows in itobs(j,2,*)
!
integer, parameter :: nwksst = 3
integer, save :: ndxsst, iosstao
integer, save, dimension(:,:,:), allocatable :: itobs
integer, save, dimension(nwksst) :: ipsstd
!
#endif
#ifdef asm_tmp
!
!     nwktmp  = number of weeks in assimilation window for temperature
!     ndxtmp  = index for temperature for assimilation time manager
!     iptmpd  = array of pointers to weeks that are in assimilation window
!     iotmpao = unit number for word addressable disk
!     jtobs   = array holding number of observations per row in jtobs(j,1,*)
!               and sum of observations in all prior rows in jtobs(j,2,*)
!
integer, parameter :: nwktmp = 5
integer, save :: ndxtmp, iotmpao
integer, save, dimension(:,:,:), allocatable :: jtobs
integer, save, dimension(nwktmp) :: iptmpd
!
#endif
#ifdef asm_sal
!
!     nwksal  = number of weeks in assimilation window for salinity
!     ndxsal  = index for salinity for assimilation time manager
!     ipsald  = array of pointers to weeks that are in assimilation window
!     iosalao = unit number for word addressable disk
!     ktobs   = array holding number of observations per row in ktobs(j,1,*)
!               and sum of observations in all prior rows in ktobs(j,2,*)
!
integer, parameter :: nwksal = 5
integer, save :: ndxsal, iosalao
integer, save, dimension(:,:,:), allocatable :: ktobs
integer, save, dimension(nwksal) :: ipsald
!
#endif
#ifdef asm_ssh
!
!     nwkssh  = number of weeks in assimilation window for SSH
!     ndxssh  = index for SSH for assimilation time manager
!     ipsshd  = array of pointers to weeks that are in assimilation window
!     iosshao = unit number for word addressable disk
!     ltobs   = array holding number of observations per row in ltobs(j,1,*)
!               and sum of observations in all prior rows in ltobs(j,2,*)
!
integer, parameter :: nwkssh = 5
integer, save :: ndxssh, iosshao
integer, save, dimension(:,:,:), allocatable :: ltobs
integer, save, dimension(nwkssh) :: ipsshd
!
!     sshc = climatological mean of model sea surface height
!     dssh = difference between current SSH and climatological mean
!
real, save, dimension(:,:), allocatable :: sshc, dssh
!
!  Coefficients for computing dynamic height as a linear function
!   of temperature or of temperature and salinity.  The model layer
!   thicknesses are also taken into account.  See lDensCoefTS.F and 
!   mkLDZCoef.c for computing these coefficients.
!
!     cdnz  = coefficients for temperature
real, save, dimension(:), allocatable :: cdnz
!
# if defined cor_sal || defined asm_sal
!     cdnzs = coefficients for salinity
real, save, dimension(:), allocatable :: cdnzs
!
# endif
#endif
!
!   Arrays used in specifying the background error covariance
!
real, save, dimension(kass2,kass2) :: covsr
!
!   cvn  = local vertical covariance matrix for temperature
!   vtmp = vertical variance for temperature
!
real, save, dimension(kass,kass) :: cvn
# ifdef fix_vv
real, save, dimension(:,:,:), allocatable :: vtmp
# else
real, save, dimension(:,:,:,:), allocatable :: vtmp
# endif
real, dimension(:,:), allocatable :: bufik
!
!   nvvrec = number of variance records (12 - periodic, 14 - monthly)
!   vvprd  = (.true./.false.) to indicate data (periodic/not periodic)
!   vstamp = time stamps marking the end of data records
!   vprec  = period lengths for the data records (days per month)
!   vvrec  = times marking the centers the data records
!   indxvv = interpolation index
!   mthdvv = interpolation method = (0,1,2,3)
!   dayvv  = model time in days after start of variance data
!
!   iprvvd = index pointing to data on disk prior to current model time
!   inxtvd = index pointing to data on disk after current model time
!   iprvvm = index pointing to data in memory prior to current model time
!   inxtvm = index pointing to data in memory after current model time
!   wprvv  = interpolation factor for prior data, so
!               data = wprvv*data(iprvvm) + (1-wprvv)*data(inxtvm)
!
!   iotvv  = unit for time dependent monthly temperature variance data
# if defined cor_sal || defined asm_sal
!   iosvv  = unit for time dependent monthly salinity variance data
# endif
!
integer, parameter :: maxvv = 20
!
# if defined cor_sal || defined asm_sal
integer, save :: nvvrec, iprvvd, inxtvd, iprvvm, inxtvm, indxvv, mthdvv, iotvv
integer, save :: iosvv
# else
integer, save :: nvvrec, iprvvd, inxtvd, iprvvm, inxtvm, indxvv, mthdvv, iotvv
# endif
real, save :: dayvv, wprvv, vprec(maxvv), vvrec(maxvv)
logical, save :: vvprd, rdvv
character(len=32), save :: vstamp0, vstamp(maxvv)
type(time_type), save :: vv_start_time, vv_end_time, v_time

# if defined cor_sal || defined asm_sal
!
!   cvnsalt = local vertical covariance matrix for salinity
!   vsal    = vertical variance for salinity
!
real, save, dimension(kass,kass) :: cvnsalt
#  ifdef fix_vv
real, save, dimension(:,:,:), allocatable :: vsal
#  else
real, save, dimension(:,:,:,:), allocatable :: vsal
#  endif
# endif
!      #endif
!
!   Arrays used in the minimizing the cost function. 
!
!   In general, arrays used by the Laplace smoother will be allocated in
!   the meridional direction with the global dimension, jmt, because that
!   routine is decomposed in the vertical direction.  The exception is s2h,
!   which facilitates the multiplication by the vertical covariance matrix
!   in a horizontal decomposition. Other arrays that are not part of the
!   Laplace smoother will be allocated for a horizontal decomposition.
!
!     elipt controls anisotropy of background error covariance (jmt)
!     wgta is normalization for Laplace smoother (lpwghts and lpsmthr) (jmt)
!
real, save, dimension(:), allocatable :: wgta, elipt
!
!     wew, wno, wso are N, E, S and W coefficients, allocated as (imt,jmt)
!     scra, temp, dpth are arrays allocated as (imt,jmt)
!
real, save, dimension(:,:), allocatable :: wew, wno, wso, scra, temp, dpth
!
!    There are two s2 arrays used in lpsmthr
!     s2v   = vertical slab version (imt,km,js:je)
!     s2h   = horizontal slab version (imt,jmt,ks:ke)
!     s2buf = buffer for transposing between s2v and s2h
!
real, save, dimension(:,:,:), allocatable :: s2v, s2h, s2buf
integer, save :: ns2b
!
!    The following arrays represent T, d, e, f, g, and h in Derber 
!    and Rosati (1989) and will be dimensioned (imt,kass2,jstask:jetask)
!    They are tagged with "_cg" to identify them with the congugate
!    gradient algorithm of the assimilation.
!
real, save, dimension(:,:,:), allocatable :: t_cg, d_cg, e_cg, f_cg, g_cg, h_cg
real, save, dimension(:,:), allocatable :: r_cg
!
!   Arrays for internal management of observations
!
!    The following arrays are for all types of observations taken together 
!     nobs  contains the number of observations per latitude row (js:je)
!     nsobs contains the running sum of observations per latitude row (js:je)
!    The following arrays are for one latitude row at a time (mobs)
!     rtm   contains the observation times relative to the current model time
!     val   contains the observation values
!     aip   contains the i positions
!     ajp   contains the j positions
!     akp   contains the k positions
!     aerr  contains the estimated error values
!     a111, a211, a221, a121 are linear interpolation coefficients
!     ils   are integer truncations of aip
!
integer, save, dimension(:), allocatable :: nobs, nsobs
real, save, dimension(:), allocatable :: rtm, val, aip, ajp, akp, aerr, ov
real, save, dimension(:), allocatable :: a111, a211, a221, a121
integer, save, dimension(:), allocatable :: ils
!
!    The following arrays are storage for all rows in a process ((js:je)*mobs)
!     srtm  contains the observation times relative to the current model time
!     sval  contains the observation values
!     sip   contains the i positions
!     sjp   contains the j positions
!     skp   contains the k positions
!     serr  contains the estimated error values
!     lass1, lass2 control the alternating storage in these arrays
!
real, save, dimension(:), allocatable :: srtm, sval, sip, sjp, skp, serr
integer, save :: lass1, lass2
!
!   dtlim = limit on model-obs temp difference beyond which invF is tapered
!
real, save :: dtlim = 5.0
# ifdef asm_sal
!
!   dslim = limit on model-obs saltdifference beyond which invF is tapered
!
real, save :: dslim = 3.0
# endif
!
# ifdef exRes
!
!   rsdys  = save residuals within rsdys of the model time
!   ioextr = file unit for saving extracted temperature residuals
!   ioexsr = file unit for saving extracted salinity residuals
!   ioexer = file unit for saving extracted eta (SSH) residuals
!   tresd  = holds filename for temperature residuals
!   sresd  = holds filename for salinity residuals
!   eresd  = holds filename for eta (SSH) residuals
!
real, save :: rsdys = 0.25
integer, save :: ioextr, ioexsr, ioexer
character(len=10), save :: tresd, sresd, eresd
!
# endif
#else
real :: assm_dum
#endif
end module assim_mod
