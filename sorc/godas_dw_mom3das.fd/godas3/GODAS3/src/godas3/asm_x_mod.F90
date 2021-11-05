module asm_x_mod
#ifdef assm_dta
!
!   External assimilation data for MOM
!
!   It is through these parameters and arrays that the model and the 
!   assimilation communicate
!
!   rassint = frequency in days at which the assimilation is triggered
!             (see "casimsw.h")
!   nassdur = number of repeat assimilations at each triggering
!   nassim  = counter for repeat assimilations
!
integer, save :: nassdur, nassim
!
!   arex = first guess at the beginning of an assimilation
!   ares = field correction at completion of an assimilation
!
real, dimension(:,:,:), allocatable, save :: arex, ares
# ifdef asm_ssh
real, dimension(:,:), allocatable, save :: etax
# endif
!
#else
real :: ares
#endif
end module asm_x_mod
