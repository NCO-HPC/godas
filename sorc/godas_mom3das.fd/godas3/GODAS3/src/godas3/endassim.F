      subroutine endassim
#ifdef assm_dta
!
!-----------------------------------------------------------------------
!    Close the assimilation files
!----------------------------------------------------------------------
!
      use assim_mod
!
      integer :: error
!
#ifdef asm_sst
!
      call relunit (iosstao)
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
#endif
#ifdef asm_tmp
!
      call relunit (iotmpao)
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
#endif
#ifdef asm_sal
!
      call relunit (iosalao)
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
#endif
#ifdef asm_ssh
!
      call relunit (iosshao)
!
      call MPI_BARRIER(MPI_COMM_WORLD, error)
!
#endif
#ifdef exRes
!
      call relunit (ioextr)
# ifdef asm_sal
      call relunit (ioexsr)
!
# endif
# ifdef asm_ssh                                                                 
      call relunit (ioexer)
!
# endif
#endif
!----------------------------------------------------------------------
#else
      real :: setadum
#endif
      end subroutine endassim
