c====================== include file "casimsw.h" ========================
c
c     holds switching information for assimilation
c     defaults and namelist settings are handled in setassim.F
c     it is "included" in switch.h
c
c     rassim  = true when within 1/2 time step of a specified interval
c               else ... false
c     rassint = interval in days for triggering a new assimilation
c
      logical rassim
      common /switcr/ rassint
      common /switcl/ rassim
c
