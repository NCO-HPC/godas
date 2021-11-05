C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE TIDEG(LUBFR,LUBFO,SUBSET,IDATE,NRPT)
 
      COMMON /MNEMONIC/ NEMOS(10),NEMS,COMPRESS
 
      CHARACTER*80 NEMOS
      CHARACTER*8  SUBSET
      LOGICAL      COMPRESS,NOTIME
      REAL*8       BMISS,CLONH,CLATH,ADATE,DATE(5)
      REAL*8       TIMEP(2),TIDES(2,255),ARR(10,255)
 
      DATA BMISS/10E10/
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  READ THE TIME PERIOD INFORMATION AND HIGH RES LON/LAT
C  -----------------------------------------------------
 
      CALL UFBINT(LUBFR,TIMEP,2,  1,IRET,'TPMI  TIMI')
      CALL UFBINT(LUBFR,TIDES,2,255,NTID,'TIDER TERC')
      CALL CLONLAT(LUBFR,CLONH,CLATH)
 
C  CHECK IF THERE IS A TIME FRAME HERE AND SETUP FIRST IDATE
C  ---------------------------------------------------------
 
      NOTIME = TIMEP(1).GE.BMISS.AND.TIMEP(2).GE.BMISS

      IF(NOTIME.AND.NTID.GT.1) THEN
         PRINT*,'>>>>>TIDEG - NO TIMES FOR MULTIPLE READINGS<<<<<'
      ENDIF

      IF(NOTIME) THEN
         ADATE = IDATE
         DHR = 0
      ELSE
         ADATE = IDATE
         DHR = (TIMEP(1)+TIMEP(2))/60.
         CALL RADDATE(ADATE,DHR,ADATE)
         DHR = TIMEP(2)/60.
      ENDIF
 
C  CREATE TIDE GAUGE REPORTS
C  -------------------------
 
      DO M=1,MAX(NTID,1)
 
      IF(TIDES(1,M).LT.BMISS.OR.TIDES(2,M).LT.BMISS.OR.M.EQ.1) THEN
         JDATE = ADATE
         DATE(1) = JDATE/1000000
         DATE(2) = MOD(JDATE/10000,100)
         DATE(3) = MOD(JDATE/100  ,100)
         DATE(4) = MOD(JDATE      ,100)
         DATE(5) = (ADATE-JDATE)*60.          

         CALL OPENMB(LUBFO,SUBSET,IDATE)
 
         CALL UFBINT(LUBFO,CLONH,1,1,IRET,'CLONH')
         CALL UFBINT(LUBFO,CLATH,1,1,IRET,'CLATH')

         CALL UFBINT(LUBFO,DATE,5,1,IRET,'YEAR MNTH DAYS HOUR MINU')
         nrpt = nrpt + 1
 
         DO N=1,NEMS
         CALL UFBINT(LUBFR,ARR,10,255,NRET,NEMOS(N))
         IF(NRET.GT.0) CALL UFBINT(LUBFO,ARR,10,NRET,IRET,NEMOS(N))
         ENDDO
 
         CALL UFBINT(LUBFO,TIDES(1,M),2,1,IRET,'TIDER TERC')
 
         CALL RADDATE(ADATE,DHR,ADATE)  ! APPLY TIME INCREMENT
      ENDIF
 
c Disable compression (due to problem porting to Cirrus; fails on wrcmps)
c  (this should be revisited after update to Cirrus BLib)
ccc   IF(.NOT.COMPRESS) CALL WRITSB(LUBFO)
ccc   IF(     COMPRESS) CALL WRITCP(LUBFO)
      CALL WRITSB(LUBFO)
 
 
      ENDDO
 
      RETURN
      END
