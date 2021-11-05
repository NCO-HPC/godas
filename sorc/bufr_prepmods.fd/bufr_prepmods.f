C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM: BUFR_PREPMODS
C   PRGMMR: WHITING          ORG: NP20        DATE: 2013-02-04
C
C ABSTRACT: GENERATES MODS (MARINE OCEAN DATA-HANDLING SYSTEM)
C   BUFR ARCHIVE FILES FROM INPUT DCOM-FORMATTED DATA DUMPS.  ONLY
C   SELECT INPUT MESSAGE TYPES DEEMED TO BE OF INTEREST TO MARINE
C   APPLICATIONS ARE CONSIDERED AS VALID.  
C
C PROGRAM HISTORY LOG:
C 2004-05-07  JWHITING    ORIGINAL VERSION FOR IMPLEMENTATION
C 2004-06-03  JWOOLLEN    UPDATED  VERSION FOR IMPLEMENTATION
C 2004-06-10  JWhiting    Changed OBNAME 001,002 to point to 
C                           DBUOYSUB rather than DBUOY.
C                         Added reports read/written output.
C 2004-06-22  JWhiting    Added report counting to SSMIPN subr.
C 2005-05-06  JWhiting    Added inline wrcmps() (fr J.Ator lib)
C 2005-05-12  JWhiting    Added REHU (rel humidity); MAX1=11
c                         seqnum needed???
c 2007-01-30  JWhiting    removed inline wrcmps() so as to be 
c                           compatible w/ BUFRLIB upgrade.
c 2007-10-04  JWhiting    enabled MODS type envsal (NC031109)
c 2009-03-11  JWhiting    disabled compression; this is due to a 
c                       problem when running code recompiled on 
c                       Cirrus (debugged to a BUFRLIB memory 
c                       bounds issue);  use of compression should 
c                       be reconsidered after implementation of a 
c                       BUFRLIB fix.
c 2011-03-08  JWhiting  development for restricted SHIPS tanking;
c                        added SHIPSU processing
c                        added restriction mnemonics for SHIPS
c 2011-03-28  JWhiting  enabled support for shpall dump group
C 2013-02-04  JWhiting  Port to WCOSS (linux) platforms:
C                        rm'd unused variables (ADATE, MI);
C                        Updated doc block and comments;
C                        Ready for implementation.
C
C USAGE:
C   INPUT FILES:
C     UNIT 20  - BUFR DATA FILE (DCOM FORMATTED)
C     UNIT 21  - MODS BUFR MNEMONIC TABLE
C
C   OUTPUT FILES:
C     UNIT 50  - BUFR DATA FILE (MODS FORMATTED)
C     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:   - MNEMONICS OBNAME RADDATE SSMIPN TIDEG  CLONLAT
C
C     LIBRARY:
C      W3LIB    - W3TAGB W3TAGE ERREXIT
C      BUFRLIB  - OPENBF MAXOUT OPENMB UFBINT
C               - WRITSB WRITCP CLOSBF
C
C   EXIT STATES:
C     COND =   0 - SUCCESSFUL RUN
C          >   0 - ABNORMAL RUN
C
C REMARKS:
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  WCOSS
C
C$$$
      PROGRAM BUFR_PREPMODS
 
      PARAMETER(MAX1=11,MAX2=2550)
 
      COMMON /MNEMONIC/ NEMOS(10),NEMS,COMPRESS
 
      CHARACTER*80 NEMOS,DBNEM,DBNEMBF
      CHARACTER*8  SUBSET,NCOSET,NCOSET0,OBNAME
      LOGICAL      COMPRESS
      REAL*8       DATES(5),ARR(MAX1,MAX2),DTCC(2)
      REAL*8       stsl(3,max2),curr(3,max2)
      REAL*8       SUBS,CLONH,CLATH,bmiss
 
      DATA LUBFR/20/
      DATA LUNDX/21/
      DATA LUBFO/50/
 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
 
      CALL W3TAGB('BUFR_PREPMODS',2013,035,50,'NP20')
 
      print *
      print * ,'-> Welcome to BUFR_PREPMODS of 2/04/2013',
     &           ' (WCOSS implementation)'

      print *

      bmiss=10e10; call setbmiss(bmiss)

      open(6,recl=150)

C  OPEN THE INPUT AND OUTPUT BUFR FILES
C  ------------------------------------

      CALL OPENBF(LUBFR,'IN ',LUBFR)
      CALL OPENBF(LUBFO,'OUT',LUNDX)
 
      call maxout(50000)

C  READ THE MESSAGES FROM THE INPUT FILE
C  -------------------------------------
 
      NRPT_NCO = 0
      NRPT_MOD = 0
 
1     DO WHILE(IREADMG(LUBFR,NCOSET,IDATE).EQ.0)
      read(ncoset(3:5),'(i3)') mtp
      read(ncoset(6:8),'(i3)') mst
      subset = obname(mtp,mst)

      if(subset.eq.'        ') then
         print*,ncoset,' >>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<'
         print*,ncoset,' >>>>>has no modsbufr equivalent<<<<<'
         print*,ncoset,' >>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<'
         print*,'prepmods error - unknown input subset encountered'
         CALL W3TAGE('BUFR_PREPMODS')
         call errexit(99)
c        call exit(99)
      endif

C  INITIALIZE THE MNNEMONIC TEMPLATE FOR THIS MESSAGE TYPE
C  -------------------------------------------------------
 
      CALL MNEMONICS(SUBSET,NEMOS,NEMS,COMPRESS)
 
      if(nems.eq.0) then
         print*,subset,' >>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<'
         print*,subset,' >>>>>has no modsbufr definitions<<<<<'
         print*,subset,' >>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<'
         print*,'prepmods error - modsbufr subset missing definitions'
         CALL W3TAGE('BUFR_PREPMODS')
         call errexit(99)
c        call exit(99)
      endif

C  READ THE SUBSETS IN EACH MESSAGE
C  --------------------------------
 
2     DO WHILE(IREADSB(LUBFR).EQ.0)
 
      NRPT_NCO = NRPT_NCO + 1
 
C  CASES WHERE REPORTS ARE SPLIT UP IN MODSBUFR
C  --------------------------------------------
 
      IF(SUBSET.EQ.'TIDEG') THEN
         NRPT_TID = 0
         CALL TIDEG(LUBFR,LUBFO,SUBSET,IDATE,NRPT_TID)
         NRPT_MOD = NRPT_MOD + NRPT_TID
         GOTO 2
      ENDIF
 
      IF(SUBSET.EQ.'SSMIPN') THEN
         NRPT_SSMI = 0
         CALL SSMIPN(LUBFR,LUBFO,SUBSET,IDATE,NRPT_SSMI)
         NRPT_MOD = NRPT_MOD + NRPT_SSMI
         GOTO 2
      ENDIF
 
C  CASES WHERE REPORTS ARE JUST COPIES INTO MODSBUFR
C  -------------------------------------------------
 
      CALL OPENMB(LUBFO,SUBSET,IDATE)
 
      CALL CLONLAT(LUBFR,CLONH,CLATH)
      CALL UFBINT(LUBFO,CLONH,1,1,IRET,'CLONH')
      CALL UFBINT(LUBFO,CLATH,1,1,IRET,'CLATH')
      NRPT_MOD = NRPT_MOD + 1
 
      DO N=1,NEMS
      CALL UFBINT(LUBFR,ARR,MAX1,MAX2,LEV2,NEMOS(N))
      if(lev2.gt.0) CALL UFBINT(LUBFO,ARR,MAX1,LEV2,IRET,NEMOS(N))
      ENDDO

c  dbouy type has two possible sources of sub-surface data
c  -------------------------------------------------------

      DBNEM='DBSS STMP SALN DROC SPOC'
      DBNEMBF='DBSS SST1 SALN DROC SPOC'

      if(subset=='DBUOY'.and.ncoset=='NC001002') then
         CALL UFBINT(LUBFR,ARR,MAX1,MAX2,LEV2,DBNEM)
         IF(LEV2>0) CALL UFBINT(LUBFO,ARR,MAX1,LEV2,IRET,DBNEM)
      endif

      if(subset=='DBUOY'.and.ncoset=='NC001103') then
         call ufbseq(lubfr,stsl,3,max2,lev0,'BBYSTSL')
         call ufbseq(lubfr,curr,3,max2,lev1,'BBYCURR')
         lev2=lev0+lev1
         if(lev2>0) then
            arr=bmiss
            do lev=1,lev0
            arr(1,lev)=stsl(1,lev)
            arr(2,lev)=stsl(2,lev)
            arr(3,lev)=stsl(3,lev)
            enddo
            do lev=1,lev1
            arr(1,lev0+lev)=curr(1,lev)
            arr(4,lev0+lev)=curr(2,lev)
            arr(5,lev0+lev)=curr(3,lev)
            enddo
            if(lev2>bmiss) then
             do lev=1,lev2
             write(6,'(5f16.2)') arr(1:5,lev)  
             enddo
             pause; print*
            endif
c           CALL UFBINT(LUBFO,ARR,MAX1,LEV2,IRET,DBNEM)
            CALL UFBINT(LUBFO,ARR,MAX1,LEV0,IRET,DBNEMBF)
         endif
      endif

c  tesac has a weird dtcc situation
c  --------------------------------
 
      IF(SUBSET.EQ.'TESAC') THEN 
         CALL UFBINT(LUBFR,DTCC,2,1,IRET,'DTCC')
         IF(DTCC(1).NE.19) IDTCC = 1
         IF(DTCC(2).NE.19) IDTCC = 2
         CALL UFBINT(LUBFO,DTCC(IDTCC),1,1,IRET,'DTCC')
      ENDIF
 
C  END OF READ AND CONVERT LOOPS
C  -----------------------------
 
      CALL WRITSB(LUBFO)

      ENDDO  !  end of readSB loop
      NCOSET0 = NCOSET
      ENDDO  !  end of readMG loop
 
      CALL CLOSBF(LUBFO)
 
C  FINISHED
C  --------
 
      PRINT'(/,80("-"))'
      PRINT'("Read    ",I7," reports from BUFR messages with ",
     $       "Table A entry: ",A8)', NRPT_NCO,NCOSET0
      PRINT'("Wrote   ",I7," reports to MODS file type ",
     $          a8)', NRPT_MOD,subset
      PRINT'(80("-"))'
 
      CALL W3TAGE('BUFR_PREPMODS')
      STOP
      END
