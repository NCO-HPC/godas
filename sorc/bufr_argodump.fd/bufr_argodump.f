!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM: bufr_argodump
!   PRGMMR: Woollen          ORG: NP20        DATE: 2018-05-15
!
! ABSTRACT: GENERATES tesac dumpfiles from NC031005 (argo) dump files
!
! PROGRAM HISTORY LOG:
! 2018-05-15  JWoollen    ORIGINAL VERSION FOR IMPLEMENTATION
!
! USAGE:
!   INPUT FILES:
!     UNIT 20  - DCOM format argo file            
!
!   OUTPUT FILES:
!     UNIT 50  - ARGO data in DCOM format tesac file  
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   SUBPROGRAMS CALLED:
!     UNIQUE:   - 
!
!     LIBRARY:
!      W3LIB    - W3TAGB W3TAGE ERREXIT
!      BUFRLIB  - OPENBF MAXOUT OPENMB UFBINT
!               - WRITSB WRITCP CLOSBF
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!          >   0 - ABNORMAL RUN
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  WCOSS
!
!$$$
      PROGRAM BUFR_ARGODUMP
 
      PARAMETER(MAX1=11,MAX2=2550)
 
      CHARACTER*8  SUBSET,SUBOUT,SUBINN
      CHARACTER*8  CARR(MAX1,MAX2) 
      EQUIVALENCE  (ARR,CARR)
      LOGICAL      COMPRESS,argos
      REAL*8       ARR(MAX1,MAX2),PP,ZZ
 
      DATA LUBFR/20/
      DATA LUBFO/50/

      DATA SUBINN/'NC031005'/ ! argos dcom message type
      DATA SUBOUT/'NC031002'/ ! tesac dcom message type
 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 
      CALL W3TAGB('bufr_argodump',2013,035,50,'NP20')
 
!  OPEN THE INPUT AND OUTPUT BUFR FILES
!  ------------------------------------

      CALL OPENBF(LUBFR,'IN ',LUBFR)
      CALL OPENBF(LUBFO,'OUT',LUBFR)
 
      call maxout(50000)

!  READ THE MESSAGES FROM THE INPUT FILE
!  -------------------------------------
 
      NRPT_NCO = 0
      NRPT_MOD = 0
 
      DO WHILE(IREADMG(LUBFR,SUBSET,IDATE).EQ.0)
      IF(SUBSET/=SUBINN) CALL BORT('ARGODUMP - INPUT NOT '//SUBINN)
      DO WHILE(IREADSB(LUBFR).EQ.0)
 
      NRPT_NCO = NRPT_NCO + 1
 
!  copy the 031005 elements into the tesac format 
!  ----------------------------------------------
 
      CALL OPENMB(LUBFO,SUBOUT,IDATE)
 
      ! copy x,y,t elements

      CALL UFBINT(LUBFR,ARR,MAX1,MAX2,IRET,'YEAR MNTH DAYS HOUR MINU CLAT CLON CLATH CLONH ')
      CALL UFBINT(LUBFO,ARR,MAX1,IRET,JRET,'YEAR MNTH DAYS HOUR MINU CLAT CLON CLATH CLONH ')

      ! copy receipt time elements

      CALL UFBINT(LUBFR,ARR,MAX1,MAX2,IRET,'RCTS RCYR RCMO RCDY RCHR RCMI IWTEMP ')  
      CALL UFBINT(LUBFO,ARR,MAX1,IRET,JRET,'RCTS RCYR RCMO RCDY RCHR RCMI IWTEMP ')  

      ! copy subsurface elements converting WPRES to DBSS

      CALL UFBINT(LUBFR,ARR,MAX1,MAX2,IRET,'SSTH SALNH WPRES')
      !ARR(3,:)=10.*(101325.+ARR(3,:))/101325. ! convert water pressure to depth below sea surface
      do i=1,iret
        PP = ARR(3,i)/10000.  !  converting ARR(3,:) from Pa to dbar
        ZZ = ((-3.434e-12*PP+1.113e-7)*PP+0.712953)*PP + 14190.7*log(1.0+1.83e-5*PP)
        ZZ = (ZZ/(980.+1.113e-4*PP))*1000.
        ARR(3,i) = ZZ
      enddo
      CALL UFBINT(LUBFO,ARR,MAX1,IRET,JRET,'STMP SALN  DBSS ')

      ! copy the WMOP (numeric) to RPID (char) and left justify RPID string 

      CALL UFBINT(LUBFR,ARR,MAX1,MAX2,IRET,'WMOP')   
      write(carr(1,1),'(i8)') nint(arr(1,1))
      call jstchr(carr(1,1),jret)

      ! append a Q if id fits argo template

      ncar=0; do i=1,8; if(carr(1,1)(i:i)/=' ') ncar=ncar+1; enddo
      if (ncar==7 .and. carr(1,1)(2:2)=='9' .and. carr(1,1)(1:1)<='7' .and. carr(1,1)(1:1)>='1') carr(1,1)='Q'//carr(1,1)
      
      ! write the finished character ID in RPID

      CALL UFBINT(LUBFO,ARR,MAX1,IRET,JRET,'RPID')  

      ! write the tesac templated subset

      call writsb(lubfo)

      NRPT_MOD=NRPT_MOD+1
 
!  END OF READ AND CONVERT LOOPS
!  -----------------------------
 
      ENDDO  !  end of readSB loop
      ENDDO  !  end of readMG loop
 
      CALL CLOSBF(LUBFO)
 
!  FINISHED
!  --------
 
      PRINT'(/,80("-"))'
      PRINT'("Read    ",I7," reports from NC031005 messages: ",a8)', NRPT_NCO
      PRINT'("Wrote   ",I7," reports to   NC031002 messages: ",a8)', NRPT_MOD
      PRINT'(80("-"))'
 
      CALL W3TAGE('bufr_argodump')
      STOP
      END
