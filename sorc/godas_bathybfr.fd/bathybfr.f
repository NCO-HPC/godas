!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM: GODAS_BATHYBFR
!   PRGMMR: STOKES           ORG: NP23        DATE: 2003-06-13
!
! ABSTRACT:  Extract subsurface temperature data from input bufr file
!   for use in Global Ocean Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2003-06-02  Diane C. Stokes
! 2004-01-20  Diane C. Stokes - modified for upcoming change to 
!                               bufrtab.031.
!                                 - removed reference to {BTOCN}.
!                                 - allow up to 65535 levels for 
!                                   subsurface data.
! 2007-02-16  Diane C. Stokes - discard obs from inst type 852, currently 
!                               known to have pressure offset errors.
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - BATHY OR TESAC DATA IN BUFR
!
!   OUTPUT FILES:
!     UNIT 51  - BATHY OR TESAC DATA IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - none
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage,errexit
!        BUFR    - openbf, ufbint,closbf
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - none
!     LIBRARY:
!       W3LIB    - w3tagb, w3tage,errexit
!        BUFR    - openbf, ufbint,closbf
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!
! REMARKS AND IMPORTANT LOCAL VARIABLES: 
!                   Reports with less than 3 levels are excluded
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM SP
!   LANGUAGE: F90
!
!
!$$$
! get subsurface data from bufr'd surface marine data 
!
      real(8) vec1(10),vsub(2,65535),dval,wtmp,xtemp
      real xz(65535),xt(65535)
      integer iw3jdn
      character*80 string1,string2,string3,string4
      character*8 subset,blank,ctemp,cdate
      character*8 cstnid
      character*50 card50
      character sid*2,dtyp*2,qkey*1
      logical lbad,wrongdate
      data fmiss/10e10/ , eps/10e-2/ , blank/'        '/
      data sid/'NW'/,qkey/' '/   
      equivalence(xtemp,ctemp)

      call w3tagb('GODAS_BATHYBFR',2003,0164,0164,'NP23')

      string1 = ' RPID YEAR MNTH DAYS HOUR MINU CLATH CLONH CLAT CLON'
      string2 = ' SST1 '
      string3 = ' DBSS STMP '
*     string4 = ' {BTOCN} '

*     read(5,*)cdate,ndays
*     read(cdate,'(i4,2i2)')iyrst,imst,idst
*     ijdst=iw3jdn(iyrst,imst,idst)
*     ijdend = ijdst + ndays - 1
*     call w3fs26(ijdend,iyrend,imend,idend,idwk,idyr)

*     print*,'Start date:',iyrst,imst,idst
*     print*,'  End date:',iyrend,imend,idend
*     print*

      ict = 0
      nct = 0
      n852 = 0
      mxlevs=0
      klt3=0
      lubfr=11
      call openbf( lubfr, 'IN', lubfr )
      DO WHILE(IREADNS(LUBFR,SUBSET,IDATE).EQ.0)
        if (subset.eq.'NC031001'.or.subset.eq.'BATHY   ')then
          dtyp='BA'
        else if (subset.eq.'NC031002'.or.subset.eq.'TESAC   ')then
          dtyp='TE'
        else
          cycle
        endif
        ict = ict + 1
        lbad=.false.


! hard code toss of instrument type 852, currently known to have pressure offset errors (2/16/2007)
        call ufbint(lubfr, dval, 1, 1, nret, 'IWTEMP')
        if(nret.gt.0.and.nint(dval).eq.852)then
          n852=n852+1
          cycle
        endif

*       call ufbint(lubfr, dval, 1, 1, nret, string4)
*       if(nret.eq.0)go to 800
*       if(dval.eq.0)go to 800

        np=0
*       call ufbint(lubfr, wtmp, 1, 1, nret, string2)
*       if(nret.eq.0.or.abs(wtmp-fmiss).le.eps)then
*         np=0
*       else
*         np=1
*         xz(np)=0.0
*         xt(np)=wtmp-273.15
*       endif
       
        call ufbint(lubfr, vsub, 2, 65535, nret, string3)
        if(nret.eq.0)go to 800
*       if(nret.ne.0)then
          if(nret.gt.mxlevs)mxlevs=nret
          do k=1,nret
            depth=vsub(1,k)
            wtmp=vsub(2,k)
            if(abs(depth-fmiss).gt.eps.and.abs(wtmp-fmiss).gt.eps)then
              np=np+1
              xz(np)=depth
              xt(np)=wtmp-273.15
            endif
          enddo
*       endif
*       if(np.eq.0)go to 800
        if(np.lt.3)then
           klt3=klt3+1
           go to 800
        endif
        
        call ufbint(lubfr, vec1, 10, 1, nret, string1)
        xtemp=vec1(1)
        yr   =vec1(2)
        xmo  =vec1(3)
        day  =vec1(4)
        hr   =vec1(5)
        xmin =vec1(6)
        if (vec1(7) .lt. fmiss) then
           clat=vec1(7)        ! CLATH used
        else if (vec1(9) .lt. fmiss) then
           clat=vec1(9)        ! CLAT used
        else
           print*,'NO VALID LAT.  SKIP THIS REPORT'
           go to 800
        endif
        if (vec1(8) .lt. fmiss) then
           clon=vec1(8)        ! CLONH used
        else if (vec1(10) .lt. fmiss) then
           clon=vec1(10)       ! CLON used
        else
           print*,'NO VALID LON.  SKIP THIS REPORT'
           go to 800
        endif

        if(clon.lt.0)clon=clon+360.
        iyr=nint(yr)
        imo=nint(xmo)
        ida=nint(day)
*       wrongdate=(iyr.ne.jyr.or.
*    .             imo.ne.jmo.or.
*    .             ida.ne.jda)
*       if(wrongdate)then
*         print '(a,3i5)','Wrong date:',iyr,imo,ida
*         go to 800
*       endif
        ihr=nint(hr)
        imin=nint(xmin)
        jdate=imo*1000000+ida*10000+ihr*100+imin

        if(abs(xtemp-fmiss).le.eps)then
          cstnid=blank
        else
          cstnid=ctemp
        endif

! no ukey for this file
      write (51) iyr,jdate,cstnid,sid,dtyp,qkey,clat,clon,np
     1          ,(xz(k),xt(k),k=1,np)

      write(61,215)iyr,jdate,cstnid,clat,clon,np,sid,dtyp,qkey
  215 format(i4,1x,i10,1x,a8,1x,f8.2,f9.2,i3,a2,a2,a1)
        write(61,'(15f8.2)') (xz(k),k=1,np)
        write(61,'(15f8.2)') (xt(k),k=1,np)

        nct=nct+1
  800 continue
      ENDDO
      call closbf(lubfr)
  900 continue
      print*, 'ict=',ict,'    nct=',nct
      print*, 'max num levs found:',mxlevs
      print*, 'num reports skipped (lt 3 levs):',klt3
      print*, 'num reports from recorder type 852 skipped:',n852
      call w3tage('GODAS_BATHYBFR')
      call errexit(0)
!     call exit(0)
      end
