!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  BUOYBFR
!   PRGMMR: Diane C. Stokes  ORG: NP23        DATE: 2013-01-28
!
! ABSTRACT:  Extract subsurface temperature data from buoy bufr file
!   for use in Global Ocean Data Assimilation System
!
! PROGRAM HISTORY LOG:
! 2003-06-02  Diane C. Stokes
! 2005-08-10  Diane C. Stokes - modified to process multi-day file 
!                               and modsbufr mnemonics (implemented 2006-08-08)
! 2006-09-12  Diane C. Stokes - added check to reject obs from suspect buoys
! 2013-01-28  Diane C. Stokes - minor modifications for WCOSS transition. 
!                               improve missing checks.
!
! USAGE:
!   INPUT FILES:
!     UNIT 11  - BUOY DATA IN BUFR
!     UNIT 31  - LIST OF TAO BUOY IDS  (includes ATLAS and PIRATA)
!     UNIT 32  - LIST OF TRITON BUOY IDS
!     UNIT 39  - LIST OF BUOY IDS TO REJECT
!
!   OUTPUT FILES:
!     UNIT 51  - BUOY DATA IN IEEE
!     UNIT 06  - UNIT 6 (STANDARD PRINTFILE)
!
!   SUBPROGRAMS CALLED FROM PROGRAM: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     UNIQUE:    - none
!     LIBRARY:
!       W3NCO    - w3fs26, iw3jdn, w3tagb, w3tage, errexit
!        BUFR    - openbf, datelen, ireadns, ufbint, closbf
!
!   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
!     UNIQUE:    - none
!     LIBRARY:
!       W3NCO    - w3fs26, iw3jdn, w3tagb, w3tage, errexit
!        BUFR    - openbf, datelen, ufbint, ireadns, closebf
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!     COND =   7 - ERROR READING TAO IDS
!     COND =   8 - ERROR READING TRITON IDS
!     COND =   9 - ERROR READING REJECT LIST
!     COND =  17 - BUOY OB DATE DOES NOT MATCH EXPECTED DATE
!
! REMARKS AND IMPORTANT LOCAL VARIABLES: 
!                   Reports with less than 3 levels are excluded
!
! ATTRIBUTES:  (LIST MACHINES FOR WHICH CODE IS USED AND CHECKED OUT)
!
!   MACHINE:  IBM iDataPlex
!   LANGUAGE: F90
!
!
!$$$

c get subsurface data from bufr'd surface marine data 
c  output in format used for vsam
c  output with modsbufr table (commented out)

      integer, parameter :: lubfw=71,lubfmt=41
      real(8) getbmiss, fmiss
      real(8) vec1(10),vsub(2,255),dval,wtmp,depth,veco(8),vsubo(2,255)
      real(8) xtemp
      real xz(255),xt(255),xtk(255)
      integer iw3jdn
      character*80 string1, stringo
      character*8 subset,blank,ctemp,cdate,subseto
      character*8 cstnid,rawrpt(255)
      character sid*2,dtyp*2,ukey*1,qkey*1
      logical lbad
      data blank/'        '/
C  ukey is U for update, R for replace
      data sid/'NW'/,ukey/'U'/,qkey/' '/   
*     data sid/'NW'/,ukey/'R'/,qkey/' '/
      equivalence(xtemp,ctemp)
      character (len=8),allocatable:: ctao(:), ctriton(:), creject(:)

      call w3tagb('GODAS_BUOYBTBFR',2013,0032,0000,'NP23')

      fmiss=getbmiss()

      subseto = 'DBUOY   '   ! For mbufr output -- Not used yet.
      string1 = 'RPID YEAR MNTH DAYS HOUR MINU CLATH CLONH CLAT CLON'
      stringo = 'RPID YEAR MNTH DAYS HOUR MINU CLATH CLONH'

      read(5,*)cdate,ndays
      read(cdate,'(i4,2i2)')iyrst,imst,idst
      ijdst=iw3jdn(iyrst,imst,idst)
      ijdend = ijdst + ndays - 1
      call w3fs26(ijdend,iyrend,imend,idend,idwk,idyr)

      print '(a,i5.4,i3.2,i3.2)','Start date:',iyrst,imst,idst
      print '(a,i5.4,i3.2,i3.2)','  End date:',iyrend,imend,idend

c read in TAO and TRITON buoy ids
      ntao=0
      do
        read(31,*,end=70)
         ntao= ntao+1
      enddo
   70 continue
      allocate(ctao(ntao))
      rewind 31
      do kid=1, ntao
        read(31,*,err=970)ctao(kid)
      enddo

      ntriton=0
      do
        read(32,*,end=80)
         ntriton=ntriton+1
      enddo
   80 continue
      allocate(ctriton(ntriton))
      rewind 32
      do kid=1, ntriton
        read(32,*,err=980)ctriton(kid)
      enddo

c read in buoy ids to reject
      nrejid=0
      open(39,status='old',action='read',iostat=ios)
      if(ios.eq.0)then
        do
          read(39,*,end=90)
          nrejid= nrejid+1
        enddo
   90   continue
        print '(i0,a)',nrejid,' BUOYS ON REJECT LIST'
        if(nrejid.gt.0)then
          allocate(creject(nrejid))
          rewind 39
          do kid=1, nrejid
            read(39,*,err=990)creject(kid)
          enddo
        endif   ! end-if-any-reject-ids
      endif     ! end-if-reject-file-opened

c  open output bufr file, using modsbufr table
*dcs* call openbf(lubfw,'OUT',lubfmt)

      ict = 0
      nct = 0
      klt3 = 0
      krejob = 0
      lubfr=11
      call openbf( lubfr, 'IN', lubfr )
      call datelen(10)
      DO WHILE(IREADNS(LUBFR,SUBSET,IDATE).EQ.0)
        if (subset.ne.'DBUOY   '.and.subset.ne.'NC001002'
     &      .and.subset.ne.'NC001103')cycle
        ict = ict + 1
        lbad=.false.
        if (subset.eq.'DBUOY   '.or.subset.eq.'NC001002') then       
          call ufbint(lubfr, vsub,  2, 255, nret, 'DBSS STMP')
        elseif (subset.eq.'NC001103') then
          call ufbint(lubfr, vsub,  2, 255, nret, 'DBSS SST1')
        endif
        if(nret.le.0) go to 800

        np=0
*       if(nret.ne.0)then
          do k=1,nret
            depth=vsub(1,k)
            wtmp=vsub(2,k)
            if(depth.ne.fmiss.and.wtmp.ne.fmiss)then
              np=np+1
              xz(np)=depth
              vsubo(1,np)=depth
*             if(wtmp.gt.200)then
                xt(np)=wtmp-273.15
                vsubo(2,np)=wtmp
*             else
*               xt(np)=wtmp
*             endif
            endif
          enddo
*       endif
*       if(np.eq.0)go to 800
        if(np.lt.3)then
           klt3=klt3+1
           go to 800
        endif

        
        call ufbint(lubfr, vec1, 10, 1, nret, string1)
        if(nret.ne.1) then
           print '(a,i0,a)','NRET =',nret,'   get next msg'
           go to 800
        endif
        xtemp=vec1(1)
        if (ibfms(xtemp).eq.1)then
          cstnid=blank
        else
          cstnid=ctemp
        endif
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
        ijd=iw3jdn(iyr,imo,ida)
        if(ijd.lt.ijdst.or.ijd.gt.ijdend)then
          print '(a,3i5)','Date out of range:',iyr,imo,ida
          go to 800
        endif
        iyoc=mod(iyr,100)  !  2000 -> 0, as req by Ming
*       icen = (iyr + 99)/100  !  century  1999->20, 2000->20, 2001->21
        ihr=nint(hr)

        imin=nint(xmin)
        jdate=imo*1000000+ida*10000+ihr*100+imin
        kdate=iyr*1000000+imo*10000+ida*100+ihr
        if(kdate.ne.idate)then
          print '(a,i0,1x,i0)', 'ERROR?  kdate ne idate:',kdate, idate
          call w3tage('GODAS_BUOYBTBFR')
          call errexit(17)
C         call exit(17)
        endif
        

        read(cstnid(3:3),'(i1)',err=199)idr

        if(nrejid.gt.0)then
c         Check if buoy on reject list...
          do kid=1,nrejid
            if(cstnid.eq.creject(kid)) then
              print*,' AUTO REJECT OB FROM BUOY ',cstnid
              krejob=krejob+1
              go to 800
            endif
          enddo
        endif

        if(idr.ge.0.and.idr.le.4)then
c         moored buoy.  probably TAO or TRITON, check id against list
c         TAO first
          do kid=1,ntao
            if(cstnid.eq.ctao(kid))then
              dtyp='TA'
              go to 200
            endif
          enddo
c         now TRITON
          do kid=1,ntriton
            if(cstnid.eq.ctriton(kid))then
              dtyp='TR'
              go to 200
            endif
          enddo
c  TAO and TRITON ids exhausted.  
          print*,'Moored buoy not in TAO or TRITON list.'
          print '(a,a,2(1x,f7.2))','dtyp set to BU for id and loc: ',
     1                 cstnid,clat,clon
          dtyp='BU'
        else
          dtyp='BU'
        endif
        go to 200
  199   continue
          print*,' non-numeric in id.  ',cstnid,'. dtyp set to BU'
          dtyp='BU'
  200   continue

c
*dcs*   call openmb(lubfw, subseto, kdate)
*dcs*   veco(1:6)=vec1(1:6)
*dcs*   veco(7)=clon
*dcs*   veco(8)=clat
*dcs*   call ufbint(lubfw, veco, 8, 1, nret, stringo)
*dcs*   if(nret.ne.1) print*,'Not all values written for:',veco
*dcs*   call ufbint(lubfw, vsubo, 2, np, nret, string3)
*dcs*   if(nret.ne.np) print*,'Not all values written for:',
*dcs*1      veco,'/',vsubo
*dcs*   call writsb(lubfw)

c no ukey for this file
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
      call closbf(lubfw)
      print*, 'ict=',ict,'    nct=',nct
      print*, 'num obs skipped for rejected id:',krejob
      print*, 'num reports skipped (lt 3 levs):',klt3
      call w3tage('GODAS_BUOYBTBFR')
      call errexit(0)
C     call exit(0)
  970 print*,'Premature end of file reading TAO ids'
      call w3tage('GODAS_BUOYBTBFR')
      call errexit(7)
C     call exit(7)
  980 print*,'Premature end of file reading TRITON ids'
      call w3tage('GODAS_BUOYBTBFR')
      call errexit(8)
C     call exit(8)
  990 print*,'Error reading file with buoy ids to reject'
      call w3tage('GODAS_BUOYBTBFR')
      call errexit(9)
C     call exit(9)
      end
