C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM: GODAS_DAILYFLX
C   PRGMMR: STOKES           ORG: NP23        DATE: 2013-01-28
C
C ABSTRACT: Extract select fields from GDAS files and process as 
C   needed for flux fields used by the Ocean Data Assimlation System 
C   for the Climate Forecast Suite.
C
C PROGRAM HISTORY LOG:
C   03-03-13  Diane C. Stokes
C   03-02-04  Diane C. Stokes - include PRATE as possible input 
C                                (replace APCP)
C  2013-01=28 Diane C. Stokes - Changes to run on WCOSS.
C
C USAGE:
C   INPUT FILES:
C     UNIT 21  - subset of gdas 00Z 6hr fcst pgrb file or
C                subset of avn 18Z 12hr fcst pgrb file
C     UNIT 22  - subset of gdas 06Z 6hr fcst pgrb file or
C                subset of avn 00Z 12hr fcst pgrb file
C     UNIT 23  - subset of gdas 12Z 6hr fcst pgrb file or
C                subset of avn 06Z 12hr fcst pgrb file
C     UNIT 24  - subset of gdas 18Z 6hr fcst pgrb file or
C                subset of avn 12Z 12hr fcst pgrb file
C     PARM     - UNIT 5 (STANDARD READ)
C
C   OUTPUT FILES:
C     UNIT 61  - zonal wind stress      (-UFLX)
C     UNIT 62  - meridional wind stress (-VFLX)
C     UNIT 63  - sensible heat flux     (-SHTFL)
C     UNIT 64  - latent heat flux       (-LHTFL)
C     UNIT 65  - precip rate            (APCP/(6*60*60))
C     UNIT 66  - net longwave flux      (DLWRF-ULWRF)
C     UNIT 67  - net shortwave flux     (DSWRF-USWRF)
C     UNIT 06  - standard printout
C
C   SUBPROGRAMS CALLED FROM PROGRAM: (CALLED FROM ANYWHERE IN CODES)
C     UNIQUE:    
C     LIBRARY:
C       W3NCO    - getgb, w3tagb, w3tage, errexit
C       BALIB    - baopenr
C
C   SUBPROGRAMS CALLED FROM MAIN: (LIST ALL CALLED FROM MAIN)
C     UNIQUE:    
C     LIBRARY:
C       W3NCO    - getgb, errexit
C       BALIB    - baopenr
C
C   EXIT STATES:
C     COND =   0 - SUCCESSFUL RUN
C          =   6 - ENVIRONMENT VARIABLE FOR INPUT FILE DOES NOT EXIST
C          =   7 - INPUT FILENAME TOO LONG FOR VARIABLE LENGTH
C          =   8 - NON-SPECIFIC ERROR GETTING INPUT FILENAME
C          =   9 - UNEXPECTED STATUS GETTING INPUT FILENAME
C
C REMARKS: AVN f12 used if GDAS f06 missing.
C
C ATTRIBUTES:
C   MACHINE:  IBM SP
C   LANGUAGE: f90
C
C$$$

      implicit none
      integer,parameter :: imax=360, jmax=181
      integer,parameter :: ijmax=imax*jmax
      integer,parameter :: nfld=10,nfldo=7
      integer jpds(200),jgds(200),ivar(nfld),ilev(nfld)
      integer kpds(200),kgds(200)
      real rinc(5)
      integer idat(8),jdat(8)
      integer lugb,length,iestatus
      integer ifld,jskp,nij,kmsg,iretb,iretg,kret(nfld)
      integer jyr,jyr2,jmo,jda,jhr,jcen,icyc,kcyc
      integer jyr_avn,jyr2_avn,jmo_avn,jda_avn,jhr_avn,jcen_avn

      logical lb(ijmax)
      real, dimension (imax,jmax,nfld) :: flds
      real, dimension (imax,jmax) :: fld1,wkarray,dailyavg
      character(4),dimension(4)::cycle=(/'t00z','t06z','t12z','t18z'/)
      character(6)::clun
      character(80)::cgrib
      type grbinfo
         character (len=5) :: cvar
         integer :: ivar
         integer :: ilev
         character (len=30) :: cdesc
      end type grbinfo
      type (grbinfo), dimension(nfld) :: infld
      type outinfo
         real, dimension(imax,jmax) :: dailysum=0.
         integer :: ncycle=0, kflg=0
         integer :: iuout
      end type outinfo
      type (outinfo) taux,tauy,senflx,latflx,prate,netlwflx,netswflx


      call w3tagb('GODAS_DAILYFLX',2013,0032,0000,'NP23')

      read(5,'(i4,i2,i2)')jyr,jmo,jda
      print '(3(x,i0))',jyr,jmo,jda
      jyr2 = mod(jyr-1,100)+1
      jcen  = (jyr + 99)/100

      infld(1)=grbinfo('UFLX', 124,1,'Zonal momentum flux [N/m^2]')
      infld(2)=grbinfo('VFLX', 125,1,'Meridional momentum flux [N/m^2]')
      infld(3)=grbinfo('SHTFL',122,1,'Sensible heat flux [W/m^2]')
      infld(4)=grbinfo('LHTFL',121,1,'Latent heat flux [W/m^2]')
      infld(5)=grbinfo('DLWRF',205,1,'Downward long wave flux [W/m^2]')
      infld(6)=grbinfo('ULWRF',212,1,'Upward long wave flux [W/m^2]')
      infld(7)=grbinfo('USWRF',211,1,'Upward solar rad flux [W/m^2]')
      infld(8)=grbinfo('DSWRF',204,1,'Downward solar rad flux [W/m^2]')
      infld(9)=grbinfo('PRATE',  59,1,'Precipitation Rate [kg/m^2/s]')
      infld(10)=grbinfo('APCP',  61,1,'Total precipitation [kg/m^2]')

      taux%iuout=61
      tauy%iuout=62
      senflx%iuout=63
      latflx%iuout=64
      prate%iuout=65
      netlwflx%iuout=66
      netswflx%iuout=67

      do 700 icyc=1,4

        jhr=(icyc-1)*6
        lugb=20+icyc

        write(clun,'(a,i2.2)')"FORT",lugb
        call get_environment_variable(clun, cgrib, length, iestatus)
        select case(iestatus)
          case(0)
            continue
          case(1)
            print '(a,x,i0)','NO FILENAME ASSOCIATED WITH FORT',lugb
            print*,'SKIP CYCLE ',cycle(icyc)
            go to 700
          case(-1)
            print '(a,a,a,x,i0,a)','env variable ',trim(clun),
     1        ' is set to string of', length,
     2        ' characters which does not fit in cgrib.'
            call errexit(7)
!           call exit(7)
          case(3)
            print*,'non-specific error(s) from GET_ENVIRONMENT_VARIABLE'
            call errexit(8)
!           call exit(8)
          case default
            print*,'unexpected status from GET_ENVIRONMENT_VARIABLE'
            call errexit(9)
!           call exit(9)
        end select

        if(len_trim(cgrib).le.0)then
          print '(a,x,i0)','NO FILENAME ASSOCIATED WITH FORT',lugb
          print*,'SKIP CYCLE ',cycle(icyc)
          go to 700
        endif
        print '(a,x,i0,a,a)','name for lugb=',lugb,' is ',trim(cgrib)
        call baopenr(lugb,cgrib,iretb)

        if(iretb.eq.0)then
          print*,'OPENED FILE FOR ',cycle(icyc),':  ',trim(cgrib)
c  check that file is for correct date/time...
          jgds=-1
          jpds=-1
          jskp=-1
          call getgb(lugb,0,ijmax,jskp,jpds,jgds,nij,kmsg,kpds,kgds,
     1        lb,fld1,iretg)
          if(iretg.eq.0)then
            print '(a,22(1x,i0))','kgds: ',kgds(1:22)
            if(kpds(16).eq.10)then
              if(kpds(14).eq.6)then
                print*,'We have f06'
                kcyc=icyc
              else if(kpds(14).eq.12)then
                print*,'We have f12'
                kcyc=icyc-1
                if(kcyc.eq.0)kcyc=4
              endif
            else
              if(kpds(14).eq.0.and.kpds(15).eq.6)then
                print*,'We have f06'
                kcyc=icyc
              else if(kpds(14).eq.6.and.kpds(15).eq.12)then
                print*,'We have f12'
                kcyc=icyc-1
                if(kcyc.eq.0)kcyc=4
              else
                print 205,'KPDS:',kpds(1:25)
                print*,'UNEXECTED TIME RANGE IN ',trim(cgrib)
                print*,'SKIPPING THIS CYCLE'
                go to 700
              endif
            endif
          else
            print*,' GETGB FAILED FOR ',cycle(icyc),':  ',trim(cgrib)
            print '(a,x,i0)',' return code from getgb is:',iretg
            print 205,'KPDS:',kpds(1:25)
            print*,'SKIPPING THIS CYCLE'
            go to 700
          endif
        else
          print*,' BAOPEN FAILED FOR ',cycle(icyc),':  ',trim(cgrib)
          print*,'SKIPPING THIS CYCLE'
          go to 700
        endif

  100   continue

        flds=-9999.            !initialize array

        do ifld=1,10
          if(ifld.eq.1)then
            jskp=-1
          else
            jskp=0
          endif

          print*,' Getting ',infld(ifld)%cvar,':  ',infld(ifld)%cdesc
          jgds=-1
          jpds=-1
          jpds(5)=infld(ifld)%ivar
          jpds(6)=infld(ifld)%ilev
          jpds(7)=0
          if(kcyc.eq.icyc)then
            jpds(8)=jyr2
            jpds(9)=jmo
            jpds(10)=jda
            jpds(11)=jhr
            jpds(14)=0
            jpds(15)=6
            jpds(21)=jcen
          else
            rinc=(/0.,-6.,0.,0.,0./)
            idat=(/jyr,jmo,jda,0,jhr,0,0,0/)
            call w3movdat(rinc,idat,jdat)
            jyr_avn=jdat(1)
            jmo_avn=jdat(2)
            jda_avn=jdat(3)
            jhr_avn=jdat(5)
            jyr2_avn=mod(jyr_avn-1,100)+1
            jcen_avn=(jyr_avn+99)/100
            jpds(8)=jyr2_avn
            jpds(9)=jmo_avn
            jpds(10)=jda_avn
            jpds(11)=(kcyc-1)*6
            jpds(14)=6
            jpds(15)=12
            jpds(21)=jcen_avn
          endif
          call getgb(lugb,0,ijmax,jskp,jpds,jgds,nij,kmsg,kpds,kgds,
     1       lb,fld1,iretg)
          kret(ifld)=iretg
          if(iretg.eq.0) then
            flds(:,:,ifld)=fld1
            print 205,'KPDS:',kpds(1:25)
            print*,fld1(180,90)
          else
            print '(a,x,i0)',
     1           'GETGB FAILED FOR JPDS.  return code=',iretg
            print 205,'JPDS:',jpds(1:25)
            print*
            print 205,'KPDS:',kpds(1:25)
          endif
      enddo

        
c get output fields, reversing latitudes from S-N to N-S
      if(kret(1).eq.0) then
        wkarray             = -flds(:,jmax:1:-1,1) 
        taux%dailysum       = taux%dailysum + wkarray
        taux%ncycle         = taux%ncycle   + 1
        taux%kflg          = ibset(taux%kflg,icyc-1)
        print*,'daily sum for taux(180,90)=',taux%dailysum(180,90)
        print '(a,i0)','ncycle for taux=',taux%ncycle
      endif
      if(kret(2).eq.0) then
        wkarray             = -flds(:,jmax:1:-1,2) 
        tauy%dailysum       = tauy%dailysum + wkarray
        tauy%ncycle         = tauy%ncycle   + 1
        tauy%kflg          = ibset(tauy%kflg,icyc-1)
      endif
      if(kret(3).eq.0) then
        wkarray             = -flds(:,jmax:1:-1,3) 
        senflx%dailysum     = senflx%dailysum + wkarray
        senflx%ncycle       = senflx%ncycle   + 1
        senflx%kflg        = ibset(senflx%kflg,icyc-1)
      endif
      if(kret(4).eq.0) then
        wkarray             = -flds(:,jmax:1:-1,4) 
        latflx%dailysum     = latflx%dailysum + wkarray
        latflx%ncycle       = latflx%ncycle   + 1
        latflx%kflg        = ibset(latflx%kflg,icyc-1)
      endif
c see if PRATE available.  If not, convert total precip to precip rate.
      if(kret(9).eq.0) then
        wkarray             = flds(:,jmax:1:-1,9)
        prate%dailysum      = prate%dailysum + wkarray
        prate%ncycle        = prate%ncycle   + 1
        prate%kflg         = ibset(prate%kflg,icyc-1)
        print '(a,i0)','kflg for prate=',prate%kflg     !dcs
      else if(kret(10).eq.0) then
        wkarray             = flds(:,jmax:1:-1,10)/float(6*60*60)
        prate%dailysum      = prate%dailysum + wkarray
        prate%ncycle        = prate%ncycle   + 1
        prate%kflg         = ibset(prate%kflg,icyc-1)
        print '(a,i0)','kflg for prate=',prate%kflg     !dcs
      endif
      if(kret(5).eq.0.and.kret(6).eq.0)then 
        wkarray             = flds(:,jmax:1:-1,5)-flds(:,jmax:1:-1,6)
        netlwflx%dailysum   = netlwflx%dailysum + wkarray
        netlwflx%ncycle     = netlwflx%ncycle   + 1
        netlwflx%kflg      = ibset(netlwflx%kflg,icyc-1)
      endif
      if(kret(7).eq.0.and.kret(8).eq.0)then 
        wkarray             = flds(:,jmax:1:-1,8)-flds(:,jmax:1:-1,7)
        netswflx%dailysum   = netswflx%dailysum + wkarray
        netswflx%ncycle     = netswflx%ncycle   + 1
        netswflx%kflg      = ibset(netswflx%kflg,icyc-1)
      endif

  700 continue
      

      if((taux%ncycle).gt.0)then
        dailyavg=taux%dailysum/float(taux%ncycle)
        write(taux%iuout)jyr,jmo,jda,taux%kflg,dailyavg
      endif
      if(tauy%ncycle.gt.0)then
        dailyavg=tauy%dailysum/float(tauy%ncycle)
        write(tauy%iuout)jyr,jmo,jda,tauy%kflg,dailyavg
      endif
      if(senflx%ncycle.gt.0)then
        dailyavg=senflx%dailysum/float(senflx%ncycle)
        write(senflx%iuout)jyr,jmo,jda,senflx%kflg,dailyavg
      endif
      if(latflx%ncycle.gt.0)then
        dailyavg=latflx%dailysum/float(latflx%ncycle)
        write(latflx%iuout)jyr,jmo,jda,latflx%kflg,dailyavg
      endif
      if(prate%ncycle.gt.0)then
        dailyavg=prate%dailysum/float(prate%ncycle)
        write(prate%iuout)jyr,jmo,jda,prate%kflg,dailyavg
      endif
      if(netlwflx%ncycle.gt.0)then
        dailyavg=netlwflx%dailysum/float(netlwflx%ncycle)
        write(netlwflx%iuout)jyr,jmo,jda,netlwflx%kflg,dailyavg
      endif
      if(netswflx%ncycle.gt.0)then
        dailyavg=netswflx%dailysum/float(netswflx%ncycle)
        write(netswflx%iuout)jyr,jmo,jda,netlwflx%kflg,dailyavg
      endif

  205 format(a,25i5)
      call w3tage('GODAS_DAILYFLX')
      call errexit (0)
C     call exit (0)
      end
