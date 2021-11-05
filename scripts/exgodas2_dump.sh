#!/bin/sh
###############################################################################
#                                                                             #
# This script is the data dump for the GODAS model that runs at 18:00 UTC     #
#                                                                             #
#                                                                 April 2003  #
#                                                                             #
###############################################################################

#####################################################################
echo "------------------------------------------------"
echo "JGODAS_DUMP processing:"
echo "Get observed subsurface temperature profiles from /dcom tanks."
echo "- BUFR types bathy and tesac available in monthly tanks are dumped"
echo "using window centered on date of current GODAS run."
echo "- BUFR type dbuoy available in modsbufr.  use dumpmb."
echo "Get recent flux and SST fields for use in future runs."
echo "------------------------------------------------"
echo "History: "
echo " APR 2003 - First implementation of this new script."
echo " JUL 2006 - D.Stokes - Get older buoy data from MODSBUFR.  "
echo "                       This allows earlier cleanup of /com/godas subdirs."
echo "                     - Removed calls to warnmsg.sh and godas_errwarn.sh because"
echo "                       use of prodllsubmit (required to send the email warnings)"
echo "                       is no longer encouraged.  Added warning posts to jlogfile."
echo " FEB 2007 - D.Stokes - Add dump of Jason Altimeter data."
echo "                     - Add buoy_reject list input buoy processing code."
echo " SEP 2007 - D.Stokes - Modified to pick up flux and sst data for the entire period"
echo "                       needed for both standard and extended GODAS."
echo "                     - Since the OISST data is lagged, and day-to-day variation"
echo "                       is relatively small, the last available field is"
echo "                       repeated to pad tail end of SBC window."
echo " MAR 2008 - D.Stokes - Modified to fail (by default) if a daily flux file "
echo "                       is not produced. Prep job will likely fail if 1st or last"
echo "                       day is missing.  Will interpolate if middle day is missing"
echo " FEB 2013 - D.Stokes - Modified to run on WCOSS."
echo " MAY 2019 - H-C LEE  - Modified to run on WCOSS Phase3."
#####################################################################

set -x

cd $DATA
> errfile  # create empty errfile to avoid unnecessary syntax errors if err_chk called following dump scripts

msg="Begin GODAS DATA DUMP PROCESSING for $job on `hostname`"
postmsg "$jlogfile" "$msg"

# define executables and ush scripts

export HOMEbufr=${HOMEbufr:-${HOMEobsproc_dump}}
export USHbufr=${USHbufr:-$HOMEbufr/ush}
export DUMP=${DUMP:-$USHbufr/dumpjb}
export DUMPMB=${DUMPMB:-$USHbufr/dumpmb}
export EXECbufr=${EXECbufr:-$HOMEbufr/exec}
export EXECcombfr=${EXECcombfr:-$EXECbufr/bufr_combfr}

export EXECgodas=${EXECgodas:-$HOMEgodas/exec}
export EXECbathybfr=${EXECbathybfr:-${EXECgodas}/godas_bathybfr}
export EXECbuoybtbfr=${EXECbuoybtbfr:-${EXECgodas}/godas_buoybtbfr}
export EXECdailyflx=${EXECdailyflx:-${EXECgodas}/godas_dailyflx}

export EXECwgrib=${EXECwgrib:-$WGRIB}
export EXECcopygb=${EXECcopygb:-$COPYGB}
export EXECoverdategrb=${EXECoverdategrb:-${OVERDATEGRIB}}

export LIST=${LIST:-${HOMEobsproc_shared_bufr_dumplist}/fix/bufr_dumplist}
export TANK=${TANK:-${DCOMROOT}/prod}

set +x
echo " "
echo "############################################################"
echo "get subsurface data for model ingest window                 "
echo "set below to +/- 2 weeks from model date which is 2 weeks   "
echo "and one day prior to run date                               "
echo "############################################################"
echo " "
set -x

ndays_bathy=29
ndays=$ndays_bathy
# this run at 00Z, so today not included in data window
nmid=`expr $ndays / 2 + 1`
CDATE=`sh finddate.sh $PDY d-$nmid`
CHOUR=12
CTIME=$CDATE$CHOUR
CRAD=`expr  $ndays \* 24 / 2 - 1 `.999

echo $CDATE $ndays_bathy > model_date_info.dump
[ "$SENDCOM" = 'YES' ] && cp model_date_info.dump $COMOUT


if  [ ${SPINTEST:-"NO"} = 'NO' ]; then    # SPINTEST might be YES for development tests.  default is NO

  # tesac and bathy data stored in monthly tanks often delayed, 
  #  so dump entire analysis window for each daily run

# TESTING dcomdev
# The following code section when put into script exgodas1_dump.sh.ecf will convert argo to tesac and cat the tesacs togethe.
# So you end up with a regular looking tesac.ibm file which will flow into the rest of the GODAS as is.
# This should be the only change needed for GODAS.

  export ARGOTANK=$TANK # jsw - post implementation export ARGOTANK=$TANK
  export DUPC=off
  $DUMP  $CTIME  $CRAD  tesac       # jsw - need separate DUMP tesac call for dumping tesacs and argos together
  $DUMP  $CTIME  $CRAD  bathy

# TESTING dcomdev end of block

  dumpstat=$?
  echo status from dump is $dumpstat
  if [ $dumpstat -ne 0 ]; then
    msg="Fatal error -- bathy/tesac bufr dump failed for ${ndays} days surrounding ${CTIME}"
    postmsg "$jlogfile" "$msg"
    export pgm=dumpjb err=$dumpstat; err_chk
  fi
  
  
  set +x
  echo " "
  echo "############################################################"
  echo "extract select info from bathy and tesac bufr files         "
  echo "############################################################"
  echo " "
  set -x
  
  
  for CTYPE in bathy tesac
  do
    export pgm=godas_bathybfr
    . prep_step
  
    export FORT11="$CTYPE.ibm"
    export FORT51="$CTYPE.noqc"
    export FORT61="$CTYPE.noqc.asc"
  
    startmsg
    $EXECbathybfr  >> $pgmout 2> errfile
    export err=$?;err_chk
  
    if test "$SENDCOM" = 'YES'
    then
         cp $CTYPE.noqc $COMOUT/$CTYPE.noqc
    fi
  done
  
  
  set +x
  echo " "
  echo "############################################################"
  echo "get subsurface buoy data from modsbufr"
  echo "############################################################"
  echo " "
  set -x
  
  
  CTYPE=dbuoy
  
  ndays=$ndays_bathy
  dsub=`expr $ndays - 1`
  if [ $dsub -gt 0 ];then
    export StartDY=`sh finddate.sh $PDYm1 d-$dsub`
  else
    export StartDY=$PDYm1
  fi
  
  export EndDY=$PDYm1
  
  [ -z "$MODSONLY" ] && MODSONLY=off
  
  flist=dcom_dbfile.lst
  [ -f $flist ] && rm $flist
  
  if [ $MODSONLY = 'off' ]; then
  # Get what we can from original bufr tanks (save time spent on unnecessary prepmods processing).
  # Original dcom dumped files cannot be combined with modsbufr dumped files because internal tables differ.
  # So, must call godas_buoybtbfr for each source.
    
    CRAD=11.999
    CHOUR=12
    CDATE=$EndDY
  
    until [ $CDATE -lt $StartDY ]; do
      CTIME=$CDATE$CHOUR
  
      $DUMP  $CTIME  $CRAD  $CTYPE
      dumpstat=$?
      echo status from dump is $dumpstat
     #if [ $dumpstat -ne 0 ]; then
      if [ $dumpstat -ne 0 ] || [ ! -f ${TANK}/${CDATE}/${TANKbdir}/${TANKbouy} ]; then
        if [ $CDATE -eq $PDY ]; then   # warning only for TODAYS data.
          msg="WARNING -- Bufr buoy_dump failed for ${CDATE}"
          postmsg "$jlogfile" "$msg"
          export err=0
        else
          msg="Will try modsbufr data for $StartDY thru ${CDATE}"
          postmsg "$jlogfile" "$msg"
          break
        fi
      else 
        mv $CTYPE.ibm $CTYPE.ibm.$CDATE
        echo $CTYPE.ibm.$CDATE >> $flist
      fi
  
      CDATE=`sh finddate.sh $CDATE d-1`
    done
    set +x
    echo " "
    echo "############################################################"
    echo "combine daily bufr files into one file"
    echo "############################################################"
    echo " "
    set -x
  
    if [ -s $flist ]; then
      sort $flist > $flist.sorted
      export pgm=bufr_combfr
      . prep_step
  
      export FORT50="dbuoy.ibm"
  
      startmsg
      $EXECcombfr < $flist.sorted  >> $pgmout 2> errfile
      export err=$?;err_chk
  
      set +x
      echo " "
      echo "############################################################"
      echo "extract select info from combined dcom bufr file            "
      echo "############################################################"
      echo " "
      set -x
  
      export pgm=godas_buoybtbfr
      . prep_step
  
      ndays=$ndays_bathy
      echo $StartDY $ndays > $pgm.input
  
      export FORT11="dbuoy.ibm"
      #--- for NC001103
      export FORT31="$FIXgodas/godas_tao_bufr.txt"
      export FORT32="$FIXgodas/godas_triton_bufr.txt"
      export FORT39="$FIXgodas/godas_rejectbuoy.txt"
      export FORT51="dcom_buoy.noqc"
      export FORT61="dcom_buoy.noqc.asc"
  
      startmsg
      $EXECbuoybtbfr  < $pgm.input  >> $pgmout 2> errfile
      export err=$?;err_chk
  
      cp dcom_buoy.noqc buoy.noqc
    fi
  fi
  
  
  if [ $CDATE -ge $StartDY ]; then   # check if need data from MODSBUFR
    if [ -s $flist ];then
      mbEndDY=$CDATE
    else
      mbEndDY=$EndDY
    fi
  
    $DUMPMB  $StartDY $mbEndDY  dbuoy
    dumpstat=$?
    echo status from dump is $dumpstat
    if [ $dumpstat -ne 0 ]; then
      msg="Fatal error -- MODSBUFR buoy dump failed for ${StartDY} to ${mbEndDY}"
      postmsg "$jlogfile" "$msg"
      export pgm=dumpmb err=$dumpstat; err_chk
    fi
  
  
    set +x
    echo " "
    echo "############################################################"
    echo "combine daily bufr files into one file"
    echo "############################################################"
    echo " "
    set -x
  
    flist=mb_dbfile.lst
    [ -f $flist ] && rm $flist
    export pgm=bufr_combfr
    . prep_step
  
    CDATE=$StartDY
    until [ $CDATE -gt $mbEndDY ]; do
      if [ -s dbuoy.$CDATE ]; then
        echo dbuoy.$CDATE >> $flist
      else
        msg="WARNING:  NO MODSBUFR BUOY DUMP FOR ${CDATE}"
        postmsg "$jlogfile" "$msg"
      fi
      CDATE=`sh finddate.sh $CDATE d+1`
    done
  
    export FORT50="mb_dbuoy.ibm"
  
    startmsg
    $EXECcombfr < $flist  >> $pgmout 2> errfile
    export err=$?;err_chk
  
  
  
    set +x
    echo " "
    echo "############################################################"
    echo "extract select info from dumped bufr file                   "
    echo "############################################################"
    echo " "
    set -x
  
  
    export pgm=godas_buoybtbfr
    . prep_step
  
    ndays=$ndays_bathy
    echo $StartDY $ndays > $pgm.input
  
    export FORT11="mb_dbuoy.ibm"
    #--- for NC001002
    export FORT31="$FIXgodas/godas_tao.txt"
    export FORT32="$FIXgodas/godas_triton.txt"
    export FORT39="$FIXgodas/godas_rejectbuoy.txt"
    export FORT51="mb_buoy.noqc"
    export FORT61="mb_buoy.noqc.asc"
  
    startmsg
    $EXECbuoybtbfr  < $pgm.input  >> $pgmout 2> errfile
    export err=$?;err_chk
  
    cat mb_buoy.noqc >> buoy.noqc
  
  fi
  
  if [ ! -s buoy.noqc ];then
    msg="WARNING.  No Buoy Data for $StartDY to $EndDY"
    postmsg "$jlogfile" "$msg"
    export err=17;err_chk
  fi
  
  if test "$SENDCOM" = 'YES'
  then
      cp buoy.noqc $COMOUT
  fi

else  # SPINTEST...  get original dumped data from production
  for CTYPE in bathy tesac buoy; do                 # SPINTEST
    cp -p $COMIN/$CTYPE.noqc $COMOUT/$CTYPE.noqc    # SPINTEST
  done                                              # SPINTEST
fi    # end of if SPINTEST


## ALTIMETER DATA ##
set +x
echo " "
echo "############################################################"
echo "get altimeter bufr data for entire period."
echo "no extra processing done yet."
echo "############################################################"
echo " "
set -x

ndays=$ndays_bathy
dsub=`expr $ndays - 1`
if [ $dsub -gt 0 ];then
  export StartDY=`sh finddate.sh $PDYm1 d-$dsub`
else
  export StartDY=$PDYm1
fi
export EndDY=$PDYm1

CTYPE=njsnal

if  [ ${SPINTEST:-"NO"} = 'NO' ]; then      # not SPINTEST
$DUMPMB  $StartDY $EndDY  $CTYPE

fi    # end of if NOT SPINTEST

if test "$SENDCOM" = 'YES'
then
  CDATE=$StartDY
  until [ $CDATE -gt $EndDY ]; do
    if  [ ${SPINTEST:-"NO"} = 'NO' ]; then      # not SPINTEST
      DATAFILE=$CTYPE.$CDATE
    else                                        # SPINTEST
      DATAFILE=$COMIN/$CTYPE.$CDATE   # SPINTEST... get data from production run
    fi                                        # end if SPINTEST
    if test -f ${DATAFILE}; then
      cp $DATAFILE $COMOUT
    else
      msg="NO $CTYPE DATA FOR $CDATE."
      postmsg "$jlogfile" "$msg"
    fi
    CDATE=`sh finddate.sh $CDATE d+1`
  done
fi


set +x
echo " "
echo "############################################################"
echo "FLUXES"
echo "############################################################"
echo " "
set -x

#  This job is dumping SBC files to be used by both standard godas
#  and the extended godas.  If the the latter is run for the 13 days
#  after today's standard model run, ndysbc=14.

ndysbc=${ndysbc:-14}
nfields=`expr $ndysbc + 2`

# change this later...
lag=${lag:-15}
echo ${model_start_day:=`sh finddate.sh $PDY d-$lag`}

flux_start_date=`sh finddate.sh $model_start_day d-1`
flux_end_date=`sh finddate.sh $model_start_day d+$ndysbc`

# set list of fields to extract from pressure grib file
cat > fields.list << listEOF
UFLX:sfc:
VFLX:sfc:
SHTFL:sfc:
LHTFL:sfc:
DLWRF:sfc:
ULWRF:sfc:
USWRF:sfc:
DSWRF:sfc:
PRATE:sfc:
listEOF


FLX_DATE=$flux_start_date
until [ $FLX_DATE -gt $flux_end_date ]; do

  FLX_YRMO=`echo $FLX_DATE | cut -c1-6`

  set +x
  echo " "
  echo "############################################################"
  echo "Obtain $FLX_DATE flux fields from pressure grib files "
  echo "  and compute daily averages"
  echo "############################################################"
  echo " "
  set -x
  
  
  export pgm=godas_dailyflx
  . prep_step
  
  for cycle in 00 06 12 18
  do
    fullname=$COMINflux.$FLX_YRMO/flx.ft06.${FLX_DATE}${cycle}.grib  # use this when getting from /com/arkv
#     shortname=flx.ft06.${FLX_DATE}${cycle}.grib  # DWB PATCH
    if [ -s $fullname ];then
      cp $fullname .
#     else
#       cp ../Gflux/$shortname .  # DWB PATCH
#     fi
#     if [ -s $shortname ]; then   # DWB PATCH
      gribfile=`basename $fullname`
      select_name=${gribfile}.select
      $EXECwgrib -s $gribfile | egrep -f fields.list | $EXECwgrib -i -grib $gribfile -o ${select_name}
      err=$?
      if [ $err -eq 0 -a -s $select_name ];then
        $EXECcopygb -g3 -i"1 1" -x $select_name fluxes.$cycle.1deg.grb
      else
        msg="WARNING -- No flux fields extracted from $fullname" 
        postmsg "$jlogfile" "$msg"
        [ -f $select_name ] && mv $select_name $select_name.no_good
      fi
    else
      msg="WARNING -- Atmospheric Analysis file not available for ${FLX_DATE}${cycle}."
      postmsg "$jlogfile" "$msg"
    fi
  
    case $cycle in
      00) export FORT21=fluxes.$cycle.1deg.grb;;
      06) export FORT22=fluxes.$cycle.1deg.grb;;
      12) export FORT23=fluxes.$cycle.1deg.grb;;
      18) export FORT24=fluxes.$cycle.1deg.grb;;
       *) echo "SKIP UNEXPECTED CYCLE: $cycle";
           continue;;
      esac
  
  done
  
  export FORT61=taux.$FLX_DATE
  export FORT62=tauy.$FLX_DATE
  export FORT63=senflx.$FLX_DATE
  export FORT64=latflx.$FLX_DATE
  export FORT65=prate.$FLX_DATE
  export FORT66=netlwflx.$FLX_DATE
  export FORT67=netswflx.$FLX_DATE
  
  echo $FLX_DATE > $pgm.input
  
  startmsg
  $EXECdailyflx  < $pgm.input  >> $pgmout 2> errfile
  export err=$?;err_chk
  
  
  if test "$SENDCOM" = 'YES'
  then
    for flux in taux tauy senflx latflx prate netlwflx netswflx
    do
      if [ -s $flux.$FLX_DATE ]; then
        cp -p $flux.$FLX_DATE $COMOUT
      else
        msg="WARNING -- Check if any flx files for $FLX_DATE available in $COMINflux.$FLX_YRMO"
        postmsg "$jlogfile" "$msg"
        if  [ ${FAIL_IF_NO_FLUX:-"YES"} = 'YES' ]; then    #  allows option to continue if missing a day of fluxes
          err=8; err_chk
        fi
      fi
    done
  fi

  FLX_DATE=`sh finddate.sh $FLX_DATE d+1`
done



set +x
echo " "
echo "############################################################"
echo "SST"
echo ""
echo "  Obtain SST fields from cdas archive or sst directory."
echo "  OISST is a 7-day analysis run for days PDYm7 to PDYm1 and stamped PDY."
echo "  CDAS arkv has them stamped as originally generated."
echo "  Change internal date and label file based on middle"
echo "      of 7-day analysis period. (PDYm4)"
echo "############################################################"
echo " "
set -x


SSTGODAS_DATE=$flux_start_date
until [ $SSTGODAS_DATE -gt $flux_end_date ]; do


  SSTOI_DATE=`sh finddate.sh $SSTGODAS_DATE d+4`
  SSTOI_YRMO=`echo $SSTOI_DATE | cut -c1-6`

# Get the 12Z file.  OISST runs only once per day and gets copied to arkv 4 times.  
# The 12Z cdas file is typically the first of the day to have the "PDY" SST
  fullname=$COMINflux.$SSTOI_YRMO/sstgrb${SSTOI_DATE}12  # use this when getting from /com/arkv
# if cdas2 file missing (eg, if delayed), try original sstoi directory
  [ -s $fullname ] || fullname=$COMINsst.$SSTOI_DATE/sstoi_grb
  [ ${SPINTEST:-"NO"} = 'YES' -a $SSTOI_DATE -gt $PDYm1 ] && fullname=sstoi_grb.last   # simulate realtime conditions

# If desired file does not exist, repeat last field.  We need boundary conditions up to $PDYm1, 
#   but SST is lagged.  So, this is likely to be needed for last few days. 
  [ -s $fullname ] || fullname=sstoi_grb.last
  if [ -s $fullname ];then
    echo cp $fullname ./sstoi_grb
    cp $fullname ./sstoi_grb
    err=$?
    if [ $err -eq 0 ]; then
#    overdate.grib.sh ${SSTGODAS_DATE}00 sstoi_grb
# wrapper script overdate.grib.sh not ported to WCOSS by NCO.  replace with following snippet -dcs ###
      rm -f fort.11 fort.51
      ln -s sstoi_grb fort.11
      ln -s sstoi_grb.ovrdate fort.51
      echo ${SSTGODAS_DATE}00 | $EXECoverdategrb
      [ $? -ne 0 ] && err=`expr $err + 2`
      mv sstoi_grb.ovrdate sstoi_grb
## end overdate.grib.sh replacement snippet ###
      [ $? -ne 0 ] && err=`expr $err + 2`
    fi
    if [ $err -ne 0 ]; then
      msg="WARNING -- Problem finding or overdating SST field"
      postmsg "$jlogfile" "$msg"
    else
      if test "$SENDCOM" = 'YES'
      then
        cp sstoi_grb $COMOUT/sst.grb.$SSTGODAS_DATE
        if [ $fullname != sstoi_grb.last ]; then
          cp sstoi_grb sstoi_grb.last
          LAST_SSTOI_DATE=$SSTOI_DATE
          LAST_SSTGODAS_DATE=$SSTGODAS_DATE
        else
          set +x
          echo " "
          echo "########################################################################################################"
          echo "SST field from $SSTOI_DATE not found.  Repeat SST field from $LAST_SSTOI_DATE centered on $LAST_SSTGODAS_DATE"
          echo "########################################################################################################"
          echo " "
          set -x
        fi
      fi
    fi
  else
    msg="WARNING -- SST Analysis file not available."
    postmsg "$jlogfile" "$msg"
  fi
  SSTGODAS_DATE=`sh finddate.sh $SSTGODAS_DATE d+1`
done



#####################################################################
# GOOD RUN
set +x
echo "************** $job COMPLETED NORMALLY ON THE IBM SP"
echo "************** $job COMPLETED NORMALLY ON THE IBM SP"
echo "************** $job COMPLETED NORMALLY ON THE IBM SP"
set -x
#####################################################################

msg="$job HAS COMPLETED NORMALLY!"
postmsg "$jlogfile" "$msg"

############## END OF SCRIPT #######################

