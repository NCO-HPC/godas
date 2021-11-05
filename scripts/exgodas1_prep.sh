#!/bin/sh
###############################################################################
#                                                                             #
# This script is the prep step for the GODAS model.                           #
# Data and fields are pulled from com directories.                            #
# QC is done on subsurface temp profiles.                                     #
#                                                                             #
#                                                                   May 2003  #
#                                                                             #
###############################################################################

#####################################################################
echo "------------------------------------------------"
echo "JGODAS_PREP processing"
echo "------------------------------------------------"
echo "History: "
echo " MAY 2003 - First implementation of this new script."
echo " MAY 2004 - D.Stokes - Correct input time_mean file."
echo "                     - Clean up search for flux and SST fields. "
echo "                     - Add err_chks to exit if bathy or tesac missing."
echo " DEC 2005 - D.Stokes - Added daily averaging for bathy and tesac data."
echo "                     - Moved qc of buoy data (editprf) to run before daily averaging."
echo " JUL 2006 - D.Stokes - Processing one 29-day buoy file rather than 29 daily buoy files."
echo "                     - Removed calls to warnmsg.sh and godas_errwarn.sh because" 
echo "                       use of prodllsubmit (required to send the email warnings)" 
echo "                       is no longer encouraged.  Added warning posts to jlogfile."
echo " FEB 2007 - D.Stokes - Add processing of altimeter data."
echo " SEP 2007 - D.Stokes - Widen temporal window surface boundary conditions "
echo "                       for extended GODAS analysis.  (variable ndysbc)."
echo "                     - Processing of profile data modified to allow for deeper assimilation"
echo " JUL 2009 - C. Magee - Remove obsolete references to godas_dw."
echo " FEB 2013 - D.Behringer - Modified to run on WCOSS."
echo " MAY 2019 - H-C Lee  - Modified to run on WCOSS phase3."
#####################################################################

set -x

cd $DATA

msg="Begin GODAS PREP PROCESSING for $job on `hostname`"
postmsg "$jlogfile" "$msg"

# define executables and ush scripts

export EXECwgrib=${EXECwgrib:-$WGRIB}
export EXECgodas=${EXECgodas:-$HOMEgodas/exec}
export EXECtmSrtPrf=${EXECtmSrtPrf:-$EXECgodas/godas_tmSrtPrf}
export EXECeditPrf=${EXECeditPrf:-$EXECgodas/godas_editPrf}
export EXECavePrfDly=${EXECavePrfDly:-$EXECgodas/godas_avePrfDly}
export EXECmrgPrf=${EXECmrgPrf:-$EXECgodas/godas_mrgPrf}
export EXECmkLSAchv=${EXECmkLSAchv:-$EXECgodas/godas_mkLSAchv}
export EXECmkAsmPrf=${EXECmkAsmPrf:-$EXECgodas/godas_mkAsmPrf}
export EXECmkAsmPrfs=${EXECmkAsmPrfs:-$EXECgodas/godas_mkAsmPrfs}
export EXECmkAsmAlt=${EXECmkAsmAlt:-$EXECgodas/godas_mkAsmAlt}
export EXECgdas2g3i=${EXECgdas2g3i:-$EXECgodas/godas_gdas2g3i}
export EXECsst2g3i=${EXECsst2g3i:-$EXECgodas/godas_sst2g3i}

# Get date info from dump job or use default if not available
#   Orig version of dump script had output going to model_date_info. 
#   Migrate towards model_date_info.dump (because this step adding more vars)

if [ -s $COMIN/model_date_info.dump ]; then
  cp $COMIN/model_date_info.dump .
elif
   [ -s $COMIN/model_date_info ]; then
  cp $COMIN/model_date_info ./model_date_info.dump
fi
if [ -s model_date_info.dump ]; then
  read model_end_day ndays_bathy < model_date_info.dump
  err=$?
  if [ $err -ne 0 ]; then
     msg="Warning: Model info not input.  Set defaults"
     postmsg "$jlogfile" "$msg"
     unset model_end_day ndays_bathy
  fi
fi

#if ndays_bathy not set, use default of 29
echo ${ndays_bathy:=29}

export lag=${lag:-15}
export runlen=${runlen:-1}
export numruns=${numruns:-1}
icadj=`expr $runlen \* $numruns - 1`

# (in this run model_start_day should equal model_end_day)
echo ${model_end_day:=`sh finddate.sh $PDY d-$lag`}
if [ $icadj -ne 0 ]; then
  echo ${model_start_day:=`sh finddate.sh $model_end_day d-$icadj`}
else
  echo ${model_start_day:=$model_end_day}
fi

echo $model_start_day $model_end_day $ndays_bathy $runlen $numruns > model_date_info.prep

# This script was originally written to prepare data for the standard godas only.
# The script was later modified to process SBC's for the extended run also, but
# the above model date information pertains to the standard run.  

if test "$SENDCOM" = 'YES'
then
    [ -s model_date_info.prep ] && cp model_date_info.prep $COMOUT
fi

dsub=`expr $ndays_bathy - 1`
export StartDY=`sh finddate.sh $PDYm1 d-$dsub` 
export EndDY=$PDYm1

# set the time stamp of the time_mean field used to qc the temperature files
yrM=`echo "$model_start_day" | cut -c -4`
moM=`echo "$model_start_day" | cut -c 5-6`
dyM=`echo "$model_start_day" | cut -c 7-8`
export dteMSfx=00$yrM.$moM.$dyM

set +x
echo " "
echo "############################################################"
echo "  Copy dumped bathy, tesac and buoy data from ${COMIN}. "
echo "  If missing, print warning msg and exit with error."
echo "############################################################"
echo " "
set -x

if [ -s $COMIN/bathy.noqc ]; then
  cp $COMIN/bathy.noqc ./bathy.noqc
else
  export msg="BATHY DATA NOT FOUND"
  postmsg "$jlogfile" "$msg"
  export err=21; err_chk
fi

if [ -s $COMIN/tesac.noqc ]; then
  cp $COMIN/tesac.noqc ./tesac.noqc
else
  export msg="TESAC DATA NOT FOUND"
  postmsg "$jlogfile" "$msg"
  export err=22; err_chk
fi

if [ -s $COMIN/buoy.noqc ]; then
  cp $COMIN/buoy.noqc ./buoy.noqc
else
  export msg="BUOY DATA NOT FOUND"
  postmsg "$jlogfile" "$msg"
  export err=23; err_chk
fi

set +x
echo " "
echo " "
echo "############################################################"
echo "QC temperature profile data"
echo "############################################################"
echo " "
echo " "
echo "############################################################"
echo "     Time sort bathy"
echo "############################################################"
echo " "
set -x

export pgm=godas_tmSrtPrf
. prep_step

ln -sf bathy.noqc fort.11
ln -sf bathy.srt fort.51

startmsg
$EXECtmSrtPrf  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Edit bathy"
echo "############################################################"
echo " "
set -x

echo "Copy time_mean.$dteMSfx.nc file from previous run"
# used in two calls to godas_editPrf. netcdf, so no link.

if [ -s $COMBASE.$PDYm1/time_mean.$dteMSfx.nc ]; then
  cp $COMBASE.$PDYm1/time_mean.$dteMSfx.nc ./time_mean.nc
else
  msg='CANNOT FIND TIME_MEAN FILE.  EXIT.'
  postmsg "$jlogfile" "$msg"
  export err=25; err_chk
fi

export pgm=godas_editPrf
. prep_step

# no link used for time_mean.nc

ln -sf bathy.srt fort.11
ln -sf bathy.edt fort.51

startmsg

$EXECeditPrf  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Daily average bathy"
echo "############################################################"
echo " "
set -x

export pgm=godas_avePrfDly
. prep_step

ln -sf bathy.edt fort.11
ln -sf bathy.dav fort.51

startmsg

$EXECavePrfDly  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Time sort tesac"
echo "############################################################"
echo " "
set -x

export pgm=godas_tmSrtPrf
. prep_step

ln -sf tesac.noqc fort.11
ln -sf tesac.srt fort.51

startmsg
$EXECtmSrtPrf  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Edit tesac"
echo "############################################################"
echo " "
set -x

export pgm=godas_editPrf
. prep_step

# no link for time_mean.nc

ln -sf tesac.srt fort.11
ln -sf tesac.edt fort.51

startmsg

$EXECeditPrf  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Daily average tesac"
echo "############################################################"
echo " "
set -x

export pgm=godas_avePrfDly
. prep_step

ln -sf tesac.edt fort.11
ln -sf tesac.dav fort.51

startmsg

$EXECavePrfDly  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Time sort buoy"
echo "############################################################"
echo " "
set -x

export pgm=godas_tmSrtPrf
. prep_step

ln -sf buoy.noqc fort.11
ln -sf buoy.srt fort.51

startmsg
$EXECtmSrtPrf  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Edit buoy"
echo "############################################################"
echo " "
set -x

export pgm=godas_editPrf
. prep_step

# no link for time_mean.nc

ln -sf buoy.srt fort.11
ln -sf buoy.edt fort.51

startmsg

$EXECeditPrf  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Daily average buoy"
echo "############################################################"
echo " "
set -x

export pgm=godas_avePrfDly
. prep_step

ln -sf buoy.edt fort.11
ln -sf buoy.dav fort.51

startmsg

$EXECavePrfDly  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "Merge temperature profile data and write to assimilation file"
echo "############################################################"
echo " "
echo " "
echo "############################################################"
echo "     Merge bathy and tesac"
echo "############################################################"
echo " "
set -x

export pgm=godas_mrgPrf
. prep_step

ln -sf bathy.dav fort.11
ln -sf tesac.dav fort.12
ln -sf bat_tes.edt fort.51

startmsg

$EXECmrgPrf  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Merge bathy/tesac and buoy"
echo "############################################################"
echo " "
set -x

export pgm=godas_mrgPrf
. prep_step

ln -sf bat_tes.edt fort.11
ln -sf buoy.dav fort.12
ln -sf tmpa.edt fort.51

startmsg

$EXECmrgPrf  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Write temperature profiles to assimilation file"
echo "############################################################"
echo " "
set -x

export pgm=godas_mkAsmPrf
. prep_step

startmsg

ln -sf tmpa.edt fort.11
ln -sf $FIXgodas/godas_tmask.gbl fort.12
ln -sf tmpa.mom fort.51

$EXECmkAsmPrf $model_start_day  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "Salinity Profiles"
echo "############################################################"
echo " "
echo " "
echo "############################################################"
echo "     Make synthetic salinity profiles"
echo "############################################################"
echo " "
set -x

export pgm=godas_mkLSAchv
. prep_step

ln -sf tmpa.edt fort.11
ln -sf $FIXgodas/godas_ann.temp fort.12
ln -sf $FIXgodas/godas_ann.salt fort.13
ln -sf sala.syn fort.51

startmsg

$EXECmkLSAchv  >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "############################################################"
echo "     Write salinity profiles to assimilation file"
echo "############################################################"
echo " "
set -x

export pgm=godas_mkAsmPrfs
. prep_step

ln -sf sala.syn fort.11
ln -sf $FIXgodas/godas_tmask.gbl fort.12
ln -sf sala.mom fort.51

startmsg

$EXECmkAsmPrfs $model_start_day >> $pgmout 2> errfile
export err=$?;err_chk

set +x
echo " "
echo "#####################################################################"
echo "Send temperature and salinity assimilation files to COM for model run"
echo "#####################################################################"
echo " "
set -x

if test "$SENDCOM" = 'YES'
then
  for CTYPE in tmpa sala 
  do
    [ -s $CTYPE.mom ] && cp $CTYPE.mom $COMOUT   
  done
fi

# DBG ALT
if [ $DOALT = 'YES' ]; then
# DBG ALT
#  Altimeter data  #
CDATE=$StartDY
CTYPE=njsnal
until [ $CDATE -gt $EndDY ]; do
  DATAFILE=$COMIN/$CTYPE.$CDATE
  if test -f ${DATAFILE}; then
    cp $DATAFILE   bufr.$CTYPE
    export pgm=godas_getBufrJsn
    . prep_step
    echo Altim available for $CDATE.$CTYPE
    export XLFUNIT_11="bufr.${CTYPE}"
    export XLFUNIT_51="$CDATE.${CTYPE}.qc2048"
    $EXECgodas/godas_getBufrJsn $CDATE   >> $pgmout 2> errfile
    err=$?
    if [ $err -ne 0 ]; then
      msg="NON_CRITICAL ERROR $err IN $pgm PROCESSING FILE bufr.${CTYPE}.  CONTINUE WITH REMAINING FILES"
      postmsg "$jlogfile" "$msg"
      err=0
    fi
  fi
  CDATE=`sh finddate.sh $CDATE d+1`
done


echo Convert daily files to cycle files.
ls *.${CTYPE}.qc2048 > DList
#
export pgm=godas_dlyAlt2cyc
. prep_step
export XLFUNIT_11=DList
$EXECgodas/godas_dlyAlt2cyc   >> $pgmout 2> errfile
export err=$?; err_chk   


echo Average along track and calibrate.
export pgm=godas_trkAveCyc
. prep_step
export XLFUNIT_21="$FIXgodas/godas_calAltBufr.mom"
export XLFUNIT_61=JsnTbl
set -A file_list `ls -1 jb_c*.txt`
n=0
while [ n -lt ${#file_list[*]} ]; do
  echo "${file_list[$n]} to ca${file_list[$n]}"
  export XLFUNIT_11="${file_list[$n]}"
  export XLFUNIT_51="ca${file_list[$n]}"
  $EXECgodas/godas_trkAveCyc   >> $pgmout 2> errfile
  err=$? 
  if [ $err -ne 0 ]; then
    msg="NON_CRITICAL ERROR $err IN $pgm PROCESSING FILE ${file_list[$n]}.  CONTINUE WITH REMAINING FILES"
    postmsg "$jlogfile" "$msg"
    err=0
  fi
  ((n=$n+1))
done

echo ###  Create weekly assimilation files and consolidate 
echo ###   them into a 5-week file.
export pgm=$EXECmkAsmAlt
. prep_step
export XLFUNIT_11=JsnTbl
export XLFUNIT_21="$FIXgodas/godas_tmask.gbl"
export XLFUNIT_51=ssha.mom
$EXECmkAsmAlt $model_start_day   >> $pgmout 2> errfile
export err=$? 

  if [ $err -ne 0 ]; then
    msg="NON_CRITICAL ERROR $err IN $pgm.  CONTINUE WITHOUT SSH DATA"
    postmsg "$jlogfile" "$msg"
    err=0
  else
    if test "$SENDCOM" = 'YES'
    then
      cp ssha.mom $COMOUT/ssha.mom
    fi
  fi

#  End Altimeter processing  #
# DBG ALT
fi
# DBG ALT

# hardwire these values for now...  This job is generating SBC files to be 
#  used by both standard godas and 1-day delay godas.  
#  The latter is run for the 13 days after today's standard model run.
#  So, ndysbc=14.

ndysbc=${ndysbc:-14}
nfields=`expr $ndysbc + 2`
  
set +x
echo " "
echo "############################################################"
echo "FLUXES"
echo "############################################################"
echo " "
echo "############################################################"
echo "     Retrieve flux fields for $nfields days surrounding run dates "
echo "     from $COMIN, where they have been placed by the dump job."
echo "     Look for surrounding days if any of the fields are missing."
echo "############################################################"
echo " "

set -x

flux_start_date=`sh finddate.sh $model_start_day d-1`
flux_end_date=`sh finddate.sh $model_start_day d+$ndysbc`

cat <<cntrlEOF > namelist.cntrlatm
        &cntrlatm  ndysbc=$ndysbc,
        /
cntrlEOF

for flux in taux tauy senflx latflx prate netlwflx netswflx
do

  cp $COMIN/${flux}.* .
  cat ${flux}.* > $flux

  case $flux in
    senflx)   f11=senflx;   f12=/dev/null; fbase=fsh;  inp=sh ;;
    latflx)   f11=latflx;   f12=/dev/null; fbase=flh;  inp=lh ;;
    netlwflx) f11=netlwflx; f12=/dev/null; fbase=ful;  inp=lw ;;
    netswflx) f11=netswflx; f12=/dev/null; fbase=fsr;  inp=sw ;;
    prate)    f11=latflx;   f12=prate;     fbase=fslt; inp=sf ;;
    taux)     f11=taux;     f12=/dev/null; fbase=taux; inp=tx ;;
    tauy)     f11=tauy;     f12=/dev/null; fbase=tauy; inp=ty ;;
    *) echo ERROR IN FLUX CASE.;  export err=15; err_chk;;
  esac


  set +x
  echo " "
  echo "############################################################"
  echo "     Make flux file ${fbase}.mom "
  echo "############################################################"
  echo " "
  set -x

  export pgm=godas_gdas2g3i
  . prep_step

  ln -sf "$f11" fort.11
  ln -sf "$f12" fort.12

  ln -sf "$FIXgodas/godas_mask.gdas" fort.13
  if test "$f11" = 'taux' -o "$f11" = 'tauy'
  then
    ln -sf "$FIXgodas/godas_umask.gbl" fort.14
  else
    ln -sf "$FIXgodas/godas_tmask.gbl" fort.14
  fi
  ln -sf "namelist.cntrlatm" fort.31
  ln -sf "$fbase.mom" fort.51

  startmsg

  $EXECgdas2g3i $model_start_day "$inp" >> $pgmout 2> errfile
  export err=$?;err_chk


  if test "$SENDCOM" = 'YES'
  then
    cp $fbase.mom $COMOUT/$fbase.mom
  fi

done

set +x
echo " "
echo "############################################################"
echo "SST"
echo "############################################################"
echo " "
echo "############################################################"
echo "     Retrieve sst fields for $nfields days surrounding run dates "
echo "     from $COMIN, where they have been placed by the dump job."
echo "############################################################"
set -x

cp $COMIN/sst.grb.* .

declare -a gribfile
gribfile=(`ls sst.grb.*`)
echo ${gribfile[0]}
$EXECwgrib -s ${gribfile[0]} | egrep TMP:sfc | $EXECwgrib -i -grib ${gribfile[0]} -o sst.grb

nfield=1
while [ $nfield -lt $nfields ]
do
  echo ${gribfile[$nfield]}
  $EXECwgrib -s ${gribfile[0]} | egrep TMP:sfc | $EXECwgrib -i -grib ${gribfile[$nfield]} -o sst.grb -append
  nfield=`expr $nfield + 1`
done

$EXECwgrib sst.grb -4yr -s -d all -ieee -o sst.i3e > sst.info
cat sst.info

export pgm=godas_sst2g3i
. prep_step

  ln -sf sst.i3e fort.11
  ln -sf sst.info fort.12
  ln -sf $FIXgodas/godas_tmask.gbl fort.14
  ln -sf namelist.cntrlatm fort.31
  ln -sf sst.mom fort.51

startmsg

$EXECsst2g3i $model_start_day >> $pgmout 2> errfile
export err=$?;err_chk

if test "$SENDCOM" = 'YES'
then
  cp sst.mom $COMOUT/sst.mom
fi

#####################################################################
# GOOD RUN
set +x
echo "************** $job COMPLETED NORMALLY"
echo "************** $job COMPLETED NORMALLY"
echo "************** $job COMPLETED NORMALLY"
set -x
#####################################################################

msg="$job HAS COMPLETED NORMALLY!"
postmsg "$jlogfile" "$msg"

############## END OF SCRIPT #######################
