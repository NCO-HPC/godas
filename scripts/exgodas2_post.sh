#!/bin/sh
############################ EXGODAS_POST.SH.SMS ##############################
#                                                                             #
# This script performs GODAS post processing                                  #
#                                                                             #
#                                                              February 2008  #
#                                                                             #
###############################################################################

#####################################################################
set +xa
echo "------------------------------------------------"
echo "JGODAS_POST processing:"
echo "Generate pentad and/or monthly fields from GODAS data."
echo "Either start sum or add to current running sum for pentad and monthly fields."
echo "Compute pentad or monthly average when full time period is available".
echo "------------------------------------------------"
echo "History: September 2005 - First development run of new script."
echo "             March 2008 - Implementation into production."
echo "             February 2012 - Adapted to WCOSS."
set -xa
#####################################################################

cd $DATA

msg="Begin GODAS POST PROCESSING for $job on `hostname` for PDY=$PDY"
postmsg "$jlogfile" "$msg"

# define executables
export EXECgodas=${EXECgodas:-$HOMEgodas/exec}
export EXECDly2MnthNc=${EXECDly2MnthNc:-${EXECgodas}/godas_Dly2MnthNc}
export EXECDly2PntdNc=${EXECDly2PntdNc:-${EXECgodas}/godas_Dly2PntdNc}
export EXECnc2grib=${EXECnc2grib:-${EXECgodas}/godas_nc2grib}

set +x
echo " "
echo "############################################################"
echo " get date info from anl job or use default if not available"
echo "############################################################"
set -x

if [ -s $COMIN/model_date_info.anl ]; then
  cp $COMIN/model_date_info.anl .
  read model_start_day model_end_day dum < model_date_info.anl
  export err=$?
  if [ $err -ne 0 ]; then
     msg="Warning: Model info not input.  Set defaults"
     postmsg "$jlogfile" "$msg"
     unset model_start_date model_end_day numruns
  fi
fi


lag=${lag:-15}
echo ${model_end_day:=`sh finddate.sh $PDY d-$lag`}
model_end_day_p1=`sh finddate.sh $model_end_day d+1`

yr=`echo "$model_end_day_p1" | cut -c -4`
mo=`echo "$model_end_day_p1" | cut -c 5-6`
dy=`echo "$model_end_day_p1" | cut -c 7-8`

dteSfx=00$yr.$mo.$dy
set +x
echo "Time Stamp of input time_mean file is $dteSfx"


# Copy this job's input files to the working directory

set +x
echo "Copy running sum file(s) from PREVIOUS run, if available"
set -x
[ -s $COMBASE.$PDYm1/SUM_time_mean.P.dta ] && cp $COMBASE.$PDYm1/SUM_time_mean.P.dta .
[ -s $COMBASE.$PDYm1/SUM_time_mean.M.dta ] && cp $COMBASE.$PDYm1/SUM_time_mean.M.dta .

set +x
echo "Copy time_mean.$dteSfx.nc from TODAY'S run"
set -x
if [ -s $COMIN/time_mean.$dteSfx.nc ]; then
  cp $COMIN/time_mean.$dteSfx.nc .
else
  export err=22
  msg="ERROR:  CANNOT FIND TIME_MEAN FILE.  EXIT $err"
  postmsg "$jlogfile" "$msg"
  err_chk
fi

do_month=${do_month:-1}

if [ $do_month -eq 1 ]; then
# MONTHLY

set +x
echo ""
echo "###############################################################################"
echo " Accumulate GODAS fields for monthly avg.  Makes avg if all days available"
echo "###############################################################################"
echo ""
set -x

export pgm=godas_Dly2MnthNc
. prep_step
$EXECDly2MnthNc $model_end_day_p1 >> $pgmout 2> errfile
export err_month=$?

if [ $err_month -eq  0 ]; then

  if test "$SENDCOM" = 'YES'
  then
    [ -s SUM_time_mean.M.dta ] && cp SUM_time_mean.M.dta $COMOUT
    [ -s time_mean.$dteSfx.M.nc ] && cp time_mean.$dteSfx.M.nc $COMOUT
  fi

  if [ -s time_mean.$dteSfx.M.nc ]; then
    set +x
    echo "###########################################################################"
    echo "Month complete.  Extract desired fields from pentad and convert to grib"
    echo "###########################################################################"
    set -x
    export pgm=godas_nc2grib
    . prep_step
    $EXECnc2grib time_mean.$dteSfx.M.nc >> $pgmout 2> errfile
    export err=$?
    err_chk

    # grib file gets stamped date of last day of period 
    yrG=`echo "$model_end_day" | cut -c -4`
    moG=`echo "$model_end_day" | cut -c 5-6`
    dyG=`echo "$model_end_day" | cut -c 7-8`
    dteGSfx=00$yrG.$moG.$dyG

    if test "$SENDCOM" = 'YES'
    then
      if [ -s time_mean.$dteGSfx.M.grb ]; then
        monthly_name=godas.M.${yrG}${moG}.grb
        cp time_mean.$dteGSfx.M.grb $COMOUT/$monthly_name
      fi
    fi
  fi
else
  msg="**WARNING:  ERROR $err_month IN MONTHLY STEP.  SCRIPT CONTINUING....   "
  postmsg "$jlogfile" "$msg"
fi

fi    # end if do_month


# PENTAD

set +x
echo " "
echo "###############################################################################"
echo " Accumulate GODAS fields for pentad.  Makes pentad if all days available"
echo "###############################################################################"
echo " "
set -x

export pgm=godas_Dly2PntdNc
. prep_step
$EXECDly2PntdNc $model_end_day_p1 >> $pgmout 2> errfile
export err_pentad=$?

if [ $err_pentad -eq  0 ]; then
  if test "$SENDCOM" = 'YES'; then
    [ -s SUM_time_mean.P.dta ] && cp SUM_time_mean.P.dta $COMOUT
    [ -s time_mean.$dteSfx.P.nc ] && cp time_mean.$dteSfx.P.nc $COMOUT
  fi
  if [ -s time_mean.$dteSfx.P.nc ]; then
    set +x
    echo " "
    echo "############################################################################"
    echo "Pentad complete.  Extract desired fields from pentad and convert to grib"
    echo "############################################################################"
    echo " "
    set -x
    export pgm=godas_nc2grib
    . prep_step
    $EXECnc2grib time_mean.$dteSfx.P.nc >> $pgmout 2> errfile
    export err=$?
    err_chk
  fi

# grib file gets stamped date of last day of period
  yrG=`echo "$model_end_day" | cut -c -4`
  moG=`echo "$model_end_day" | cut -c 5-6`
  dyG=`echo "$model_end_day" | cut -c 7-8`
  dteGSfx=00$yrG.$moG.$dyG

  if test "$SENDCOM" = 'YES'; then
    if [ -s time_mean.$dteGSfx.P.grb ]; then
      pentad_name=godas.P.${yrG}${moG}${dyG}.grb
      cp time_mean.$dteGSfx.P.grb $COMOUT/$pentad_name
    fi
  fi

else
  msg="**WARNING:   ERROR $err_pentad IN PENTAD STEP.  SCRIPT CONTINUING....   "
  postmsg "$jlogfile" "$msg"
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

################## END OF SCRIPT #######################
