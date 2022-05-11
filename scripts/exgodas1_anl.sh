#!/bin/sh
#########################################################################
set +xa
echo "-----------------------------------------------------------------"
echo "JGODAS_ANL - This script runs the MOM_3 ocean model for the GODAS"
echo "-----------------------------------------------------------------"
echo "History: "
echo " 08/12/2003 - Initial implementation of this script on the IBM SP (frost/snow) "
echo "  "
echo " 05/26/2004 - D.Stokes - Fix input background field data and timestamp"
echo "                         of output time_mean and restart files "
echo " 02/14/2007 - D.Stokes - Add input files for processing JASON SSHA data."
echo "                       - Remove calls to godas_errwarn.sh which no longer sends warnings"
echo "                         to programmer via email.  Post warnings to jlogfile instead."
echo " 09/14/2007 - D.Stokes - Make this script more generic to handle standard or multi-day run"
echo " 02/17/2013 - D.Behringer - Adapted this script for WCOSS"
echo " 06/13/2019 - H-C Lee  - Adapted this script for WCOSS Phase3"
set -xa

###########################################################
# Processing/Flow of script
# 1) Copy this job's input files to the working directory
# 2) For each day (number determined by variable numruns): 
# 2a)    prepare background covariances
# 2b)    run mom3
# 3) Save final output files to $COMOUT
# 4) Remove files located in this job's working directory
###########################################################

#############################################################
# Ensure that this job is pointed to the production temporary
# working directory, which is where all of this job's input 
# and output files will be saved during execution
#############################################################
cd $DATA

msg="$job HAS BEGUN on `hostname`"
postmsg "$jlogfile" "$msg"

# define executables

export EXECgodas=${EXECgodas:-$HOMEgodas/exec}
export EXECmkEvNc=${EXECmkEvNc:-$EXECgodas/godas_mkEvNc}
export EXECmkSEv=${EXECmkSEv:-$EXECgodas/godas_mkSEv}
export EXECmppnccombine=${EXECmppnccombine:-$EXECgodas/godas_mppnccombine}

# default is standard 1-day run with 15-day lag
stdlag=${stdlag:-15}       # delay relative to PDY of standard GODAS run
lag=${lag:-${stdlag:-15}}    # delay relative to PDY to end this series of runs 


# for info purposes only, get number of days of bathy data from prep job or use default if not available
#
if [ -s $COMIN/model_date_info.prep ]; then
  cp $COMIN/model_date_info.prep .
  read dummy dummy ndays_bathy dummy dummy < model_date_info.prep
  err=$?
  if [ $err -ne 0 ]; then
     msg="Warning: Model info not input.  Set defaults"
     postmsg "$jlogfile" "$msg"
     unset ndays_bathy
  fi
fi

ndays_bathy=${ndays_bathy:-29}

runlen=${runlen:-1}          # length of model run
numruns=${numruns:-1}        # number of times to cycle runs of length runlen

icadj=`expr $runlen \* $numruns - 1`  # initial condition date adjustment

echo ${model_end_day:=`sh finddate.sh $PDY d-$lag`}
if [ $icadj -ne 0 ]; then
  echo ${model_start_day:=`sh finddate.sh $model_end_day d-$icadj`}
else
  echo ${model_start_day:=$model_end_day}
fi

echo $model_start_day $model_end_day $ndays_bathy $runlen $numruns > model_date_info.anl

if test "$SENDCOM" = 'YES'
then
  cp model_date_info.anl $COMOUT
fi

yrM=`echo "$model_start_day" | cut -c -4`
moM=`echo "$model_start_day" | cut -c 5-6`
dyM=`echo "$model_start_day" | cut -c 7-8`

export dteMSfx=00$yrM.$moM.$dyM

set +x
echo "Begining GODAS for $model_end_day"
echo "IC Start Date is $model_start_day"
echo "Do $numruns Run(s) of Run_length $runlen Day(s)"
set -x

model_end_day_p1=`sh finddate.sh $model_end_day d+1`

yr=`echo "$model_end_day_p1" | cut -c -4`
mo=`echo "$model_end_day_p1" | cut -c 5-6`
dy=`echo "$model_end_day_p1" | cut -c 7-8`

dteSfx=00$yr.$mo.$dy
set +x
echo "Time Stamp of final output files will be $dteSfx"

#
# Copy this job's input files to the working directory
#

echo "Copy fixed runtime data from $FIXgodas"
set -x
cp $FIXgodas/godas_ls.dat ./ls.dat      
cp $FIXgodas/godas_sss.mom ./sss.mom
cp $FIXgodas/godas_kmtE40.dta ./kmt.dta
cp $FIXgodas/godas_htE40.dta ./ht.dta
cp $FIXgodas/godas_sshc.mom ./sshc.mom
cp $FIXgodas/godas_cdnz.mom ./cdnz.mom

set +x
echo "Copy runtime data from godas_prep job"
set -x
cp $COMIN/tmpa.mom ./tmpa.mom
cp $COMIN/sala.mom ./sala.mom
cp $COMIN/taux.mom ./taux.mom
cp $COMIN/tauy.mom ./tauy.mom
cp $COMIN/fsh.mom ./fsh.mom
cp $COMIN/flh.mom ./flh.mom
cp $COMIN/fslt.mom ./fslt.mom
cp $COMIN/ful.mom ./ful.mom
cp $COMIN/fsr.mom ./fsr.mom
cp $COMIN/sst.mom ./sst.mom
cp $COMIN/ssha.mom ./ssha.mom

set +x
echo "Copy namelist files from $PARMgodas"
set -x
cp $PARMgodas/godas_namelist.asmcntrl ./namelist.asmcntrl
cp $PARMgodas/godas_namelist.blmix ./namelist.blmix
cp $PARMgodas/godas_namelist.contrl ./namelist.contrl
cp $PARMgodas/godas_namelist.diagn ./namelist.diagn
cp $PARMgodas/godas_namelist.ictime ./namelist.ictime
cp $PARMgodas/godas_namelist.io ./namelist.io
cp $PARMgodas/godas_namelist.isopyc ./namelist.isopyc
cp $PARMgodas/godas_namelist.kppmix ./namelist.kppmix
cp $PARMgodas/godas_namelist.mbcin ./namelist.mbcin
cp $PARMgodas/godas_namelist.mixing ./namelist.mixing
cp $PARMgodas/godas_namelist.ppmix ./namelist.ppmix
cp $PARMgodas/godas_namelist.riglid ./namelist.riglid
cp $PARMgodas/godas_namelist.smagnl ./namelist.smagnl
cp $PARMgodas/godas_namelist.svd ./namelist.svd
cp $PARMgodas/godas_namelist.tsteps ./namelist.tsteps

#  This job reads SBC files to be used by both standard godas 
#  and extended godas.  If the the extended analysis is run for 
#  the 13 days after today's standard model run, ndysbc=14.

ndysbc=${ndysbc:-14}  

cat <<cntrlEOF > namelist.cntrlatm
        &cntrlatm  ndysbc=$ndysbc,
        /
cntrlEOF


# Get initial conditions 
#   For standard run, this is likely from PDYm1.  For extended run, this is likley from PDY.
set +x
echo "Copy initial restart.$dteMSfx.dta and time_mean.$dteMSfx.nc from $COMIN_IC"
set -x
if [ -s $COMIN_IC/restart.$dteMSfx.dta ]; then
  cp $COMIN_IC/restart.$dteMSfx.dta ./restart.dta
else
  msg="CANNOT FIND RESTART FILE.  EXIT."
  postmsg "$jlogfile" "$msg"
  export err=21;err_chk
fi
if [ -s $COMIN_IC/time_mean.$dteMSfx.nc ] ; then
  cp $COMIN_IC/time_mean.$dteMSfx.nc ./time_mean.nc
else
  msg="CANNOT FIND TIME_MEAN FILE.  EXIT."
  postmsg "$jlogfile" "$msg"
  export err=22;err_chk
fi



CDATE=$model_start_day
i_loop=1
until [ $CDATE -gt $model_end_day ]
do

set +x
  echo "Prepare background error variance"
set -x

# default is kass=30 for shallow assimilation, set kass=35 for deep.

  kass=${kass:-30}

  export pgm=godas_mkEvNc
  . prep_step
  $EXECmkEvNc -f time_mean.nc -k $kass -o tvv.mom >> $pgmout 2> errfile
  export err=$?;err_chk

  export pgm=godas_mkSEv
  . prep_step
  $EXECmkSEv -f tvv.mom -o svv.mom >> $pgmout 2> errfile
  export err=$?;err_chk


# ======================================================================
# Run the model
# ======================================================================

  echo $CDATE > run_date   

  set +x
  echo "==>Running..."
  set -x

  export pgm=godas_mom3das
  . prep_step
  export MP_SHARED_MEMORY=no

  export MP_EAGER_LIMIT=65536
  export MP_EUIDEVELOP=min
  export MP_EUIDEVICE=sn_all
  export MP_EUILIB=us
# export MP_MPILIB=mpich2

  export MP_USE_BULK_XFER=yes
  export MPICH_ALLTOALL_THROTTLE=0
  export MP_COLLECTIVE_OFFLOAD=no

  export KMP_STACKSIZE=1024m
  export MP_TASK_AFFINITY=core:1
  export OMP_NUM_THREADS=1

#  export NPES=32

# mpiexec -n 32 -ppn 32 $EXECmom3das >> $pgmout 2> errfile
# mpiexec -n 10 -ppn 10 $EXECmom3das >> $pgmout 2> errfile
  mpiexec -n $ncpus -ppn $ncpus $EXECmom3das >> $pgmout 2> errfile
  export err=$?;err_chk

  CDATE_p1=`sh finddate.sh $CDATE d+1`
  yr=`echo "$CDATE_p1" | cut -c -4`
  mo=`echo "$CDATE_p1" | cut -c 5-6`
  dy=`echo "$CDATE_p1" | cut -c 7-8`
  dteSfx=00$yr.$mo.$dy

  set +x
  echo " time_mean"
  set -x
  declare -a file_list
  file_list=(`ls -1 time_mean.$dteSfx.dta.nc.????`)
  if [ ${#file_list[*]} -gt 0 ]; then
    tm_comb=time_mean.$dteSfx.nc
    if [ ${#file_list[*]} -gt 1 ]; then
      $EXECmppnccombine $tm_comb ${file_list[*]}
    else
      mv ${file_list[0]} $tm_comb
    fi
    [ $CDATE -ne $model_end_day ] && cp $tm_comb ./time_mean.nc   # for next loop, if necessary
  fi
  [ $CDATE -ne $model_end_day ] && cp restart.$dteSfx.dta ./restart.dta   # for next loop, if necessary

## Output the 7-day lag files:
  if [ $i_loop -eq 7 ] 
  then
     dteSfx7=$dteSfx
     tm_comb7=$tm_comb 
  fi

  CDATE=$CDATE_p1
  i_loop=`expr $i_loop + 1`

done

if test "$SENDCOM" = 'YES'
then
  [ -s tvv.mom ] && cp tvv.mom $COMOUT
  [ -s svv.mom ] && cp svv.mom $COMOUT
  [ -s $tm_comb ] && cp $tm_comb $COMOUT
  [ -s "$tm_comb7" ] && cp $tm_comb7 $COMOUT
  [ -s restart.$dteSfx.dta ] && cp restart.$dteSfx.dta $COMOUT
  [ -s restart.$dteSfx7.dta ] && cp restart.$dteSfx7.dta $COMOUT
fi

if test "$SENDDBN" = 'YES'
then
  if [ -s "$tm_comb" ]; then
    $DBNROOT/bin/dbn_alert MODEL GODAS_NC $job $COMOUT/$tm_comb
  fi
  if [ -s "$tm_comb7" ]; then
    $DBNROOT/bin/dbn_alert MODEL GODAS_NC $job $COMOUT/$tm_comb7
  fi
fi

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

################## END OF SCRIPT #######################
