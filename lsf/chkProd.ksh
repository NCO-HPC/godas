#!/bin/ksh
#

utilscript=/nwprod/util/ush
hr=00
export cycle=t${hr}z
$utilscript/setpdy.sh > scrtch
. PDY
# echo $PDY

tdy=$PDY
adate=`sh $utilscript/finddate.sh $PDY d-14`
edate=`sh $utilscript/finddate.sh $PDY d-1`

# echo $edate

pdir=/com/godas/prod

yr=`echo ${adate} | cut -c 1-4`
mo=`echo ${adate} | cut -c 5-6`
dy=`echo ${adate} | cut -c 7-8`
gtm=${pdir}/godas.${tdy}/time_mean.00${yr}.${mo}.${dy}.nc
g2tm=${pdir}/godas2.${tdy}/time_mean.00${yr}.${mo}.${dy}.nc

echo " "
echo "godas.${tdy} ${yr}-${mo}-${dy}"
chkTmNc -f $gtm

echo " "
echo "godas2.${tdy} ${yr}-${mo}-${dy}"
chkTmNc -f $g2tm

yr=`echo ${edate} | cut -c 1-4`
mo=`echo ${edate} | cut -c 5-6`
dy=`echo ${edate} | cut -c 7-8`
gextm=${pdir}/godas_ext.${tdy}/time_mean.00${yr}.${mo}.${dy}.nc

echo " "
echo "godas_ext.${tdy} ${yr}-${mo}-${dy}"
chkTmNc -f $gextm
echo " "

rm scrtch NMCDATE ncepdate PDY

