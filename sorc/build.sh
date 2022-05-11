#!/bin/bash
#
source ../versions/build.ver
module reset
module use `pwd`
module load build_godas.lua
module list

base=`pwd`
echo $base
exec=$base/../exec

declare -a dirList

dirLst=(bufr_argodump bufr_prepmods avePrfDly bathybfr buoybtbfr dailyflx Dly2MnthNc Dly2PntdNc dlyAlt2cyc dw_mkAsmAlt dw_mkAsmPrf dw_mkAsmPrfs dw_mkEvNc editPrf gdas2g3i getBufrJsn mkAsmAlt mkAsmPrf mkAsmPrfs mkEvNc mkLSAchv mkSEv mppnccombine mrgPrf nc2grib sst2g3i tmSrtPrf trkAveCyc)
#dirLst=(bufr_argodump)

let n=0
while [ $n -lt ${#dirLst[*]} ]; do
  echo ${dirLst[$n]}
  export short_name=${dirLst[$n]}
  if [ ${short_name} == 'bufr_argodump' ] || [ ${short_name} == 'bufr_prepmods' ]
  then
    export long_name=${short_name}
  else
    export long_name=godas_${short_name}
  fi
  cd ${base}/${long_name}.fd
  make all
  make clean
  mv ${long_name} $exec
  let n=$n+1
done

#-- for godas_dw_mom3das.fd
 cd ${base}/godas_dw_mom3das.fd/compile
 ./cmp_iGM3d.sh

#-- for godas_mom3das.fd
 cd ${base}/godas_mom3das.fd/compile
 ./cmp_iGM3d.sh

