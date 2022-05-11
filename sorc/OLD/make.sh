#!/bin/bash
#
module purge
#module load envvar/1.0
#module load PrgEnv-intel/8.1.0
#module load craype/2.7.8
#module load intel/19.1.3.304
#module load cray-mpich/8.1.7
#module load w3emc/2.7.3
#module load prod_envir/2.0.4
#module load w3nco/2.4.1
#module load bufr/11.4.0
#module load bacio/2.4.1
#module load netcdf/3.6.3

. ../versions/build.ver
module load envvar/${envvar_ver}
module load PrgEnv-intel/${PrgEnv_intel_ver}
module load craype/${craype_ver}
module load intel/${intel_ver}
module load cray-mpich/${cray_mpich_ver}
module load w3emc/${w3emc_ver}
module load prod_envir/${prod_envir_ver}
module load w3nco/${w3nco_ver}
module load bufr/${bufr_ver}
module load bacio/${bacio_ver}
module load netcdf/${netcdf_ver}
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

