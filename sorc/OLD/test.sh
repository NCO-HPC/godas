#! /bin/bash

. ../build.ver

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

