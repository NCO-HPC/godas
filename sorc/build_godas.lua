--#%Module############################################################
--#
--# GODAS
--#    - Global Ocean DAS
--####################################################################
--proc ModulesHelp { } {
--  puts stderr "Sets environment variables for GODAS"
--  puts stderr "This module initializes the environment"
--  puts stderr "to build the GODAS software at NCEP"
--}

whatis("GODAS module for compilation")

-- Load Intel Compiler

load("PrgEnv-intel/"..os.getenv("PrgEnv_intel_ver"))
load("craype/"..os.getenv("craype_ver"))
load("intel/"..os.getenv("intel_ver"))
load("cray-mpich/"..os.getenv("cray_mpich_ver"))
load("w3emc/"..os.getenv("w3emc_ver"))
load("w3nco/"..os.getenv("w3nco_ver"))
load("bufr/"..os.getenv("bufr_ver"))
load("bacio/"..os.getenv("bacio_ver"))
load("netcdf/"..os.getenv("netcdf_ver"))

