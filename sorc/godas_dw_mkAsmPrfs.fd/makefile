#!/bin/sh
#
short_name=mkAsmPrfs
long_name=godas_dw_mkAsmPrfs

#export LIBS= /gpfs/dell1/nco/ops/nwprod/lib/w3emc/v2.3.0/ips/18.0.1/impi/18.0.1/libw3emc_v2.3.0_4.a /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a /gpfs/dell1/nco/ops/nwprod/lib/bufr/v11.2.0/ips/18.0.1/libbufr_v11.2.0_4_64.a /gpfs/dell1/nco/ops/nwprod/lib/bacio/v2.0.2/ips/18.0.1/libbacio_v2.0.2_4.a
export LIBS = ${W3EMC_LIB4} ${W3NCO_LIB4} ${BUFR_LIB4} ${BACIO_LIB4}
#export LIBS = ${W3EMC_LIB4} ${W3NCO_LIB4} ${BACIO_LIB4} /apps/ops/prod/libs/intel/19.1.3.304/bufr/11.4.0/lib/libbufr_4.a
export FC=ifort
FFLAGS=-O3 -free -check bounds -convert big_endian -traceback

all: $(long_name)
tmUtil_mod.o: tmUtil_mod.f
	$(FC) $(FFLAGS) -c tmUtil_mod.f

$(short_name).o: $(short_name).f
	$(FC) $(FFLAGS) -c $(short_name).f

OBJ = tmUtil_mod.o $(short_name).o

$(long_name): $(OBJ)
	$(FC) -o $(long_name) $(OBJ) $(LIBS)

clean:
	rm *.o *.mod
