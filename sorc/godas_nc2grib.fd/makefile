SHELL=	/bin/sh
#
#
short_name = nc2grib
long_name = godas_$(short_name)

FC = ifort
#LIBS = -L/nwprod/lib -lw3emc_4 -lw3nco_4 -lbufr_4_64 -lbacio_4
#LIBS = /gpfs/dell1/nco/ops/nwprod/lib/w3emc/v2.3.0/ips/18.0.1/impi/18.0.1/libw3emc_v2.3.0_4.a /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a /gpfs/dell1/nco/ops/nwprod/lib/bufr/v11.2.0/ips/18.0.1/libbufr_v11.2.0_4_64.a /gpfs/dell1/nco/ops/nwprod/lib/bacio/v2.0.2/ips/18.0.1/libbacio_v2.0.2_4.a
#NCDF = -L/usrx/local/prod/packages/ips/18.0.1/netcdf/3.6.3/lib -lnetcdf
#INC = -I/usrx/local/prod/packages/ips/18.0.1/netcdf/3.6.3/include
LIBS = ${W3EMC_LIB4} ${W3NCO_LIB4} ${BUFR_LIB4} ${BACIO_LIB4} 
#LIBS = ${W3EMC_LIB4} ${W3NCO_LIB4} ${BACIO_LIB4} /apps/ops/prod/libs/intel/19.1.3.304/bufr/11.4.0/lib/libbufr_4.a
NCDF = -L${NETCDF_LIB} -lnetcdf
INC = -I${NETCDF_INC}
FFLAGS = -O3 -fixed -convert big_endian -traceback

all: $(long_name)
$(short_name).o: $(short_name).f
	$(FC) $(FFLAGS) $(INC) -c $(short_name).f

OBJ = $(short_name).o

$(long_name): $(OBJ)
	$(FC) -o $(long_name) $(OBJ) $(NCDF) $(LIBS)

clean:
	rm *.o

