SHELL=	/bin/sh
#
#
short_name = mkEvNc
long_name = godas_$(short_name)

FC = ifort
#LIBS = -L/nwprod/lib -lw3emc_4 -lw3nco_4
#LIBS = /gpfs/dell1/nco/ops/nwprod/lib/w3emc/v2.3.0/ips/18.0.1/impi/18.0.1/libw3emc_v2.3.0_4.a /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a
#NCDF = -L/usrx/local/prod/packages/ips/18.0.1/netcdf/3.6.3/lib -lnetcdf
#INC = -I/usrx/local/prod/packages/ips/18.0.1/netcdf/3.6.3/include
LIBS = ${W3EMC_LIB4} ${W3NCO_LIB4}
NCDF = -L${NETCDF_LIB} -lnetcdf
INC = -I${NETCDF_INC}
FFLAGS = -O3 -free -convert big_endian -traceback

all: $(long_name)
$(short_name).o: $(short_name).f
	$(FC) $(FFLAGS) $(INC) -c $(short_name).f

OBJ = $(short_name).o

$(long_name): $(OBJ)
	$(FC) -o $(long_name) $(OBJ) $(NCDF) $(LIBS)

clean:
	rm *.o

