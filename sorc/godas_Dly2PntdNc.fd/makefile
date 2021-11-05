SHELL=	/bin/sh
#
#
short_name = Dly2PntdNc
long_name = godas_$(short_name)

FC = ifort
#LIBS = /gpfs/dell1/nco/ops/nwprod/lib/w3emc/v2.3.0/ips/18.0.1/impi/18.0.1/libw3emc_v2.3.0_4.a /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a
LIBS = ${W3EMC_LIB4} ${W3NCO_LIB4}
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

