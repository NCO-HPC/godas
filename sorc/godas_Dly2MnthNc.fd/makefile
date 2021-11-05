SHELL=	/bin/sh
#
#
short_name = Dly2MnthNc
long_name = godas_$(short_name)

FC = ifort
LIBS = ${W3EMC_LIB4} ${W3NCO_LIB4}
#NCDF = -L${NETCDF_LIB} -lnetcdf -lnetcdff
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

