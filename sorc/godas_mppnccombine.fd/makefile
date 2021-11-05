SHELL=	/bin/sh
#
#
short_name = mppnccombine
long_name = godas_$(short_name)

CC = icc
#NCDF = -L/usrx/local/prod/packages/ips/18.0.1/netcdf/3.6.3/lib -lnetcdf
#INC = -I/usrx/local/prod/packages/ips/18.0.1/netcdf/3.6.3/include
NCDF = -L${NETCDF_LIB} -lnetcdf
INC = -I${NETCDF_INC}
CFLAGS = -O3

all: $(long_name)
$(short_name).o: $(short_name).c
	$(CC) $(CFLAGS) $(INC) -c $(short_name).c

OBJ = $(short_name).o

$(long_name): $(OBJ)
	$(CC) -o $(long_name) $(OBJ) $(NCDF)

clean:
	rm *.o

