SHELL=	/bin/sh
#
#
SRCS=	buoybtbfr.f

OBJS=	buoybtbfr.o

# Tunable parameters
#
# FC        Name of the Fortran compiling system to use
# LDFLAGS   Options of the loader
# FFLAGS    Options of the compiler
# DEBUG     Options of the compiler included for debugging
# LIBS      List of libraries
# CMD       Name of the executable
#
FC =        ifort
LDFLAGS = 
#LIBS =     /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a /gpfs/dell1/nco/ops/nwprod/lib/bufr/v11.2.0/ips/18.0.1/libbufr_v11.2.0_4_64.a
LIBS = ${W3NCO_LIB4} ${BUFR_LIB4}
#LIBS = ${W3NCO_LIB4} /apps/ops/prod/libs/intel/19.1.3.304/bufr/11.4.0/lib/libbufr_4.a

CMD =      godas_buoybtbfr
DEBUG =    
#DEBUG = -ftrapuv  -check all  -fp-stack-check  -fstack-protector
FFLAGS =   -O3 -convert big_endian -fp-model strict  -auto -g -traceback $(DEBUG)


# Lines from here down should not need to be changes.  They are
# the actual rules which make uses to build CMD.
#
all:        $(CMD)

$(CMD):      $(OBJS)
	$(FC)  $(LDFLAGS) -o $(@) $(OBJS)  $(LIBS)

clean:
	-rm -f $(OBJS)

clobber:        clean
	-rm -f $(CMD)

void:   clobber
	-rm -f $(SRCS) makefile

install:   
	mv $(CMD) ../../exec/

