SHELL=/bin/sh
#
SRCS= dailyflx.f

OBJS= dailyflx.o

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
#LIBS = /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a /gpfs/dell1/nco/ops/nwprod/lib/bacio/v2.0.2/ips/18.0.1/libbacio_v2.0.2_4.a
LIBS = ${W3NCO_LIB4} ${BACIO_LIB4}
CMD =      godas_dailyflx
DEBUG =    
#DEBUG = -ftrapuv  -check all  -fp-stack-check  -fstack-protector
FFLAGS =   -O3 -convert big_endian -fp-model strict  -auto -g -traceback $(DEBUG)


# Lines from here down should not need to be changes.  They are
# the actual rules which make uses to build CMD.
#
all:        $(CMD)

$(CMD):      $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	-rm -f $(OBJS)

clobber:        clean
	-rm -f $(CMD) $(CMD).prof

void:   clobber
	-rm -f $(SRCS) makefile

install:
	mv $(CMD) ../../exec

