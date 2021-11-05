SHELL=	/bin/sh
#
#
SRCS=	mrgPrf.f

OBJS=	mrgPrf.o

FC =        ifort
LDFLAGS = 
#LIBS = /gpfs/dell1/nco/ops/nwprod/lib/w3emc/v2.3.0/ips/18.0.1/impi/18.0.1/libw3emc_v2.3.0_4.a /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a
LIBS = ${W3EMC_LIB4} ${W3NCO_LIB4}

CMD =      godas_mrgPrf
DEBUG =    
FFLAGS =   -O3 -free -convert big_endian -traceback

# Lines from here down should not need to be changes.  They are
# the actual rules which make uses to build CMD.
#
all:        $(CMD)

$(CMD):      $(OBJS)
	$(FC)  $(LDFLAGS) -o $(@) $(OBJS)  $(LIBS)

clean:
	-rm -f $(OBJS)

