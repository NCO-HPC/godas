# ###########################################################################
# makefile (bufr_argodump) -  23 May 2018 
# 
#                use to convert argo type 031005 to tesac format 031002
# ----
#
# Source files & objects
#
SRCS	= bufr_argodump.f  
OBJS	= bufr_argodump.o


FC =            ifort
#LIBS=           -L/nwprod/lib -lw3nco_4 -lw3emc_4 -lbufr_v10.2.3_4_64
#LIBS=           -L/gpfs/dell1/nco/ops/nwprod /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a /gpfs/dell1/nco/ops/nwprod/lib/bufr/v11.2.0/ips/18.0.1/libbufr_v11.2.0_4_64.a
#LIBS=   -L${NWROOT} ${W3NCO_LIB4} ${BUFR_LIB4} ${W3EMC_LIB4} 
#LIBS=   -L${W3NCO_LIB4} ${BUFR_LIB4} ${W3EMC_LIB4} 
LIBS=   ${W3NCO_LIB4} ${BUFR_LIB4} ${W3EMC_LIB4} 
#LIBS=   -L${NWROOT} ${W3NCO_LIB4} ${W3EMC_LIB4} /apps/ops/prod/libs/intel/19.1.3.304/bufr/11.4.0/lib/libbufr_4.a
#CHKOUT = 	-ftrapuv  -check all  -fp-stack-check  -fstack-protector
DEBUG = 	-g -traceback
FFLAGS =        -O2 $(DEBUG) $(CHKOUT) -free
LDFLAGS=

CMD = 		bufr_argodump

all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(CMD) $(SRCS) $(LIBS)

clean:
	-rm -f $(OBJS)

clobber:	clean
	-rm -f $(CMD)    

install:
	mv $(CMD) ../../exec/

# -- end of Makefile

