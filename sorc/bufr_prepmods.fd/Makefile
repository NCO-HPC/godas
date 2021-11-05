# makefile (bufr_prepmods) -  4 Feb 2013
# 
#  4 Feb 2013 JWhiting cp'd fr IBM/AIX ver of 22 Mar 2012 (dated 11 Nov 2003)
#              config's for use on WCOSS linux (ifort cmplr, new lib specs)
#                rm'd profiler refs
#                parameterized cmpl debug options (-g, -traceback))
#                set optimization to default (-O2)
#                rm'd unneeded compiler options (list, source, etc)
#                set to use default variable sizes (-i4, -r4)
#                 added developmental cmpl options (CHKOUT)
# ----
#
# Source files & objects
SRCS	= bufr_prepmods.f clonlat.f mnemonics.f obname.f raddate.f ssmipn.f tideg.f
OBJS	= $(SRCS:.f=.o) 


FC =            ifort
#LIBS=           -L/nwprod/lib -lw3nco_4 -lw3emc_4 -lbufr_v10.2.3_4_64   
#LIBS= -L/gpfs/dell1/nco/ops/nwprod /gpfs/dell1/nco/ops/nwprod/lib/w3nco/v2.0.6/ips/18.0.1/libw3nco_v2.0.6_4.a /gpfs/dell1/nco/ops/nwprod/lib/bufr/v11.2.0/ips/18.0.1/libbufr_v11.2.0_4_64.a
#LIBS= -L${NWROOT} ${W3NCO_LIB4} ${BUFR_LIB4}
#LIBS= -L${W3NCO_LIB4} ${BUFR_LIB4}
LIBS= ${W3NCO_LIB4} ${BUFR_LIB4}
#LIBS= -L${NWROOT} ${W3NCO_LIB4} /apps/ops/prod/libs/intel/19.1.3.304/bufr/11.4.0/lib/libbufr_4.a
#LIBS=           -L/nwprod/lib -lw3nco_4 -lw3emc_4 -lbufr_v10.2.5_s_64
#CHKOUT = 	-ftrapuv  -check all  -fp-stack-check  -fstack-protector
DEBUG = 	-g -traceback
FFLAGS =        -O2 $(DEBUG) $(CHKOUT)
LDFLAGS=

CMD = 		bufr_prepmods


all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)

clean:
	-rm -f $(OBJS)

clobber:	clean
	-rm -f bufr_prepmods

install:
	mv $(CMD) ../../exec/

# -- end of Makefile

