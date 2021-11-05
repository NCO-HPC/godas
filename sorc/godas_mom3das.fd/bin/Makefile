# Makefile created by mkmf.mom3 $Id: Makefile 100938 2019-06-17 17:50:05Z hyun-chul.lee@noaa.gov $ 


.DEFAULT:
	-touch $@
all: a.out
SRC =
OBJ =
clean: neat
	-rm -f .cppdefs $(OBJ) a.out
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
a.out: $(OBJ)
	$(LD) $(OBJ) -o a.out $(LDFLAGS)
