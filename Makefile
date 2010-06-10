# Path to FFTW includes
FFTWINCDIR = /usr/local/include
# Path to FFTW libraries
FFTWLIBDIR = /usr/local/lib
# How to link with the FFTW libs
FFTWLINK = -L$(FFTWLIBDIR) -lfftw3f

# Other include directory (for CFITSIO, libsla, which is in PRESTO)
OTHERINCLUDE = 
# Other link directory (for CFITSIO, libsla, which is in PRESTO)
OTHERLINK = -L/usr/local/lib -lcfitsio -L$(PRESTO)/lib -lsla 

# Source directory
SRCDIR = $(shell pwd)

# git commit-hash
GITHASH = $(shell git rev-parse HEAD)

# Which C compiler
CC = gcc
CFLAGS = -I$(FFTWINCDIR) -DSRCDIR=\"$(SRCDIR)\"\
	-DGITHASH=\"$(GITHASH)\"\
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64\
	-O3 -Wall -W
#	-g -Wall -W
CLINKFLAGS = $(CFLAGS)

# When modifying the CLIG files, the is the location of the clig binary
CLIG = /usr/local/bin/clig
# Rules for CLIG generated files
%_cmd.c : %_cmd.cli
	$(CLIG) -o $*_cmd -d $<

OBJS = chkio.o vectors.o sla.o transpose.o wapp.o rescale.o\
	wapp_head_parse.o wapp_y.tab.o write_psrfits.o

wapp2psrfits: wapp2psrfits.o wapp2psrfits_cmd.o $(OBJS)
	$(CC) $(CLINKFLAGS) -o $@ wapp2psrfits.o wapp2psrfits_cmd.o $(OBJS) $(FFTWLINK) $(OTHERLINK) -lm

# Default indentation is K&R style with no-tabs,
# an indentation level of 4, and a line-length of 85
indent:
	indent -kr -nut -i4 -l85 *.c
	rm *.c~

clean:
	rm -f *.o *~ *#

cleaner: clean
	rm -f wapp2psrfits
