PREFIX := $(DESTDIR)$(PREFIX)
ifndef $(PREFIX)
	PREFIX = $(DESTDIR)/usr/local
endif

# system variables
SRCDIR = $(
BINDIR = $(PREFIX)/bin
DATADIR = $(PREFIX)/share/tinytim
DOCDIR = $(PREFIX)/share/doc/tinytim
INCLUDEDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib

# project variables
HEADERS = ../include
OPTIMIZE_LEVEL = -O3
CFLAGS = -I. -I$(HEADERS) -I$(INCLUDEDIR) $(OPTIMIZE_LEVEL) -Wall -Wextra
LDFLAGS = -L$(LIBDIR) -lm
USE_FASTMATH = 0
USE_THREADING = 0

# handle user-defined build flags
ifneq ($(USE_FASTMATH), 0)
	OPTIMIZE_LEVEL = -O2
	CFLAGS := $(CFLAGS) -ffast-math -ffloat-store
	LDFLAGS := $(LDFLAGS) -lc
endif

ifneq ($(USE_THREADING), 0)
	CFLAGS := $(CFLAGS) -D_REENTRANT -DTT_THREADED
	LDFLAGS := -lpthread $(LDFLAGS)
endif

data = $(wildcard ../data/*)
docs = $(wildcard ../docs/*)

makemaps_OBJS = image.o \
		makemaps.o \
		misc.o
tiny1_OBJS = compgrid.o \
	     entry.o \
	     inputs.o \
	     misc.o \
	     outputs.o \
	     rdpupil.o \
	     spectrum.o \
	     system.o \
	     tiny1.o \
	     v2v3.o \
	     zernike.o
tiny2_OBJS = acs.c \
	     convotf.o \
	     cpupil.o \
	     distort.o \
	     entry.o \
	     fitsio.o \
	     fft.o \
	     halo.o \
	     image.o \
	     inputs.o \
	     inter.o \
	     intpsf.o \
	     map.o \
	     misc.o \
	     monopsf.o \
	     opd.o \
	     outputs.o \
	     param.o \
	     polypsf.o \
	     psf.o \
	     pupil.o \
	     rdpupil.o \
	     rotate.o \
	     system.o \
	     tiny2.o \
	     v2v3.o \
	     vignet.o \
	     zernike.o
tiny3_OBJS = acs.o \
	     distort.o \
	     entry.o \
	     fft.o \
	     fitsio.o \
	     halo.o \
	     image.o \
	     inputs.o \
	     misc.o \
	     param.o \
	     rdpupil.o \
	     system.o \
	     tiny3.o \
	     v2v3.o \
	     zernike.o

.c.o:
	$(CC) $(CFLAGS) $(LDFLAGS) -c $<

all: config.h makemaps tiny1 tiny2 tiny3

config.h: config.h.in
	sed "s|@CONFIG_DEFAULT_DIR@|${DATADIR}|g" $< > $@

makemaps: $(makemaps_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(makemaps_OBJS)

tiny1: $(tiny1_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(tiny1_OBJS)

tiny2: $(tiny2_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(tiny2_OBJS)

tiny3: $(tiny3_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(tiny3_OBJS)

.PHONY:
install-dirs:
	mkdir -p $(BINDIR)
	mkdir -p $(DATADIR)
	mkdir -p $(DOCDIR)

.PHONY:
install-bin: makemaps tiny1 tiny2 tiny3
	install -D -m 755 $^ $(BINDIR)

.PHONY:
install-data: $(data)
	install -D -m 644 $^ $(DATADIR)
	cd $(DATADIR) && $(BINDIR)/makemaps

.PHONY:
install-docs: $(docs)
	install -D -m 644 $^ $(DOCDIR)

install: all install-dirs install-docs install-bin install-data

clean:
	rm -f *.o \
		tiny1 \
		tiny2 \
		tiny3 \
		makemaps \
		config.h
