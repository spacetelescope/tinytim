LFLAGS_NORMAL = "LFLAGS = -lm"

CFLAGS_SOLARIS = "CFLAGS = -fsingle -O"

CFLAGS_LINUX = "CFLAGS = -O3 -ffast-math -ffloat-store"
CFLAGS_MACOSX = "CFLAGS = -O3 -ffast-math -ffloat-store -DMACOSX"

CFLAGS_INTEL = "CFLAGS = -O3"
LFLAGS_INTEL = "LFLAGS = -lsvml -limf -lm"

CFLAGS_OTHER = "CFLAGS = -f -O"
CFLAGS_DEC = "CFLAGS = -f -O"
CFLAGS_HP = "CFLAGS = -Ae -O"

CFLAGS_THREADED_SOLARIS = "CFLAGS = -mt -fsingle -O -DTT_THREADED -D_REENTRANT -D_POSIX_C_SOURCE=199506L"
LFLAGS_THREADED_SOLARIS = "LFLAGS = -lm -lthread -lpthread" 

CFLAGS_THREADED_LINUX = "CFLAGS = -O3 -ffast-math -ffloat-store -D_REENTRANT -DTT_THREADED"
LFLAGS_THREADED_LINUX = "LFLAGS = -lm -pthread" 

CFLAGS_THREADED_MACOSX = "CFLAGS = -O3 -ffast-math -ffloat-store -DMACOSX -D_REENTRANT -DTT_THREADED"
LFLAGS_THREADED_MACOSX = "LFLAGS = -lm" 

CFLAGS_THREADED_LINUX_ICC = "CFLAGS = -O3 -parallel -D_REENTRANT -DTT_THREADED"
LFLAGS_THREADED_LINUX_ICC = "LFLAGS = -lsvml -limf -lm -pthread" 

INTEL_CC = "TINYTIMCC = icc"
DEFAULT_CC = "TINYTIMCC = cc"
 
nothing: 
	@echo
	@echo "To compile Tiny Tim, use one of the following as appropriate :"
	@echo
	@echo "          make solaris"
	@echo "          make linux"
	@echo "          make macosx"
	@echo "          make hp"
	@echo
	@echo "For multiprocessor systems, use :"
	@echo
	@echo "          make threadedsolaris"
	@echo "          make threadedlinux"
	@echo "          make threadedmacosx"
	@echo
	@echo "If you are compiling on something else, then you must edit the"
        @echo "Makefile according to the manual and enter:"
	@echo
	@echo "          make other"
	@echo

solaris:
	@make -f make.t1 $(CFLAGS_SOLARIS) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_SOLARIS) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_SOLARIS) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_SOLARIS) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@echo

hp:
	@make -f make.t1 $(CFLAGS_HP) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_HP) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_HP) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_HP) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@echo

decstation: 	
	@make -f make.t1 $(CFLAGS_DEC) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_DEC) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_DEC) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_DEC) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@echo

linux: 
	@make -f make.t1 $(CFLAGS_LINUX) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_LINUX) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_LINUX) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_LINUX) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@echo 

macosx: 
	@make -f make.t1 $(CFLAGS_MACOSX) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_MACOSX) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_MACOSX) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_MACOSX) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@echo
 
threadedsolaris: 
	@make -f make.t1 $(CFLAGS_THREADED_SOLARIS) $(LFLAGS_THREADED_SOLARIS) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_THREADED_SOLARIS) $(LFLAGS_THREADED_SOLARIS) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_THREADED_SOLARIS) $(LFLAGS_THREADED_SOLARIS) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_THREADED_SOLARIS) $(LFLAGS_THREADED_SOLARIS) $(DEFAULT_CC)
	@echo

threadedlinux: 
	@make -f make.t1 $(CFLAGS_THREADED_LINUX) $(LFLAGS_THREADED_LINUX) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_THREADED_LINUX) $(LFLAGS_THREADED_LINUX) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_THREADED_LINUX) $(LFLAGS_THREADED_LINUX) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_THREADED_LINUX) $(LFLAGS_THREADED_LINUX) $(DEFAULT_CC)
	@echo

threadedmacosx: 
	@make -f make.t1 $(CFLAGS_THREADED_MACOSX) $(LFLAGS_THREADED_MACOSX) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_THREADED_MACOSX) $(LFLAGS_THREADED_MACOSX) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_THREADED_MACOSX) $(LFLAGS_THREADED_MACOSX) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_THREADED_MACOSX) $(LFLAGS_THREADED_MACOSX) $(DEFAULT_CC)
	@echo

threadedlinuxintel:
	@make -f make.t1 $(CFLAGS_THREADED_LINUX_ICC) $(LFLAGS_THREADED_LINUX_ICC) $(INTEL_CC)
	@make -f make.t2 $(CFLAGS_THREADED_LINUX_ICC) $(LFLAGS_THREADED_LINUX_ICC) $(INTEL_CC)
	@make -f make.t3 $(CFLAGS_THREADED_LINUX_ICC) $(LFLAGS_THREADED_LINUX_ICC) $(INTEL_CC)
	@make -f make.mp $(CFLAGS_THREADED_LINUX_ICC) $(LFLAGS_THREADED_LINUX_ICC) $(INTEL_CC)
	@echo

intel: 
	@make -f make.t1 $(CFLAGS_INTEL) $(LFLAGS_INTEL) $(INTEL_CC)
	@make -f make.t2 $(CFLAGS_INTEL) $(LFLAGS_INTEL) $(INTEL_CC)
	@make -f make.t3 $(CFLAGS_INTEL) $(LFLAGS_INTEL) $(INTEL_CC)
	@make -f make.mp $(CFLAGS_INTEL) $(LFLAGS_INTEL) $(INTEL_CC)
	@echo

other:
	@make -f make.t1 $(CFLAGS_OTHER) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t2 $(CFLAGS_OTHER) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.t3 $(CFLAGS_OTHER) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@make -f make.mp $(CFLAGS_OTHER) $(LFLAGS_NORMAL) $(DEFAULT_CC)
	@echo

clean:
	rm -f *.o tiny1 tiny2 tiny3 makemaps *.map *.new

