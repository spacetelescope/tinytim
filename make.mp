OBJ3 = makemaps.o misc.o image.o

.c.o: 
	$(TINYTIMCC) $(CFLAGS) -c $<

makemaps: $(OBJ3)
	$(TINYTIMCC) $(OBJ3) $(LFLAGS) -o makemaps
	./makemaps

image.o:	image.c tinytim.h

makemaps.o:	makemaps.c tinytim.h

misc.o:		misc.c tinytim.h
