CFLAGS = -O -DDIMACS

OBJS = netgen.o index.o random.o

.PHONY: all clean

netgen: $(OBJS)
	$(CC) -o $@ $(OBJS)

netgen.o: netgen.c netgen.h

index.o: index.c netgen.h

random.o: random.c

clean : 
	rm -f $(OBJS) netgen

