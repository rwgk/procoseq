CC=cc
CCFLAGS=-O#-bbinder:/usr/lib/bind -bmaxdata:2000000000

DEFTYPES=#-DNORMAL_TYPES
OLIMIT=#-Olimit 1000

what:
	@echo 'make what?'

procoseq: procoseq.o fperiod.o longer.o inumber.o
	$(CC) $(CCFLAGS) -o procoseq procoseq.o fperiod.o longer.o inumber.o -lm

procoseq.o: procoseq.c
	$(CC) $(CCFLAGS) $(DEFTYPES) -c procoseq.c $(OLIMIT)

rcs: rcs.o fperiod.o longer.o inumber.o
	$(CC) $(CCFLAGS) -o rcs rcs.o fperiod.o longer.o inumber.o

rcs.o: rcs.c
	$(CC) $(CCFLAGS) $(DEFTYPES) -c rcs.c $(OLIMIT)

abc: abc.o fperiod.o longer.o inumber.o
	$(CC) $(CCFLAGS) -o abc abc.o fperiod.o longer.o inumber.o

pascal_triangle: pascal_triangle.o longer.o
	$(CC) $(CCFLAGS) -o pascal_triangle pascal_triangle.o longer.o

.c.o:
	$(CC) $(CCFLAGS) $(DEFTYPES) -c $*.c

clean:
	rm -f procoseq rcs abc pascal_triangle core
	rm -f procoseq.o rcs.o abc.o pascal_triangle.o
	rm -f fperiod.o longer.o inumber.o
