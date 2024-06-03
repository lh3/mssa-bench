CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O3
CXXFLAGS=	$(CFLAGS)
DFLAGS=	
OBJS=		ksa.o ksa64.o libsais.o libsais64.o sais.o divsufsort.o qsufsort.o ssort.o dc3.o is.o
PROG=		mssa-bench
INCLUDES=
LIBS=		-lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

mssa-bench:$(OBJS) mssac.o
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) mssac.o -o $@ $(LIBS)

ksa.o:ksa.c
		$(CC) -c $(CFLAGS) -o $@ $<

ksa64.o:ksa.c
		$(CC) -c $(CFLAGS) -D_KSA64 -o $@ $<

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cc)

# DO NOT DELETE

libsais.o: libsais.h
libsais64.o: libsais.h libsais64.h
mssac.o: libsais64.h ketopt.h kseq.h
