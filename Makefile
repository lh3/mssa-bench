CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2
CXXFLAGS=	$(CFLAGS)
DFLAGS=	
OBJS=		ksa.o sais.o divsufsort.o qsufsort.o ssort.o dc3.o is.o
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

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cc)

# DO NOT DELETE

mssac.o: ketopt.h kseq.h
