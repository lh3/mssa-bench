CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
OBJS=		ksa.o ksa64.o libsais.o libsais64.o
EXE=		mssa-bench
INCLUDES=
LIBS=		-lz

ifneq ($(openmp),)
	CPPFLAGS=-DLIBSAIS_OPENMP
	CFLAGS+=-fopenmp
endif

.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@
.cc.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

mssa-bench:$(OBJS) mssac.o
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(OBJS) mssac.o -o $@ $(LIBS)

ksa.o:ksa.c
	$(CC) -c $(CFLAGS) -o $@ $<

ksa64.o:ksa.c
	$(CC) -c $(CFLAGS) -D_KSA64 -o $@ $<

clean:
	rm -fr gmon.out *.o ext/*.o a.out $(EXE) *~ *.a *.dSYM session*

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c *.cc)

# DO NOT DELETE

libsais.o: libsais.h
libsais64.o: libsais.h libsais64.h
mssac.o: libsais.h libsais64.h ketopt.h kseq.h
