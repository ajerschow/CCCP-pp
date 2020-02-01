objects = CCCP++.o parse.o 
#CFLAGS =  -g -w -DVERBOSE
#CFLAGS =  -w -g
CFLAGS =  -w 
linkflags = -lm -lstdc++
CC = gcc

all: CCCP++

parse.cc: create_parse_cc Makefile
	perl create_parse_cc > parse_pp
	flex++ -oparse.cc parse_pp

CCCP++: Makefile $(objects) $(headers) parse.cc
	$(CC) $(CFLAGS) $(objects) -o $@ $(libflags) $(linkflags)

clean:
	rm *.o

bkup:
	tar -cvzf ../CCCP++bkup.tgz README Makefile *.h *.cc create_parse_cc *.par

dist:
	tar -cvzf ../CCCP++.tgz Makefile gcd.h parse.cc CCCP++.cc ParsedVars.h create_parse_cc mqmas?.par pfgbpste1.par

test: CCCP++
	./CCCP++ < mqmas1.par
	./CCCP++ < mqmas2.par
	./CCCP++ < mqmas3.par
	./CCCP++ < pfgbpste1.par

.cc.o:
	$(CC) -c $(CFLAGS) $*.cc $(libflags)

