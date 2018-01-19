# Makefile snarfed from https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile

CPPFLAGS=-g -O3 -pthread -I/sw/include/root  -std=gnu++11
LDFLAGS=-g -O3
LDLIBS=
CC=gcc
CXX=g++
RM=rm -f
SRCS=barnes_hut.cpp  tree.cpp 
OBJS=$(subst .cpp,.o,$(SRCS))

all : galaxy.exe tests.exe

clean :
	${RM} ${OBJS}

depend: .depend

.depend: $(SRCS) galaxy.cpp tests.cpp
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;
	
galaxy.exe: $(OBJS) galaxy.o
	${CXX} $(LDFLAGS) -o galaxy.exe galaxy.o ${OBJS} ${LDLIBS}
	
tests.exe: $(OBJS) tests.o
	${CXX} $(LDFLAGS) -o tests.exe tests.o ${OBJS} ${LDLIBS}

distclean: clean
	$(RM) *~ .depend

include .depend
