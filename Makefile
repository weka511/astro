# Makefile snarfed from https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile

CPPFLAGS=-g -O3 -pthread -I/sw/include/root  -std=gnu++11
LDFLAGS=-g -O3
LDLIBS=
CC=gcc
CXX=g++
RM=rm -f
SRCS=barnes_hut.cpp galaxy.cpp tree.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all : galaxy.exe

clean :
	${RM} ${OBJS}

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;
	
galaxy.exe: $(OBJS)
	${CXX} $(LDFLAGS) -o galaxy.exe ${OBJS} ${LDLIBS}

distclean: clean
	$(RM) *~ .depend

include .depend
