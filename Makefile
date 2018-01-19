CPPFLAGS=-g -O3 -pthread -I/sw/include/root 
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
	
galaxy.exe: $(OBJS)
	${CXX} $(LDFLAGS) -o galaxy.exe ${OBJS} ${LDLIBS}

barnes_hut.o: barnes_hut.cpp tree.h barnes_hut.h

galaxy.o: galaxy.cpp barnes_hut.h tree.h

tree.o: tree.h tree.cpp
