# Makefile snarfed from https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile

CPPFLAGS=-g -O3 -pthread -I/sw/include/root  -std=gnu++11
LDFLAGS=-g -O3
LDLIBS=
CC=gcc
CXX=g++
RM=rm -f
SRCS=barnes_hut.cpp  tree.cpp 
OBJS=$(subst .cpp,.o,$(SRCS))
TESTS=tests.exe
MAIN=galaxy.exe 
TARGETS=$(MAIN) $(TESTS) 

all : $(TARGETS)

clean :
	${RM} *.o

tests : $(TESTS)
	./$(TESTS)
	
depend: .depend

.depend: $(SRCS) galaxy.cpp tests.cpp
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;
	
$(MAIN): $(OBJS) galaxy.o
	${CXX} $(LDFLAGS) -o $(MAIN) galaxy.o ${OBJS} ${LDLIBS}
	
$(TESTS): $(OBJS) tests.o
	${CXX} $(LDFLAGS) -o $(TESTS) tests.o ${OBJS} ${LDLIBS}

distclean: clean
	$(RM) *~ .depend

include .depend
