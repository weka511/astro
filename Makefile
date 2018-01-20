# Makefile snarfed from https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile

CPPFLAGS=-g -O3 -pthread -I/sw/include/root  -std=gnu++11
LDFLAGS=-g -O3
LDLIBS=
CC=gcc
CXX=g++
RM=rm -f
SRCS=barnes_hut.cpp tree.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
TESTS=tests.exe 
TEST_CASES=test-tree.cpp
TEST_OBJS=$(subst .cpp,.o,$(TEST_CASES))
MAIN=galaxy.exe 
TARGETS=$(MAIN) $(TESTS) 

all : $(TARGETS)

clean :
	${RM} *.o

tests : $(TESTS)
	./$(TESTS)
	
depend: .depend

.depend: $(SRCS) $(TARGETS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;
	
$(MAIN): $(OBJS) galaxy.o
	${CXX} $(LDFLAGS) -o $(MAIN) galaxy.o ${OBJS} ${LDLIBS}
	
$(TESTS): $(OBJS) tests.o $(TEST_OBJS)
	${CXX} $(LDFLAGS) -o $(TESTS) tests.o $(TEST_OBJS) ${OBJS} ${LDLIBS}

distclean: clean
	$(RM) *~ .depend

include .depend
