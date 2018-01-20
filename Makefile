# Makefile snarfed from https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile

CPPFLAGS=-g -O3 -pthread -I/sw/include/root  -std=gnu++11
LDFLAGS=-g -O3
LDLIBS=
CC=gcc
CXX=g++
RM=rm -f
SRCS=barnes_hut.cpp tree.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
TEST_MAIN=test-main.exe 
TEST_CASES=test-tree.cpp
TEST_CASES_OBJS=$(subst .cpp,.o,$(TEST_CASES))
TEST_OBJS=$(TEST_CASES_OBJS) test-main.o
MAIN=galaxy.exe 
TARGETS=$(MAIN) $(TESTS) 

all : $(TARGETS)

clean :
	${RM} *.o

tests : $(TEST_MAIN)
	./$(TEST_MAIN)
	
depend: .depend

.depend: $(SRCS) $(TEST_CASES) galaxy.cpp test-main.cpp
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;
	
$(MAIN): $(OBJS) galaxy.o
	${CXX} $(LDFLAGS) -o $(MAIN) galaxy.o ${OBJS} ${LDLIBS}
	
$(TEST_MAIN): $(OBJS) $(TEST_OBJS)
	${CXX} $(LDFLAGS) -o $(TEST_MAIN) $(TEST_OBJS) ${OBJS} ${LDLIBS}

distclean: clean
	$(RM) *~ .depend

include .depend
