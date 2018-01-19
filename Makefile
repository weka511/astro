CPPFLAGS=-g -pthread -I/sw/include/root 
LDFLAGS=-g
LDLIBS=-L/sw/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
       -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lz -L/sw/lib -lfreetype -lz \
       -Wl,-framework,CoreServices -Wl,-framework,ApplicationServices -pthread -Wl,-rpath,/sw/lib/root \
       -lm -ldl

all : barnes_hut.exe

SRCS=barnes_hut.cpp tree.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

barnes_hut.exe: $(OBJS)
	g++ $(LDFLAGS) -o barnes_hut.exe ${OBJS} 

barnes_hut.o: barnes_hut.cpp tree.h
	g++ $(CPPFLAGS) -c barnes_hut.cpp

tree.o: tree.h tree.cpp
	g++ $(CPPFLAGS) -c tree.cpp
	
#barnes_hut.exe: barnes_hut.cpp tree.h
#	g++ -O3 -g -o barnes_hut barnes_hut.cpp
