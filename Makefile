all : barnes_hut.exe
barnes_hut.exe: barnes_hut.cpp tree.h
	g++ -O3 -g -o barnes_hut barnes_hut.cpp
