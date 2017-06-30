all : barnes_hut.exe
barnes_hut.exe: barnes_hut.cpp
	g++ -O3 -g -o barnes_hut barnes_hut.cpp
