CC=g++
CFLAGS=-Wall -std=c++11 -O3

all: uniqueLoci_build

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ -fopenmp

uniqueLoci_build: uniqueLoci_build.o arg_parse.o CM.o
	$(CC) $(CFLAGS) uniqueLoci_build.o arg_parse.o CM.o -o uniqueLoci_build -fopenmp

clean:
	rm *.o