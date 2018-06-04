CC=g++
CFLAGS=-Wall -std=c++11 -O3

all: uniqueLoci

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $^ -fopenmp

uniqueLoci: uniqueLoci.o arg_parse.o CM.o
	$(CC) $(CFLAGS) uniqueLoci.o arg_parse.o CM.o -o uniqueLoci -fopenmp

clean:
	rm *.o