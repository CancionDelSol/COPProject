CC=gcc
CFLAG= -Wall -I. -O3

TARGETS=fft

all: $(TARGETS)

# Run all
runAll: all
	clear
	./main

# Original unoptimized version
fft: main.o microtime.o
	$(CC) -o $@ $^

fft.o: main.cpp microtime.h
	$(CC) $(CFLAG) -c $<

microtime.o: microtime.c microtime.h
	$(CC) $(CFLAG) -c $<

clean:
	rm -f *.o *~ core $(TARGETS)