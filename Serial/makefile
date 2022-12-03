CC=gcc
CFLAG= -Wall -I. -O3

TARGETS=fft

all: $(TARGETS)

# Run all
runAll: all
	clear
	./fft

# Original unparallelized version
fft: main.o microtime.o
	$(CC) -o $@ $^ -lm

fft.o: main.cpp microtime.h
	$(CC) $(CFLAG) -c $<

microtime.o: microtime.c microtime.h
	$(CC) $(CFLAG) -c $<

clean:
	rm -f *.o *~ core $(TARGETS)