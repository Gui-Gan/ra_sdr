CC = gcc
CFLAGS = -Wall
CFLAGS += -DVERBOSE_MODE
LIBS= -lrtlsdr -lfftw3 -lm

all: ra_sdr

ra_sdr: ra_sdr.o 
	gcc $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o
	rm -f ra_sdr
