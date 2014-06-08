CC = gcc
CFLAGS = -lrtlsdr -lfftw3 -lm -Wall
CFLAGS += -DVERBOSE_MODE
LIBS=  

all: ra_sdr

ra_sdr: ra_sdr.o 
	gcc $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o
	rm -f ra_sdr
