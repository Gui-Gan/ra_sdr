ra_sdr
======

Programm to use a DVB-T USB Dongle as an SDR for radioastronomical test.
It reads data for Tint time and compute the average FFT(1024) and stores in file.
to compile: gcc testrtl.c -o testrtl -lrtlsdr -lfftw3 -lm -Wall
based on testrtl.c from http://www.m0dts.co.uk/files/simple_rtlsdr_fft.c
and RAFFT2.c from http://y1pwe.co.uk/RAProgs/index.html

Example: 
./ra_sdr -f 40000000 -s 2000000 -g 10 -i 12 -v 1 -d 0 data_out.txt

Not Fully tested....
