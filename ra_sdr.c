/*
The MIT License (MIT)

Copyright (c) 2014 Guillermo Gancio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
// based on testrtl.c from http://www.m0dts.co.uk/files/simple_rtlsdr_fft.c
// and RAFFT2.c from http://y1pwe.co.uk/RAProgs/index.html

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>  
#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _WIN32
#include <unistd.h>
#else
#include <Windows.h>
#include <io.h>
#include <fcntl.h>
#include "getopt/getopt.h"
#endif
#include "rtl-sdr.h"

#define DEF_int_t 1 
#define DEFAULT_SAMPLE_RATE 2048000
#define DEF_Freq 30000000
#define DEF_Gain 207
#define FFTs 1024
#define DEF_debug 0
#define DEFAULT_ASYNC_BUF_NUMBER 32
#define DEFAULT_BUF_LENGTH (16 * 16384)
#define MINIMAL_BUF_LENGTH 512
#define MAXIMAL_BUF_LENGTH (256 * 16384)
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
#define PI (6.28318530717959/2.0)

float yy,rr,rrr,rav,amp;
int xx,r,pts,a,int_t=DEF_int_t,debug=DEF_debug;
int flag=0,da[FFTs*2];
long int file_end,count,p_num;
double dats[FFTs*4],datr[FFTs*4],sig_pow,nois_pow;
int retval,opt;
int devices;
static rtlsdr_dev_t *dev;
int n;
uint8_t *buf;   //unsigned 8bit int - I didn't know what it was!, the _t must be 'type'
char * fftresult;
uint32_t frequency = DEF_Freq;
int gain = DEF_Gain;
static uint32_t bytes_to_read = DEF_int_t*2;
uint32_t out_block_size = DEFAULT_BUF_LENGTH;
static int do_exit = 0;
float sample_rate_aux=DEFAULT_SAMPLE_RATE;
FILE *file;
char *filename = NULL;
FILE *file1;
char *filename1 = NULL;
void four(double[],int,int);
void sum_dat(void);
void out_dat(void);

void usage(void)
{
	fprintf(stderr,
		"ra_sdr, an I/Q recorder for RTL2832 based DVB-T receivers, for RadioAstronomy use. Ver 1.0, Not Fully Tested.\n\n"
		"Usage:\t -f frequency_to_tune_to (default 30000000 [Hz)]\n"
		"\t[-s samplerate (default: 2048000 Hz)]\n"
		"\t[-d device_index (default: 0)]\n"
		"\t[-g gain (default: 20.7)]\n"
		"\t[-i 'kind of' integration Time in seconds (default: 1)]\n"
                "\t[-v vervose (0 or 1, default: 0)]\n"
		"\tfilename \n\n");
	exit(1);
}

static void rtlsdr_callback(unsigned char *buf, uint32_t len, void *ctx)
{
	if (ctx) 
		{
		if (do_exit)
			return;
		if ((bytes_to_read > 0) && (bytes_to_read < len)) 
			{
			len = bytes_to_read;
			do_exit = 1;
			rtlsdr_cancel_async(dev);
			}
		/*take fourier transform*/
		pts=FFTs;
		p_num=(len+1)/pts/2;
		int s,ss;
		count=0;
		for(ss=0;ss<p_num;ss++)
		  {
			for(s=0;s<2*pts;s++)
				{
				dats[s]=(float)buf[s+(count*2*FFTs)]+0.0; 
				if(dats[s]>127)
					dats[s]=(dats[s]-127.5)/128.0;
				else 
					dats[s]=(dats[s]-127.5)/128.0;
				}
				count++;
			four(dats-1,FFTs,-1);
			sum_dat();
			}
		if (len!= len) 
			{
			fprintf(stderr, "Short write, samples lost, exiting!\n");
			rtlsdr_cancel_async(dev);
			}
		if (bytes_to_read > 0)
		bytes_to_read -= len;
		}
}


int main(int argc, char **argv)
	{
	int sample_aux;
	double aux;
	uint32_t dev_index = 0;
	uint32_t out_block_size = DEFAULT_BUF_LENGTH;
	time_t start, stop;
	while ((opt = getopt(argc, argv, "d:f:g:s:b:i:v:t::")) != -1) {
		switch (opt) {
		case 'd':
			dev_index = atoi(optarg);
			break;
		case 'f':
			frequency = (uint32_t)atof(optarg);
			break;
		case 'g':
			gain = (int)(atof(optarg) * 10); /* tenths of a dB */
			break;
		case 's':
			sample_rate_aux = (uint32_t)atof(optarg);
			break;
		case 'i':
			int_t=atoi(optarg);
			break;
                case 'v':
                       debug=atoi(optarg);
                        break;
		default:
			usage();
			break;
		}
	}

	if (argc <= optind) {
		usage();		
	} else {
		filename = argv[optind];
	}
	time(&start);
	sample_aux=(int)(int_t/(1.0/sample_rate_aux));
	aux=pow(2,(log10(sample_aux)/log10(2)));
	sample_aux=(int)aux;
	bytes_to_read=sample_aux*2;
	if(debug)printf("Sample_aux %d Samples %d block size %zu\n",sample_aux,bytes_to_read,out_block_size);
	file = fopen("tmp", "w");
	file1 = fopen(filename, "w");
	if (!file1)
	{
			fprintf(stderr, "Failed to open %s\n", filename1);
			return(1);
	}
	devices = rtlsdr_get_device_count();

	for(n=0;n<devices;n++)
	{
		if(debug)printf("Device %d: %s\n\n",devices,rtlsdr_get_device_name(n));
	}
	if (devices>0){
		retval = rtlsdr_open(&dev, dev_index);	
		if(debug)printf("Open= %d\n",retval);
	}else{
		if(debug)printf("No Devices found...!");
		exit(0);
	}
	buf = malloc((out_block_size) * sizeof(uint8_t));
	//configure rtlsdr settings
	retval = rtlsdr_set_sample_rate(dev, sample_rate_aux);
	if (retval < 0)
		fprintf(stderr, "WARNING: Failed to set sample rate.\n");

	retval = rtlsdr_set_center_freq(dev, frequency);
	if (retval < 0)
		fprintf(stderr, "WARNING: Failed to set center freq.\n");
	else
		fprintf(stderr, "Tuned to %u Hz.\n", frequency);
	retval = rtlsdr_set_tuner_gain_mode(dev, 1);
	if (retval < 0)
		fprintf(stderr, "WARNING: Failed to enable automatic gain.\n");
	retval = rtlsdr_set_tuner_gain(dev, gain);
	if (retval < 0)
		fprintf(stderr, "WARNING: Failed to set tuner gain.\n");
	else
		fprintf(stderr, "Tuner gain set to %f dB.\n", gain/10.0);
	//Grab some samples
	if(debug)printf("Reading samples in async mode...\n");
	retval = rtlsdr_reset_buffer(dev);
	retval = rtlsdr_read_async(dev, rtlsdr_callback, (void *)file,DEFAULT_ASYNC_BUF_NUMBER, out_block_size);
	rtlsdr_close(0);
	// stop timer
	out_dat();
    time(&stop);
    if(debug)printf("Finished in about %.0f seconds. \n", difftime(stop, start));

	fclose(file1);
	exit(0);
}

/*fast fourier transform routine*/
void four( double data[], int nn,int isign)

{
 int n,mmax,m,j,istep,i,a;
 double wtemp,wr,wpr,wpi,wi,theta;
 double tempr,tempi;
 n=nn<<1;
 j=1;
 for(i=1;i<n;i+=2)
    {
    if(j>i){
    SWAP(data[j],data[i]);
    SWAP(data[j+1],data[i+1]);
    }
    m=n>>1;
    while(m>=2 && j>m){
j-=m;
m>>=1;
}
    j+=m;
    }
 mmax=2;

 while(n>mmax)
   {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr=-2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for(m=1;m<mmax;m+=2){
      for(i=m;i<=n;i+=istep){
 j=i+mmax;
 tempr=wr*data[j]-wi*data[j+1];
 tempi=wr*data[j+1]+wi*data[j];
 data[j]=data[i]-tempr;
 data[j+1]=data[i+1]-tempi;
 data[i]+=tempr;

 data[i+1]+=tempi;
if(j<0)j=0;

     }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
 }
    mmax=istep;
    }
 if(isign==1){for(a=0;a<2*pts;a++){
 data[a]=data[a]/pts;}}
}

void sum_dat(void)
{
int tt;
for(tt=0;tt<pts;tt++){
datr[tt]=datr[tt]+(float) (dats[2*tt]*dats[2*tt]+dats[2*tt+1]*dats[2*tt+1]);
}
}

void out_dat(void)
{
int tt;
float opp;
if(debug)printf("OUT to FILE\n");
for(tt=0;tt<pts;tt++)
	{
	if(tt<(pts/2))
		opp=(float)datr[tt+pts/2];
	else
		opp=(float)datr[tt-pts/2];
	fprintf(file1,"%d    %3.3f\n",(tt-pts/2),(float)(opp/p_num));
	}
}
