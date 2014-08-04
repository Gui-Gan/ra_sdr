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
// with code from testrtl.c from http://www.m0dts.co.uk/files/simple_rtlsdr_fft.c
// and RAFFT2.c from http://y1pwe.co.uk/RAProgs/index.html

//Includes
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <time.h>  
#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//Defines
#ifndef _WIN32
#include <unistd.h>
#else
#include <Windows.h>
#include <io.h>
#include <fcntl.h>
#include "getopt/getopt.h"
#endif
#include "rtl-sdr.h"
#define EEPROM_SIZE 256
#define MAX_STR_SIZE 256
#define STR_OFFSET 0x09
#define LONGITUDE -58.45
#define DEF_int_t 1 
#define DEFAULT_SAMPLE_RATE 2.048e6
#define DEF_Freq 30000000
#define DEF_Gain 207
#define FFT 1024
#define PPM 1
#define DEF_debug 0
#define DEFAULT_ASYNC_BUF_NUMBER 32
#define DEFAULT_BUF_LENGTH DEFAULT_SAMPLE_RATE*2
#define MINIMAL_BUF_LENGTH 512
#define MAXIMAL_BUF_LENGTH (256 * 16384)
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
#define PI (6.28318530717959/2.0)

//Global Vars
typedef struct rtlsdr_config {
	uint16_t vendor_id;
	uint16_t product_id;
	char manufacturer[MAX_STR_SIZE];
	char product[MAX_STR_SIZE];
	char serial[MAX_STR_SIZE];
	int have_serial;
	int enable_ir;
	int remote_wakeup;
	} rtlsdr_config_t;

float longitud=LONGITUDE;
float yy,rr,rrr,rav,amp;
int xx,r,pts,a,int_t=DEF_int_t,debug=DEF_debug;
int flag=0;//,da[FFTs*2];
long int file_end,count,p_num;
double sig_pow,nois_pow; //dats[FFTs*4],datr[FFTs*4],
int retval,opt,aux_rep=1;
int devices;
static rtlsdr_dev_t *dev;
int n;
uint8_t *buf;   //unsigned 8bit int - I didn't know what it was!, the _t must be 'type'
double *buf2;
long double *fft_out;
char * fftresult;
uint32_t frequency = DEF_Freq,aux_frequency=0;
int gain = DEF_Gain,aux_gain=0;
static uint32_t bytes_to_read = DEFAULT_SAMPLE_RATE;
uint32_t out_block_size = DEFAULT_BUF_LENGTH;
static int do_exit = 0;
float sample_rate_aux=DEFAULT_SAMPLE_RATE,aux_sample_rate=0;
FILE *file;
char *filename = NULL;
FILE *file2;
char *filename2 = NULL;
FILE *file_raw;
char *filename_raw = NULL;
int ppm_aux=PPM;
time_t start, start2, stop, stop2;
int  s_hour, s_min, s_sec;
double s_t;
int n_read,reps=1,len=0,fft_aux=FFT;

//Functions
void zenith_sideraltime(void);
void tpow(void);
void dump_config(rtlsdr_config_t *conf);
void dump_config(rtlsdr_config_t *conf);
int get_string_descriptor(int pos, uint8_t *data, char *str);
int parse_eeprom_to_conf(rtlsdr_config_t *conf, uint8_t *dat);
void usage(void);


int main(int argc, char **argv)
	{
	uint8_t buf1[EEPROM_SIZE];
	rtlsdr_config_t conf;
	int mr,mi;
	float f_init,f_step;
	int raw_out=0,fft_rep=0,aux=0,i=0,ii=0,n=0,tt=0;
	float bw_channel=0.0,t_channel=0.0;
	uint32_t dev_index = 0;
	time_t now = time(NULL);
	struct tm *t = localtime(&now);
	while ((opt = getopt(argc, argv, "d:f:g:s:i:p:r::t:v::")) != -1) {
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
			reps=atoi(optarg);
			break;
        case 'p':
            ppm_aux=atoi(optarg);
            break;
		case 'r':
			raw_out=1;
			break;
		case 't':
		   fft_aux=atoi(optarg);
		   break;
		case 'v':
		   debug=1;
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
	if(reps <1){
		usage();		
	}
	aux_rep=reps;
	time(&start);
//File handlers	
	len = (strlen(filename));
	if(debug)printf("file 1 %s\n",filename);
	filename2 = malloc(sizeof(char) * (len+4));
	memcpy(filename2, filename,len-4);
	filename2[len-4] = 0;
	strcat(filename2,"_TP.txt");
	if(debug)printf("file 2 %s\n",filename2);
        file = fopen(filename, "w");
        if (!file)
        {
                        fprintf(stderr, "Failed to open %s\n", filename);
                        return(1);
        }

	file2 = fopen(filename2, "w");
        if (!file2)
        {
                        fprintf(stderr, "Failed to open %s\n", filename2);
                        return(1);
        }
	if(raw_out)
	{
		filename_raw = malloc(sizeof(char) * (len+5));
		memcpy(filename_raw, filename,len-4);
		filename_raw[len-4] = 0;
		strcat(filename_raw,"_raw.bin");
		if(debug)printf("file 3 %s\n",filename_raw);
	        file_raw = fopen(filename_raw, "wb");
	        if (!file_raw)
	        {
                        fprintf(stderr, "Failed to open %s\n", filename_raw);
                        return(1);
	        }
	}
//RTL device checks and config	
	devices = rtlsdr_get_device_count();
	for(n=0;n<devices;n++)
	{
		if(debug)printf("Device %d: %s\n\n",devices,rtlsdr_get_device_name(n));
	}
	if (devices>0){
		retval = rtlsdr_open(&dev, dev_index);	
		if(debug)printf("Open= %d\n",dev_index);
	}else{
		if(debug)printf("No Devices found...!\n\n");
		exit(0);
	}
	retval = rtlsdr_read_eeprom(dev, buf1, 0, EEPROM_SIZE);
	if (retval < 0) 
	{
		if (retval == -3)
			{
			if(debug)printf("No EEPROM has been found.\n");
			}
			else
			{
			if(debug)printf("Failed to read EEPROM, err %i.\n", r);
			}
		exit(0);
	}
	parse_eeprom_to_conf(&conf, buf1);
	if(debug)dump_config(&conf);
//configure rtlsdr settings
	retval = rtlsdr_set_freq_correction(dev, ppm_aux);
        if (retval < 0)
                fprintf(stderr, "WARNING: Failed to set PPM..\n");
        else
                fprintf(stderr, "PPM set to %i .\n", ppm_aux);

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
		fprintf(stderr, "WARNING: Failed to set tuner gain..\n");
	else
		fprintf(stderr, "Tuner gain set to %f dB.\n", gain/10.0);
//Get configuration
	aux_frequency=rtlsdr_get_center_freq(dev);
	aux_gain=rtlsdr_get_tuner_gain(dev);	
	aux_sample_rate=rtlsdr_get_sample_rate(dev);
	ppm_aux=rtlsdr_get_freq_correction(dev);
	zenith_sideraltime(); //Get actual sideral time
//write file headers. RAW file dont contain header
//Original PSD file	
	fprintf(file,"#FRE[Hz] %zu\n",aux_frequency );
	fprintf(file,"#GAIN[dB] %f\n",aux_gain/10.0 );
	fprintf(file,"#SAMP[Hz] %f\n",aux_sample_rate );
	fprintf(file,"#DATE %d%02d%02d\n",t->tm_year+1900,t->tm_mon+1,t->tm_mday);
	fprintf(file,"#TIME %02d%02d%02d\n",t->tm_hour,t->tm_min,t->tm_sec);
	fprintf(file,"#REPS[n] %d\n",reps);
	fprintf(file,"#FFT %d\n",fft_aux);
	fprintf(file,"#LST[ZENITH] %02d%02d%02d\n",s_hour,s_min,s_sec);
	fprintf(file,"#LSA[ZENITH] %03.03f\n",s_t);    
	fprintf(file,"#PPM %d\n",ppm_aux);
	fprintf(file,"#SERIAL %s\n",conf.serial);
//Total Power File	
	fprintf(file2,"#FREQ[Hz] %zu\n",aux_frequency );
	fprintf(file2,"#GAIN[dB] %f\n",aux_gain/10.0 );
	fprintf(file2,"#SAMP[Hz] %f\n",aux_sample_rate );
	fprintf(file2,"#DATE %d%02d%02d\n",t->tm_year+1900,t->tm_mon+1,t->tm_mday);
	fprintf(file2,"#TIME %02d%02d%02d\n",t->tm_hour,t->tm_min,t->tm_sec);
	fprintf(file2,"#REPS[n] %d\n",reps);
	fprintf(file2,"#LST[ZENITH] %02d%02d%02d\n",s_hour,s_min,s_sec);
	fprintf(file2,"#LSA[ZENITH] %03.03f\n",s_t);
	fprintf(file2,"#PPM %d\n",ppm_aux);
	fprintf(file2,"#SERIAL %s\n", conf.serial);
	fprintf(stderr, "Reading samples in sync mode...\n");
	do_exit=0;
	if(debug)printf("reps %d\n",reps);
	while (!do_exit) 
		{
		bw_channel=sample_rate_aux/fft_aux;
		t_channel=1/bw_channel;
		fft_rep=(int)ceil(DEF_int_t/t_channel)*fft_aux;
		if(debug)printf("FFT value is %d real reps %d\n",fft_aux,fft_rep);
		buf = malloc(((fft_rep*2)*reps) * sizeof(uint8_t));	
		time(&start2);
		retval = rtlsdr_reset_buffer(dev);
		retval = rtlsdr_read_sync(dev, buf, ((fft_rep*2)*reps), &n_read);
		time(&stop2);
		if(debug)printf("Finished reads in about %.0f seconds. \n", difftime(stop2, start2));	
		if(debug)printf("Sample_rate_aux %f bytes_to_read %d block size %zu n_reads %d\n",sample_rate_aux,bytes_to_read,out_block_size,n_read);
		time(&start2);
		tpow();
		time(&stop2);
		if(debug)printf("Finished TP write in about %.0f seconds. \n", difftime(stop2, start2));
		if (retval < 0) 
			{
			fprintf(stderr, "WARNING: sync read failed.\n");
			break;
			}
		if ((bytes_to_read > 0) && (bytes_to_read < (uint32_t)n_read)) 
			{
				if(debug)printf("End sync read\n");
				n_read = bytes_to_read;
				do_exit = 1;
				if(raw_out)
				{
					time(&start2);
					if (fwrite(buf, 1, n_read, file_raw) != (size_t)n_read) 
						{
						fprintf(stderr, "Short write, samples lost, exiting!\n");
						break;
						}
					time(&stop2);
								if(debug)printf("Finished RAW write in about %.0f seconds. \n", difftime(stop2, start2));
				}
				//New FFT
				//Configure FFTW to convert the samples in time domain to frequency domain
				fftw_plan my_plan;
				fftw_complex *in, *out;
				in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fft_aux);
				out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fft_aux);
				my_plan = fftw_plan_dft_1d(fft_aux, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
				//convert buffer from IQ to complex ready for FFTW, seems that rtlsdr outputs IQ data as IQIQIQIQIQIQ so ......
				if(debug)printf("n_reads %d\n",n_read);
				aux=sample_rate_aux/fft_aux;
				if(debug)printf("Repetitions for FFT %d \n",aux);
				buf2= malloc((fft_aux*aux*reps) * sizeof(double));
				n=0;
				tt=0;
				for(i=0;i<aux*reps;i++)
					{
					for (ii=0;ii<(fft_aux);ii++)
						{
						__real__  in[ii] = (buf[n]-127.5)/128.0;   //sample is 127 for zero signal,so 127 +/-127  //buf[n]-127
						__imag__ in[ii] = (buf[n+1]-127.5)/128.0;   //aux1[i]=(aux1[i]-127.5)/128.0;  			//buf[n]-127
						n=n+2;
						}
					//Convert the comples samples to complex frequency domain
					fftw_execute(my_plan);
					for (ii=0;ii<fft_aux;ii++)
						{
						mr = __real__ out[ii] * __real__ out[ii];
						mi = __imag__ out[ii] * __imag__ out[ii];
						buf2[tt] = sqrt(mr+mi);
						tt++;
						}						
					}
				if(debug)printf("FFT loops i %d ii %d n %d tt %d\n",i,ii,n,tt);
				fft_out=malloc((fft_aux) * sizeof(long double));
				tt=0;
				for(i=0;i<fft_aux;i++)
					{
					for(ii=0;ii<aux*reps;ii++)
						{
						fft_out[i]=fft_out[i]+buf2[((ii)*fft_aux)+i];
						}
					fft_out[i]=fft_out[i]/(aux*reps);
					}
				f_init=(aux_frequency)-(sample_rate_aux/2); //Frequency index for file write.
				f_step=sample_rate_aux/fft_aux;
				if(debug)printf("finit %f Fstep %f \n",f_init,f_step);
				for(i=fft_aux/2;i<fft_aux;i++)
					{
					if(finite(fft_out[i]))
						{
						fprintf(file,"%3.6f\t%3.3Lf\n",f_init/1e6,fft_out[i]); 
						}
					else
						{
						fprintf(file,"%3.6f\t%3.3Lf\n",f_init/1e6,fft_out[i+1]); 
						}
					f_init=f_init+f_step;
					}
				for(i=0;i<fft_aux/2;i++)
					{
					if(finite(fft_out[i]))
						{
						fprintf(file,"%3.6f\t%3.3Lf\n",f_init/1e6,fft_out[i]); 
						}
					else
						{
						fprintf(file,"%3.6f\t%3.3Lf\n",f_init/1e6,fft_out[i-1]); 
						}
					f_init=f_init+f_step;
					}
			}
		}
	rtlsdr_close(0);
	// stop timer
	free (buf);
	free (buf2);
    time(&stop);
    if(debug)printf("Finished in about %.0f seconds. \n", difftime(stop, start));
//Close Files
	fclose(file);
	fclose(file2);
	if(raw_out)fclose(file_raw);
	exit(0);
}

void zenith_sideraltime(void)
{
	int year, month, day;
	int  hour, min, sec;
	int TZ=-3;
	//Long + at East and - at West of Greenwich
	double julian;
	float jdsince1900;
	double t, tsg;
	//const 	jd1900jan0 = 1900 JAN 0 Julian date
	double jd1900jan0= 2415020.0;
	//GMT Time
	time_t now=time(NULL);
	struct tm *tm=localtime(&now);
	year = tm->tm_year + 1900;
	month =tm->tm_mon + 1;
	day= tm->tm_mday;
	hour=(tm->tm_hour);
	min=tm->tm_min;
	sec= tm->tm_sec ;
	hour=(tm->tm_hour)-TZ;
	if (hour > 23 )
		hour=hour-24;
	if (hour <3)
		day=day+1;
	// Julian day at 0h GMT of Greenwich
	julian = (4712+year)*365.25;
	if (julian == floor(julian))
	julian = julian-1;
	else
	julian = floor(julian);
	julian=julian-13;
	//Day of Year
	julian=julian + day;
	if (month > 1)
	julian=julian+31;
	if (month > 2)
	julian=julian+28;
	if (month > 3)
	julian=julian+31;
	if (month > 4)
	julian=julian+30;
	if (month > 5)
	julian=julian+31;
	if (month > 6)
	julian=julian+30;
	if (month > 7)
	julian=julian+31;
	if (month > 8)
	julian=julian+31;
	if (month > 9)
	julian=julian+30;
	if (month > 10)
	julian=julian+31;
	if (month > 11)
	julian=julian+30;
	if(year == (floor(year/4)*4) && (month > 1))
	julian=julian + 1;
	julian=julian - .5;
	jdsince1900=julian-jd1900jan0;
	t=jdsince1900/36525;
	tsg=6.64606 + 2400.05126 * t + 2.58055e-5 * t * t;
	tsg=tsg+longitud/15.0;
	tsg = tsg + (hour+min/60.0+sec/3600.0)*1.002737909265;
	t = tsg - floor(tsg/24)*24;
	s_hour = floor(t);
	s_min = floor((t-s_hour)*60);
	s_sec = ((t - s_hour - s_min/60.0)*3600);
	s_t=t*15;
	if(debug)printf("\n Local sidereal time %02d h, %02d m, %02d  s", s_hour, s_min, s_sec);
	if(debug)printf("\n Local sidereal Anlge %03.03f ", s_t);
	if(debug)printf("\n ");
/* end */
}

void tpow(void)
{
int i=0,ii=0,n=0;
double long avg_tp=0;
double *aux1;   
aux1 = malloc((out_block_size/reps) * sizeof(double));
double *aux2;   
aux2 = malloc(((out_block_size/reps)/2) * sizeof(double));
for(ii=0;ii<reps;ii++)
	{
	for(i=0;i<out_block_size/reps;i++)
		{
		aux1[i]=(float)buf[n]+0.0;
		if(aux1[i]>127)
				aux1[i]=(aux1[i]-127.5)/128.0;
			else
					aux1[i]=(aux1[i]-127.5)/128.0;
		n++;
		}
	for(i=0;i<(out_block_size/reps)/2;i++)
		{
		aux2[i]=sqrt(pow(aux1[i],2)+pow(aux1[i+1],2));
		avg_tp=avg_tp+aux2[i];
		i++;
		}
	avg_tp=avg_tp/((out_block_size/reps)/2);
	if(debug)printf("Total power relative value avg %Lf - %d of %d\n",avg_tp,i,out_block_size);
	if(finite(avg_tp))
			fprintf(file2,"%Lf\n",avg_tp);
	}
free(aux1);
free(aux2);
}

void dump_config(rtlsdr_config_t *conf)
{
fprintf(stderr, "__________________________________________\n");
fprintf(stderr, "Vendor ID:\t\t0x%04x\n", conf->vendor_id);
fprintf(stderr, "Product ID:\t\t0x%04x\n", conf->product_id);
fprintf(stderr, "Manufacturer:\t\t%s\n", conf->manufacturer);
fprintf(stderr, "Product:\t\t%s\n", conf->product);
fprintf(stderr, "Serial number:\t\t%s\n", conf->serial);
fprintf(stderr, "Serial number enabled:\t");
fprintf(stderr, conf->have_serial ? "yes\n": "no\n");
fprintf(stderr, "IR endpoint enabled:\t");
fprintf(stderr, conf->enable_ir ? "yes\n": "no\n");
fprintf(stderr, "Remote wakeup enabled:\t");
fprintf(stderr, conf->remote_wakeup ? "yes\n": "no\n");
fprintf(stderr, "__________________________________________\n");
}

int parse_eeprom_to_conf(rtlsdr_config_t *conf, uint8_t *dat)
{
int pos;
if ((dat[0] != 0x28) || (dat[1] != 0x32))
fprintf(stderr, "Error: invalid RTL2832 EEPROM header!\n");
conf->vendor_id = dat[2] | (dat[3] << 8);
conf->product_id = dat[4] | (dat[5] << 8);
conf->have_serial = (dat[6] == 0xa5) ? 1 : 0;
conf->remote_wakeup = (dat[7] & 0x01) ? 1 : 0;
conf->enable_ir = (dat[7] & 0x02) ? 1 : 0;
pos = get_string_descriptor(STR_OFFSET, dat, conf->manufacturer);
pos = get_string_descriptor(pos, dat, conf->product);
get_string_descriptor(pos, dat, conf->serial);
return 0;
}

int get_string_descriptor(int pos, uint8_t *data, char *str)
{
int len, i, j = 0;
len = data[pos];
if (data[pos + 1] != 0x03)
fprintf(stderr, "Error: invalid string descriptor!\n");
for (i = 2; i < len; i += 2)
str[j++] = data[pos + i];
str[j] = 0x00;
return pos + i;
}

void usage(void)
{
	fprintf(stderr,
		"ra_sdr, an I/Q recorder for RTL2832 based DVB-T receivers, for RadioAstronomy use. Ver 2.0, Sync Mode Only, Not Fully Tested.\n\n"
		"Usage:\t -f frequency_to_tune_to (default 30000000 [Hz)]\n"
		"\t[-s samplerate (int default: 2048000 Hz)]\n"
		"\t[-d device_index (int default: 0)]\n"
		"\t[-g gain (float default: 20.7)]\n"
		"\t[-p PPM set (int default: 0)]\n"
		"\t[-r RAW outpur to file with extension .bin]\n"
		"\t[-i 'kind of' integration Time in seconds (int default: 1, steps of 1sec)]\n"
		"\t[-t FFT value pow of 2 1024-? (int default: 1024)]\n"
        "\t[-v vervose ]\n"
		"\tfilename \n\n");
	exit(1);
}
