/*======================================================*/
/* Author: Chao-wang Huang                                                                  */
/* Date: Wednesday, February 04, 2009                                                */
/* A Bit-based Message Passing Algorithm for ICI Channel is simulated */
/* ICI Power Distribution                                                                         */
/*======================================================*/

#include "stdafx.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <conio.h>
#include <limits.h>
#include <float.h>
#include <cstdlib>
#include <iostream>
using namespace std;
#include <fstream>
#include <cstdio>
#include <assert.h>

void AWGN_noise(float, double, double *);
int Error_count(int, int);
void JakesFading(double, double, double, int, double *, double);
void Multipath_fade_pattern(double ***, int);
//void dfour1(double *, unsigned long, int);

//#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 		3.14159265358979
#define num_packet 10000						// number of packets simulated
const int N_fft =  1024;
const int CP_length = N_fft/8;				// Length of Cyclic Prefix
const double sample_rate = 11.2e6;		/* Sampling Frequency */
/* power Delay Profile of ITU Vehicular A Channel */
const int path_num = 6;
const int path_delay[6] = {0,3,8,12,19,28};	// Delay (sample) {0,310,710,1090,1730,2510} (ns)
const double path_gain[6] = {0.485,0.3853,0.0611,0.0485,0.0153,0.0049};	// Power Gain (linear scale)
const int delay_spread = 28;				// Maximum path delay (samples)
const double vc = 350.0;						/* speed of vehicle in km/hr */
const double C = 3.0e8;						/* light speed */
const double fc = 2.5e9;						/* carrier frequency */
const double OFDM_sym_duration = 91.43e-6;		/* OFDMA Symbol Duration without Cyclic Prefix */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
double *ts, *theta_n;
int *delay;

int main(void)
{
	time_t  ti, start, end;
	int d, i, k, l, p, t;
	double *gain, ***fade, ***F_l, **H_I, **H_Q, *ICI_leak, *ICI_ave;
	FILE *ICI_dist1, *ICI_dist2;

	start = time(NULL);
	printf("ICI Power Distribution in SISO ICI Channel\n");
	printf("Number of packets of simulation = %d\n\n", num_packet);
	printf("This program is running. Don't close, please!\n\n");

	ts = new double[path_num];
	theta_n = new double[path_num];
	gain = new double[path_num];
	delay = new int[path_num];
	fade = new double**[path_num];
	for(l=0; l<path_num; l++)
		fade[l] = new double*[N_fft+CP_length];
	for(l=0; l<path_num; l++)
		for(t=0; t<N_fft+CP_length; t++)
			fade[l][t] = new double[2];
	F_l = new double**[path_num];
	for(i=0; i<path_num; i++)
		F_l[i] = new double*[N_fft];
	for(i=0; i<path_num; i++)
		for(l=0; l<N_fft; l++)
			F_l[i][l] = new double[2];
	H_I = new double*[N_fft];
	for(l=0; l<N_fft; l++)
		H_I[l] = new double[N_fft];
	H_Q = new double*[N_fft];
	for(l=0; l<N_fft; l++)
		H_Q[l] = new double[N_fft];
	ICI_leak = new double[N_fft];
	ICI_ave = new double[N_fft];						

	// Channel Weighting Gain and Channel path delay
	for(i=0; i<path_num; i++)
	{
		delay[i] = path_delay[i];
		gain[i] = sqrt(path_gain[i]);
	}

	srand((unsigned) time(&ti));

/*======================*/
/* main simulation loop            */
/*======================*/
	fopen_s(&ICI_dist1, "ICI_distribution.log", "w");
	fopen_s(&ICI_dist2, "ICI_dist_ave.log", "w");
	// Initial time of Rayleigh fading pattern
	for(l=0; l<path_num; l++)
	{
		ts[l] = 100.0 + l*100.0;
		theta_n[l] = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
	}

	for(t=0; t< N_fft; t++)
		ICI_leak[t] = ICI_ave[t] = 0.0;

	p = 0;
	do
	{
		// Multi-path Rayleigh fading pattern
		//for(l=0; l<path_num; l++)
		//	theta_n[l] = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
		Multipath_fade_pattern(fade, path_num);

		// Channel Frequency Response
		for(l=0; l<path_num; l++)
			for(k=0; k<N_fft; k++)
			{
				F_l[l][k][0] = F_l[l][k][1] = 0.0;
				for(t=0; t<N_fft; t++)
				{
					F_l[l][k][0] += gain[l]*fade[l][t+CP_length-delay[l]][0] * cos(2*Pi*t*k/(float)N_fft) + gain[l]*fade[l][t+CP_length-delay[l]][1] * sin(2*Pi*t*k/(float)N_fft);
					F_l[l][k][1] += gain[l]*fade[l][t+CP_length-delay[l]][1] * cos(2*Pi*t*k/(float)N_fft) - gain[l]*fade[l][t+CP_length-delay[l]][0] * sin(2*Pi*t*k/(float)N_fft);
				}
			}

		// Channel Frequency Response & ICI Channel Coefficients
		for(t=0; t<N_fft; t++)
			for(d=0; d<N_fft; d++)
			{
				H_I[t][d] = H_Q[t][d] = 0.0;
				for(l=0; l<path_num; l++)
				{
					H_I[t][d] += F_l[l][d][0] * cos(2*Pi*delay[l]*(t-d)/(float)N_fft) + F_l[l][d][1] * sin(2*Pi*delay[l]*(t-d)/(float)N_fft);
					H_Q[t][d] += F_l[l][d][1] * cos(2*Pi*delay[l]*(t-d)/(float)N_fft) - F_l[l][d][0] * sin(2*Pi*delay[l]*(t-d)/(float)N_fft);
				}

				H_I[t][d] /= (float)N_fft;
				H_Q[t][d] /= (float)N_fft;
			}

		// ICI Power Distribution (subcarrier 0)
		for(t=0; t< N_fft; t++)
			ICI_leak[t] += (pow(H_I[t][t], 2.0) + pow(H_Q[t][t], 2.0));

		// ICI Power Distribution (Average)
		for(d=0; d<=20; d++)
			for(t=0; t<N_fft; t++)
				ICI_ave[d] += (pow(H_I[t][d], 2.0) + pow(H_Q[t][d], 2.0));

		for(d=1023; d>1023-20; d--)
			for(t=0; t<N_fft; t++)
				ICI_ave[d] += (pow(H_I[t][d], 2.0) + pow(H_Q[t][d], 2.0));
	
		p++;
	} while(p < num_packet);

	for(t=0; t< N_fft; t++)
	{
		ICI_leak[t] /= num_packet;
		ICI_ave[t] /= (num_packet*N_fft);
		fprintf(ICI_dist1, "%f\n", 10*log10(ICI_leak[t]));
		fprintf(ICI_dist2, "%f\n", 10*log10(ICI_ave[t]));
	}
	fclose(ICI_dist1);
	fclose(ICI_dist2);

	delete ts;
	delete theta_n;
	delete gain;
	delete delay;	
	for(l=0; l<path_num; l++)
		for(t=0; t<N_fft+CP_length; t++)
			delete fade[l][t];
	for(l=0; l<path_num; l++)
		delete fade[l];
	delete fade;
	for(i=0; i<path_num; i++)
		for(l=0; l<N_fft; l++)
			delete F_l[i][l];
	for(i=0; i<path_num; i++)
		delete F_l[i];
	delete F_l;
	for(l=0; l<N_fft; l++)
		delete H_I[l];
	delete H_I;
	for(l=0; l<N_fft; l++)
		delete H_Q[l];
	delete H_Q;
	delete ICI_leak;
	delete ICI_ave;

	end = time(NULL);
	printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
	printf("This program is ended. Press any key to continue.\n");
	getchar();

	return 0;
}

void AWGN_noise(float mu, double variance, double *noise)
{
//	const  float Pi = 3.14159265358979;
   double u1, u2;
   do
   {
   	u1 = (double)rand()/(double)RAND_MAX;
      u2 = (double)rand()/(double)RAND_MAX;
   }
   while(u1 == 0.0 || u2 == 0.0);

   *(noise+0) = (sqrt(-2.0*log(u1))*cos(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
   *(noise+1) = (sqrt(-2.0*log(u1))*sin(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
}

int Error_count(int x, int y)
{
	if(x == y)
   	return 0;
   else
   	return 1;
}

void JakesFading(double f_c/*Hz*/, double v/*m/s*/, double t/*s*/, int type, double *fade, double theta_n)
{
	//const double C = 3.0e8;     // (m/s)
   //const float Pi = 3.14159265358979;
   int n, Np, N_o = 32;
   double lamda, w_m, beta_n, w_n, alpha, T_c2, T_s2/*, theta_n*/;

   lamda = C/f_c;     // wave length (meter)
   w_m = 2.0*Pi*v/lamda;    // maximum Doppler frequency
   Np = 2*(2*N_o+1);

   switch(type)
   {
   	case 1:
   		alpha = 0.0;
         T_c2 = (double)N_o;
         T_s2 = (double)N_o + 1.0;
         break;
      case 2:
      	alpha = 0.0;
         T_c2 = (double)N_o + 1.0;
         T_s2 = (double)N_o;
         break;
      case 3:
      	alpha = Pi/4.0;
         T_c2 = (double)N_o + 0.5;
         T_s2 = (double)N_o + 0.5;
         break;
      default:
      	printf("\nInvalid type selection for Jake's fading channel model.\n");
         break;
   }

   if(v < 0.0)
	{
		printf("Warning!! The vehicle speed is invalid.\n");
		cout << "Vehicle speed: " << v << "km/h" << endl;
		getchar();
   		//*(fade+0) = 1.0;
		//*(fade+1) = 0.0;
   }
   else
   {
   	*(fade+0) = sqrt(1.0/T_c2)*cos(alpha)*cos(w_m*t);
      *(fade+1) = sqrt(1.0/T_s2)*sin(alpha)*cos(w_m*t);

      for(n = 1; n <= N_o; n++)
      {
      	switch(type)
         {
         	case 1:
            	beta_n = (double)n*Pi/((double)N_o+1.0);
               break;
            case 2:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
            case 3:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
         	default:
            	break;
         }
         w_n = w_m*cos(2.0*Pi*(double)n/(double)Np);
            //theta_n = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
			//theta_n = 0.0;
         *(fade+0) += sqrt(2.0/T_c2)*cos(beta_n)*cos(w_n*t+theta_n);
         *(fade+1) += sqrt(2.0/T_s2)*sin(beta_n)*cos(w_n*t+theta_n);
		}
	}
}

void Multipath_fade_pattern(double ***fade, int path_num)
{
	int l, t;
	double gain[2];

	for(l=0; l<path_num; l++)
		for(t=0; t<N_fft+CP_length; t++)
		{
			JakesFading(fc, vc*1000/3600.0, ts[l], 2, &gain[0], theta_n[l]);
			ts[l] += OFDM_sym_duration / (float)N_fft;
			fade[l][t][0] = gain[0];
			fade[l][t][1] = gain[1];
		}
}
/*
void dfour1(double data[],unsigned long nn,int isign)
//double data[];
//unsigned long nn;
//int isign;
{
	unsigned long n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;
   double tempr,tempi;

   n=nn << 1;
   j=1;
   for (i=1;i<n;i+=2)
   {
   	if (j > i)
      {
      	SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
      }
      m=n >> 1;
      while (m >= 2 && j > m)
      {
      	j -= m;
         m >>= 1;
      }
      j += m;
   }
   mmax=2;
   while (n > mmax)
   {
   	istep=mmax << 1;
      theta=isign*(6.28318530717959/mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2)
      {
      	for (i=m;i<=n;i+=istep)
         {
         	j=i+mmax;
            tempr=wr*data[j]-wi*data[j+1];
            tempi=wr*data[j+1]+wi*data[j];
            data[j]=data[i]-tempr;
            data[j+1]=data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
      }
      mmax=istep;
   }
}
#undef SWAP
*/