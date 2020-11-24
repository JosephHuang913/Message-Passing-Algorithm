/*========================================================*/
/* Author: Chao-wang Huang                                                                       */
/* Date: Friday, April 23, 2010                                                                    */
/* A Bit-based Message Passing Algorithm for MIMO Channel is simulated */
/* 4 by 4 QPSK                                                                                           */
/* MIMO Detector: Log-domain MPA with Sum-Product Rule                       */
/*========================================================*/

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
void JakesFading(double, double, double, int, double *);
void Multipath_fade_pattern(double ****, int);
void dfour1(double *, unsigned long, int);

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 	3.14159265358979
#define Num_packet 125					// number of packets simulated
const int Mc = 2;								// modulation order (2 for QPSK)
const int T = 4;									// Number of transmit antenna
const int R = 4;									// Number of receive antenna
const int Iteration = 10;
const int N_fft =  1024;
const int CP_length = N_fft/8;				// Length of Cyclic Prefix
const double sample_rate = 11.2e6;		/* Sampling Frequency */
/* Power Delay Profile of ITU Vehicular A Channel */
const int path_num = 6;
const int path_delay[6] = {0, 3, 8, 12, 19, 28};	// Delay (sample) {0,310,710,1090,1730,2510} (ns)
const double path_gain[6] = {0.485, 0.3853, 0.0611, 0.0485, 0.0153, 0.0049};	// Power Gain (linear scale)
const int delay_spread = 28;				// Maximum path delay (samples)
const double vc = 250.0;						/* speed of vehicle in km/hr */
const double C = 3.0e8;						/* light speed */
const double fc = 2.5e9;						/* carrier frequency */
//const double OFDM_sym_duration = 91.43e-6;		/* OFDMA Symbol Duration without Cyclic Prefix */
const double OFDM_sym_duration = 102.86e-6;		/* OFDMA Symbol Duration with Cyclic Prefix */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
const double fdxT = Doppler * OFDM_sym_duration;
double ***ts;
#define N_o		16
double theta_n[N_o];

int main(void)
{
	time_t  ti, start, end;
	int i, j, l, p, r, x, z, t, **data_bit, **Ak;
	int *err_count, *delay, k, n, s, *data;
	double snr, noise_pwr, noise[2], **Y, *err_rate, Ps, Eb_No, ***Yk;
	double P, sum, sum0, sum1, pdf, ****R_pn, ****Q_np, ***R_tmp, **Out;
	double **sym_I, **sym_Q, *I, *Q, *gain, ****fade, ***H_f, *fft;
	FILE *ber, *records;

	start = time(NULL);
	printf("BER Performance of Bit-based Message Passing in MIMO Channel\n");
	printf("Number of transmit antenna is %d\n", T);
	printf("Number of receive antenna is %d\n", R);
	printf("Vehicular speed: %f km/h\n", vc);
	printf("Maximum Doppler frequency: %f\n", Doppler);
	printf("Channel fading rate: %f\n", fdxT);
	printf("Maximum number of bits of simulation = %d\n\n", Mc*T*N_fft*Num_packet);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Record_Bit-based_Message_Passing.log", "a");
	fprintf(records, "BER Performance of Bit-based Message Passing in MIMO Channel\n");
	fprintf(records, "Number of transmit antenna is %d\n", T);
	fprintf(records, "Number of receive antenna is %d\n", R);
	fprintf(records, "Vehicular speed: %f km/h\n", vc);
	fprintf(records, "Maximum Doppler frequency: %f\n", Doppler);
	fprintf(records, "Channel fading rate: %f\n", fdxT);
	fprintf(records, "Maximum number of bits of simulation = %d\n\n", Mc*T*N_fft*Num_packet);
	fprintf(records, "Eb/No     BER\n");
	fflush(records);

	data_bit = new int*[Mc*T];
	for(i=0; i<Mc*T; i++)
		data_bit[i] = new int[N_fft];
	I = new double[T];
	Q = new double[T];
	Y = new double*[R];
	for(i=0; i<R; i++)
   		Y[i] = new double[2];
	R_pn = new double***[N_fft];
	for(t=0; t<N_fft; t++)
		R_pn[t] = new double**[R];
	for(t=0; t<N_fft; t++)
		for(i=0; i<R; i++)
   			R_pn[t][i] = new double*[Mc*T];
	for(t=0; t<N_fft; t++)
		for(i=0; i<R; i++)
   			for(j=0; j<Mc*T; j++)
   				R_pn[t][i][j] = new double[2];
	R_tmp = new double**[R];
	for(i=0; i<R; i++)
   		R_tmp[i] = new double*[Mc*T];
	for(i=0; i<R; i++)
   		for(j=0; j<Mc*T; j++)
   			R_tmp[i][j] = new double[2];
	Q_np = new double***[N_fft];
	for(t=0; t<N_fft; t++)
		Q_np[t] = new double**[Mc*T];
	for(t=0; t<N_fft; t++)
		for(i=0; i<Mc*T; i++)
   			Q_np[t][i] = new double*[R];
	for(t=0; t<N_fft; t++)
		for(i=0; i<Mc*T; i++)
   			for(j=0; j<R; j++)
   				Q_np[t][i][j] = new double[2];
	Out = new double*[Mc*T];
	for(i=0; i<Mc*T; i++)
   		Out[i] = new double[2];
	err_count = new int[Iteration];
	err_rate = new double[Iteration];
	data = new int[Mc*T];
	ts = new double**[R];
	for(i=0; i<R; i++)
		ts[i] = new double*[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			ts[i][j] = new double[path_num];
	gain = new double[path_num];
	delay = new int[path_num];
	sym_I = new double*[T];
	for(i=0; i<T; i++)
		sym_I[i] = new double[N_fft];
	sym_Q = new double*[T];
	for(i=0; i<T; i++)
		sym_Q[i] = new double[N_fft];
	Ak = new int*[Mc*T];
	for(i=0; i<Mc*T; i++)
		Ak[i] = new int[N_fft];
	fade = new double***[R];
	for(i=0; i<R; i++)
		fade[i] = new double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			fade[i][j] = new double*[path_num];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				fade[i][j][l] = new double[2];
	Yk = new double**[R];
	for(i=0; i<R; i++)
   		Yk[i] = new double*[N_fft];
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			Yk[i][l] = new double[2];
	I = new double[T];
	Q = new double[T];
	H_f = new double**[R];
	for(i=0; i<R; i++)
   		H_f[i] = new double*[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_f[i][j] = new double[2*N_fft+1];
	fft = new double[2*N_fft+1];

	// Channel Weighting Gain and Channel path delay
	for(i=0; i<path_num; i++)
	{
		delay[i] = path_delay[i];
		gain[i] = sqrt(path_gain[i]);
	}

	srand((unsigned) time(&ti));

/*======================*/
/* main simulation loop */
/*======================*/
	fopen_s(&ber, "ber_Message_Passing.log", "w");
	for(snr=15; snr<=15; snr+=2)
	{
   		for(s=0; s<Iteration; s++)
   			err_count[s] = 0;

		// noise power calculation
		Eb_No = (double)snr;
		Ps = 1 * T;
		noise_pwr = 0.5*Ps*R/(T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption
		printf("Eb_No = %f\n", Eb_No);

		// Initial time of Rayleigh fading pattern
		for(i=0; i<R; i++)
			for(j=0; j<T; j++)
				for(l=0; l<path_num; l++)
					ts[i][j][l] = 1000.0 + i*1000.0 + j*100.0 + l*10.0;

		p = 0;
		do
		{
			cout << "Packet: " << p << endl;
			for(l=0; l<N_fft; l++)
   				for(i=0; i<T; i++)
	   			{
					// Generate random information bit stream
					if(rand()/(float)RAND_MAX>=0.5)
						data_bit[2*i][l] = 1;
					else
						data_bit[2*i][l] = 0;

					if(rand()/(float)RAND_MAX>=0.5)
						data_bit[2*i+1][l] = 1;
					else
						data_bit[2*i+1][l] = 0;

					// QPSK Mapping
					sym_I[i][l] = (2*data_bit[2*i][l]-1)/sqrt(2.0);
         			sym_Q[i][l] = (2*data_bit[2*i+1][l]-1)/sqrt(2.0);
				}

			// Random Phase of Rayleigh fading pattern
			if(p%10 == 0)
				for(i=0; i<N_o; i++)
					theta_n[i] = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase

			// Multi-path Rayleigh fading pattern without ICI
			Multipath_fade_pattern(fade, path_num);

			// Channel Frequency Response
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
				{
					for(t=0; t<=2*N_fft; t++)
						fft[t] = 0.0;

					// Time Domain Impulse Response
					for(l=0; l<path_num; l++)	
					{
						fft[2*delay[l]+1] = gain[l] * fade[i][j][l][0];
						fft[2*delay[l]+2] = gain[l] * fade[i][j][l][1];
					}

					// FFT (Frequency Response)
					dfour1(&fft[0], N_fft, -1);

					for(t=0; t<N_fft; t++)
					{
						H_f[i][j][2*t] = fft[2*t+1];
						H_f[i][j][2*t+1] = fft[2*t+2];
					}
				}

			// Frequency domain transmission
			for(l=0; l<N_fft; l++)
				for(i=0; i<R; i++)
       			{
					Yk[i][l][0] = Yk[i][l][1] = 0.0;
         			for(j=0; j<T; j++)
            		{
         				Yk[i][l][0] += (sym_I[j][l] * H_f[i][j][2*l] - sym_Q[j][l] * H_f[i][j][2*l+1]);
           				Yk[i][l][1] += (sym_I[j][l] * H_f[i][j][2*l+1] + sym_Q[j][l] * H_f[i][j][2*l]);
           			}

					//AWGN channel
					AWGN_noise(0, noise_pwr, &noise[0]);
	      	   		Yk[i][l][0] += noise[0];
		     		Yk[i][l][1] += noise[1];
   				}

/*==================================================*/
/*  B i t - b a s e d    M e s s a g e    P a s s i n g    A l g o r i t h m  */
/*==================================================*/

			// Initilization of a-priori probability
			for(t=0; t<N_fft; t++)
				for(i=0; i<Mc*T; i++)
   					for(j=0; j<R; j++)
      					for(k=0; k<2; k++)
   							Q_np[t][i][j][k] = 0.5;

			for(s=0; s<Iteration; s++)
			{
				for(t=0; t<N_fft; t++)
				{
	       			for(k=0; k<2; k++)		// Probability of 1 and 0
			    		for(i=0; i<R; i++)
							for(j=0; j<Mc*T; j++)
               				{
			               		R_tmp[i][j][k] = 0.0;
	               				for(x=0; x<pow(2.0,Mc*T-1); x++)
   								{
      	           					l = 0;
         							for(n=0; n<Mc*T; n++)
           	         					if(n!=j)
           								{
                  	      					data[n] = (x>>l) & 1;
                     						l++;
                        				}

									data[j] = k;
   									for(r=0; r<T; r++)
      								{
				                   		I[r] = (2*data[2*r]-1)/sqrt(2.0);
            							Q[r] = (2*data[2*r+1]-1)/sqrt(2.0);
           							}

                   					Y[i][0] = Y[i][1] = 0.0;
									for(r=0; r<T; r++)
									{
				                       	Y[i][0] += (I[r] * H_f[i][r][2*t] - Q[r] * H_f[i][r][2*t+1]);
										Y[i][1] += (I[r] * H_f[i][r][2*t+1] + Q[r] * H_f[i][r][2*t]);
									}

                     				P = 1.0;
                   					for(n=0; n<Mc*T; n++)
				               			if(n!=j)
						     				P *= Q_np[t][n][i][data[n]];
                    
									sum = 0.0;
									for(z=0; z<2; z++)
                   						sum += pow(Yk[i][t][z]-Y[i][z],2.0);

									pdf = exp(-sum/noise_pwr)/(Pi*noise_pwr);
								   	R_tmp[i][j][k] += pdf * P;
                  				}
               				}

					for(i=0; i<R; i++)
            			for(j=0; j<Mc*T; j++)
						{
         					R_pn[t][i][j][1] = R_tmp[i][j][1]/(R_tmp[i][j][0] + R_tmp[i][j][1]);
       						R_pn[t][i][j][0] = R_tmp[i][j][0]/(R_tmp[i][j][0] + R_tmp[i][j][1]);
						}
				}

				for(t=0; t<N_fft; t++)
					for(j=0; j<Mc*T; j++)
           				for(i=0; i<R; i++)
						{
           					sum0 = sum1 = 1.0;
           					for(n=0; n<R; n++)
           						if(n != i)
								{
                   					sum1 *= R_pn[t][n][j][1];
									sum0 *= R_pn[t][n][j][0];
								}

							Q_np[t][j][i][1] = sum1/(sum1 + sum0);
							Q_np[t][j][i][0] = sum0/(sum1 + sum0);
						}

				// Decision Output
				for(t=0; t<N_fft; t++)
				{
					for(j=0; j<Mc*T; j++)
					{
		           		sum0 = sum1 = 1.0;
						for(i=0; i<R; i++)
						{
							sum1 *= R_pn[t][i][j][1];
							sum0 *= R_pn[t][i][j][0];
						}

						Out[j][1] = sum1/(sum1 + sum0);
						Out[j][0] = sum0/(sum1 + sum0);

	       				if(Out[j][1] >= Out[j][0])
		        			Ak[j][t] = 1;
			    		else
				       		Ak[j][t] = 0;

					   	err_count[s] += Error_count(data_bit[j][t], Ak[j][t]);
					}
				}
			}

			p++;
		} while(err_count[Iteration-1] <= 10000 && p < Num_packet);

		// Statistics and records
		cout << "Error Rate = ";
		for(s=0; s<Iteration; s++)
		{
      		err_rate[s] = err_count[s] / (double)(Mc*T*N_fft*p);
      		printf("%e, ", err_rate[s]);
		}
		cout << endl;

		fprintf(ber, "%f ", Eb_No);
		fprintf(records, "%f ", Eb_No);
		for(s=0; s<Iteration; s++)
		{
			fprintf(ber, "%e ", err_rate[s]);
			fprintf(records, "%e ", err_rate[s]);
		}
		fprintf(ber, "\n");
		fprintf(records, "\n");
		fflush(records);
		fflush(ber);
	}
  
	for(i=0; i<Mc*T; i++)
		delete data_bit[i];
	delete data_bit;
	delete data;
	delete I;
	delete Q;
	for(i=0; i<R; i++)
   		delete Y[i];
	delete Y;
	for(t=0; t<N_fft; t++)
		for(i=0; i<R; i++)
   			for(j=0; j<Mc*T; j++)
   				delete R_pn[t][i][j];
	for(t=0; t<N_fft; t++)
		for(i=0; i<R; i++)
   			delete R_pn[t][i];
	for(t=0; t<N_fft; t++)
		delete R_pn[t];
	delete R_pn;
	for(i=0; i<R; i++)
   		for(j=0; j<Mc*T; j++)
   			delete R_tmp[i][j];
	for(i=0; i<R; i++)
   		delete R_tmp[i];
	delete R_tmp;
	for(t=0; t<N_fft; t++)
		for(i=0; i<Mc*T; i++)
   			for(j=0; j<R; j++)
   				delete Q_np[t][i][j];
	for(t=0; t<N_fft; t++)
		for(i=0; i<Mc*T; i++)
   			delete Q_np[t][i];
	for(t=0; t<N_fft; t++)
		delete Q_np[t];
	delete Q_np;
	for(i=0; i<Mc*T; i++)
   		delete Out[i];
	delete Out;
	delete err_count;
	delete err_rate;
	for(i=0; i<Mc*T; i++)
		delete Ak[i];
	delete Ak;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete ts[i][j];
	for(i=0; i<R; i++)
		delete ts[i];
	delete ts;
	delete gain;
	delete delay;
	for(i=0; i<T; i++)
		delete sym_I[i];
	delete sym_I;
	for(i=0; i<T; i++)
		delete sym_Q[i];
	delete sym_Q;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				delete fade[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete fade[i][j];
	for(i=0; i<R; i++)
			delete fade[i];
	delete fade;
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			delete Yk[i][l];
	for(i=0; i<R; i++)
   		delete Yk[i];
	delete Yk;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete H_f[i][j];
	for(i=0; i<R; i++)
		   delete H_f[i];
	delete H_f;
	delete fft;
  
	end = time(NULL);
	printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
	fprintf(records, "Total elapsed time: %.0f(sec)\n\n", difftime(end,start));
	fclose(ber);
	fclose(records);
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

void JakesFading(double f_c/*Hz*/, double v/*m/s*/, double t/*s*/, int type, double *fade)
{
	//const double C = 3.0e8;     // (m/s)
   //const float Pi = 3.14159265358979;
   int n, Np/*, N_o = 8*/;
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

   //if(v == 0.0)
   //{
   //	*(fade+0) = 1.0;
     // *(fade+1) = 0.0;
   //}
   //else
   //{
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
			*(fade+0) += sqrt(2.0/T_c2)*cos(beta_n)*cos(w_n*t+theta_n[n-1]);
			*(fade+1) += sqrt(2.0/T_s2)*sin(beta_n)*cos(w_n*t+theta_n[n-1]);
		}
	//}
}

void Multipath_fade_pattern(double ****fade, int path_num)
{
	int i, j, l;
	double gain[2];

	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
			{
				JakesFading(fc, vc*1000/3600.0, ts[i][j][l], 2, &gain[0]);
				ts[i][j][l] += OFDM_sym_duration;
				fade[i][j][l][0] = gain[0];
				fade[i][j][l][1] = gain[1];
			}
}

void dfour1(double data[], unsigned long nn, int isign)
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
