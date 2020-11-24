/*========================================================*/
/* Author: Chao-wang Huang                                                                       */
/* Date: Friday, April 23, 2010                                                                    */
/* A Bit-based Message Passinlg Algorithm for MIMO Channel is simulated */
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

void AWGN_noise(float, long double, long double *);
int Error_count(int, int);
void JakesFading(long double, long double, long double, int, long double *);
void Multipath_fade_pattern(long double ****, int);
void dfour1(long double *, unsigned long, int);

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 	3.14159265358979
#define Num_packet 125					// number of packets simulated
//const int Mc = 4;								// modulation order (4 for 16QAM)
const int Mc = 2;								// modulation order (2 for QPSK)
const int T = 4;									// Number of transmit antenna
const int R = 4;									// Number of receive antenna
const int Iterate = 10;
const int N_fft =  1024;
const int CP_length = N_fft/8;				// Length of Cyclic Prefix
const long double sample_rate = 11.2e6;		/* Sampling Frequency */
/* Power Delay Profile of ITU Vehicular A Channel */
const int path_num = 6;
const int path_delay[6] = {0, 3, 8, 12, 19, 28};	// Delay (sample) {0,310,710,1090,1730,2510} (ns)
const long double path_gain[6] = {0.485, 0.3853, 0.0611, 0.0485, 0.0153, 0.0049};	// Power Gain (linear scale)
const int delay_spread = 28;				// Maximum path delay (samples)
const long double vc = 250.0;						/* speed of vehicle in km/hr */
const long double C = 3.0e8;						/* light speed */
const long double fc = 2.5e9;						/* carrier frequency */
//const long double OFDM_sym_duration = 91.43e-6;		/* OFDMA Symbol Duration without Cyclic Prefix */
const long double OFDM_sym_duration = 102.86e-6;		/* OFDMA Symbol Duration with Cyclic Prefix */
const long double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
const long double fdxT = Doppler * OFDM_sym_duration;
long double ***ts;

int main(void)
{
	time_t  ti, start, end;
	int i, j, q, l, m, p, r, s1, x, z, t, d, **data_bit, **Ak;
	int *err_count, *delay, *data;
	long double Eb_No, snr, noise_pwr, noise[2], Y[2], ***Yk, *err_rate, Ps;
	long double sum, pdf, ***L_r_pn, ***L_q_np, **sym_I, **sym_Q, *I, *Q;
	long double *gain, ****fade, ***H_f, *fft;
	long double **L_q, L_sum, L_tmp[2];
	FILE *ber, *records;

	start = time(NULL);
	printf("BER Performance of Bit-based Message Passinlg in MIMO Channel\n");
	printf("Number of transmit antenna is %d\n", T);
	printf("Number of receive antenna is %d\n", R);
	printf("Vehicular speed: %f km/h\n", vc);
	printf("Maximum Doppler frequency: %f\n", Doppler);
	printf("Channel fading rate: %f\n", fdxT);
	printf("Maximum number of bits of simulation = %d\n\n", Mc*T*N_fft*Num_packet);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Record_Bit-based_Message_Passing.log", "a");
	fprintf(records, "BER Performance of Bit-based Message Passinlg in MIMO Channel\n");
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
	data = new int[Mc*T];
	ts = new long double**[R];
	for(i=0; i<R; i++)
		ts[i] = new long double*[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			ts[i][j] = new long double[path_num];
	gain = new long double[path_num];
	delay = new int[path_num];
	sym_I = new long double*[T];
	for(i=0; i<T; i++)
		sym_I[i] = new long double[N_fft];
	sym_Q = new long double*[T];
	for(i=0; i<T; i++)
		sym_Q[i] = new long double[N_fft];
	Ak = new int*[Mc*T];
	for(i=0; i<Mc*T; i++)
		Ak[i] = new int[N_fft];
	fade = new long double***[R];
	for(i=0; i<R; i++)
		fade[i] = new long double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			fade[i][j] = new long double*[path_num];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				fade[i][j][l] = new long double[2];
	Yk = new long double**[R];
	for(i=0; i<R; i++)
   		Yk[i] = new long double*[N_fft];
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			Yk[i][l] = new long double[2];
	L_r_pn = new long double**[N_fft];
	for(q=0; q<N_fft; q++)
		L_r_pn[q] = new long double*[R];
	for(q=0; q<N_fft; q++)
		for(i=0; i<R; i++)
   			L_r_pn[q][i] = new long double[Mc*T];
	L_q_np = new long double**[N_fft];
	for(q=0; q<N_fft; q++)
		L_q_np[q] = new long double*[Mc*T];
	for(q=0; q<N_fft; q++)
		for(i=0; i<Mc*T; i++)
   			L_q_np[q][i] = new long double[R];
	L_q = new long double*[N_fft];
	for(q=0; q<N_fft; q++)
		L_q[q] = new long double[Mc*T];
	err_count = new int[Iterate];
	err_rate = new long double[Iterate];
	I = new long double[T];
	Q = new long double[T];
	H_f = new long double**[R];
	for(i=0; i<R; i++)
   		H_f[i] = new long double*[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_f[i][j] = new long double[2*N_fft+1];
	fft = new long double[2*N_fft+1];

	// Channel Weighting Gain and Channel path delay
	for(i=0; i<path_num; i++)
	{
		delay[i] = path_delay[i];
		gain[i] = sqrtl(path_gain[i]);
	}

	srand((unsigned) time(&ti));

/*======================*/
/* main simulation loop            */
/*======================*/
	fopen_s(&ber, "BER_MPA.log", "w");
	for(snr=0; snr<=16; snr+=2)
	{
   		for(s1=0; s1<Iterate; s1++)
   			err_count[s1] = 0;

   		// noise power calculation
		Eb_No = (double)snr;
		Ps = 1 * T;
		noise_pwr = 0.5*Ps*R/(T*powl(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption
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
					sym_I[i][l] = (2*data_bit[2*i][l]-1)/sqrtl(2.0);
         			sym_Q[i][l] = (2*data_bit[2*i+1][l]-1)/sqrtl(2.0);
/*
					// 16-QAM Mapping	
					if(data_bit[l][Mc*i] == 0 && data_bit[l][Mc*i+2] == 0)
   		  				sym_I[i][l] = 1.0/sqrtl(10.0);
					else if(data_bit[l][Mc*i] == 0 && data_bit[l][Mc*i+2] == 1)
   			  			sym_I[i][l] = 3.0/sqrtl(10.0);
					else if(data_bit[l][Mc*i] == 1 && data_bit[l][Mc*i+2] == 1)
		  		 		sym_I[i][l] = -3.0/sqrtl(10.0);
				  	else
  		   				sym_I[i][l] = -1.0/sqrtl(10.0);

		   			if(data_bit[l][Mc*i+1] == 0 && data_bit[l][Mc*i+3] == 0)
  	   					sym_Q[i][l] = 1.0/sqrtl(10.0);
					else if(data_bit[l][Mc*i+1] == 0 && data_bit[l][Mc*i+3] == 1)
   				  		sym_Q[i][l] = 3.0/sqrtl(10.0);
				     else if(data_bit[l][Mc*i+1] == 1 && data_bit[l][Mc*i+3] == 1)
	   			 		sym_Q[i][l] = -3.0/sqrtl(10.0);
					 else
      		 			sym_Q[i][l] = -1.0/sqrtl(10.0);*/
				}

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
       			}
			
         	//AWGN channel
         	for(l=0; l<N_fft; l++)
				for(i=0; i<R; i++)
				{
   					AWGN_noise(0, noise_pwr, &noise[0]);
	      	   		Yk[i][l][0] += noise[0];
		     		Yk[i][l][1] += noise[1];
			    }

/*========================================================*/
/*  Log-Domain Bit-based Message Passinlg Algorithm                               */
/*========================================================*/
			
			// MIMO Detector: LLR Initialization
			for(t=0; t<N_fft; t++)
				for(j=0; j<T; j++)
					for(m=0; m<Mc; m++)
						for(i=0; i<R; i++)
							L_q_np[t][j*Mc+m][i] = 0.0;

				for(s1=0; s1<Iterate; s1++)
				{
					for(t=0; t<N_fft; t++)
					{
						// Extrinsic information of function nodes
			    		for(i=0; i<R; i++)
							for(j=0; j<T; j++)
								for(m=0; m<Mc; m++)
           						{
									for(q=0; q<2; q++)		// Probability of 1 and 0
									{
               							L_tmp[q] = 0.0;
		           						for(x=0; x<powl(2.0,Mc*T-1.0); x++)
	   									{
		               						l = 0;
			  								for(r=0; r<Mc*T; r++)
				   		         				if(r != j*Mc+m)
       											{
							  			  			data[r] = (x>>l) & 1;
							 						l++;
												}

											data[j*Mc+m] = q;
			   								for(r=0; r<T; r++)
				  							{
												// QPSK Mapping
												I[r] = (2*data[r*Mc]-1)/sqrtl(2.0);
         										Q[r] = (2*data[r*Mc+1]-1)/sqrtl(2.0);
/*
												// 16-QAM Mapping	
												if(data[Mc*r] == 0 && data[Mc*r+2] == 0)
							   		  				I[r] = 1.0/sqrtl(10.0);
												else if(data[Mc*r] == 0 && data[Mc*r+2] == 1)
   				  									I[r] = 3.0/sqrtl(10.0);
												else if(data[Mc*r] == 1 && data[Mc*r+2] == 1)
   		 											I[r] = -3.0/sqrtl(10.0);
			  									else
  		   											I[r] = -1.0/sqrtl(10.0);
	
						   						if(data[Mc*r+1] == 0 && data[Mc*r+3] == 0)
  					   								Q[r] = 1.0/sqrtl(10.0);
												else if(data[Mc*r+1] == 0 && data[Mc*r+3] == 1)
   					  								Q[r] = 3.0/sqrtl(10.0);
											     else if(data[Mc*r+1] == 1 && data[Mc*r+3] == 1)
							   				 		Q[r] = -3.0/sqrtl(10.0);
												 else
		       										Q[r] = -1.0/sqrtl(10.0);*/
		     								}
	
					         				Y[0] = Y[1] = 0.0;
											for(r=0; r<T; r++)
											{
												Y[0] += (I[r] * H_f[i][r][2*t] - Q[r] * H_f[i][r][2*t+1]);
               									Y[1] += (I[r] * H_f[i][r][2*t+1] + Q[r] * H_f[i][r][2*t]);
											}
										
								 			L_sum = 0.0;
											if( s1 != 0)
	                     						for(r=0; r<T; r++)
													for(d=0; d<Mc; d++)
			             								if( (r*Mc+d != j*Mc+m) && (data[r*Mc+d] == 1) )
				       										L_sum += L_q_np[t][r*Mc+d][i];
																
											if (L_sum < -50.0)
												L_sum = -50.0;
											if (L_sum > 50.0)
												L_sum = 50.0;

											sum = 0.0;
											for(z=0; z<2; z++)
							        			sum += powl(Yk[i][t][z]-Y[z],2.0);

											pdf = expl(-sum/noise_pwr)/(Pi*noise_pwr);
								     		L_tmp[q] += (pdf * expl(L_sum));
				          				}
									}

									L_r_pn[t][i][j*Mc+m] = logl(L_tmp[1]/L_tmp[0]);

									if (L_r_pn[t][i][j*Mc+m] < -50.0)
										L_r_pn[t][i][j*Mc+m] = -50.0;
									if (L_r_pn[t][i][j*Mc+m] > 50.0)
										L_r_pn[t][i][j*Mc+m] = 50.0;
								}
					}

					// Extrinsic information of bit nodes
					for(t=0; t<N_fft; t++)
						for(j=0; j<T; j++)
							for(m=0; m<Mc; m++)
							{
								for(i=0; i<R; i++)
								{
									L_q_np[t][j*Mc+m][i] = 0.0;
				   					for(r=0; r<R; r++)
		       							if(r != i)
											L_q_np[t][j*Mc+m][i] += L_r_pn[t][r][j*Mc+m];
										
									if (L_q_np[t][j*Mc+m][i] < -50.0)
										L_q_np[t][j*Mc+m][i] = -50.0;
									if (L_q_np[t][j*Mc+m][i] > 50.0)
										L_q_np[t][j*Mc+m][i] = 50.0;
				 				}
							}
		
					// For Decision
					for(t=0; t<N_fft; t++)
						for(j=0; j<T; j++)
							for(m=0; m<Mc; m++)
							{
								L_q[t][j*Mc+m] = 0.0;
				   				for(r=0; r<R; r++)		
									L_q[t][j*Mc+m] += L_r_pn[t][r][j*Mc+m];
								
								if(L_q[t][j*Mc+m] >= 0.0)
									Ak[j*Mc+m][t] = 1;
								else
									Ak[j*Mc+m][t] = 0;

								err_count[s1] += Error_count(data_bit[j*Mc+m][t], Ak[j*Mc+m][t]);
				 			}
				}

			p++;
		} while(err_count[Iterate-1] <= 10000 && p < Num_packet);

		// Statistics and records
		cout << "Error Rate = ";
		for(s1=0; s1<Iterate; s1++)
		{
      		err_rate[s1] = err_count[s1] / (double)(Mc*T*N_fft*p);
      		printf("%e, ", err_rate[s1]);
		}
		cout << endl;

		fprintf(ber, "%f ", Eb_No);
		fprintf(records, "%f ", Eb_No);
		for(s1=0; s1<Iterate; s1++)
		{
			fprintf(ber, "%e ", err_rate[s1]);
			fprintf(records, "%e ", err_rate[s1]);
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
	for(q=0; q<N_fft; q++)
		for(i=0; i<R; i++)
   			delete L_r_pn[q][i];
	for(q=0; q<N_fft; q++)
   		delete L_r_pn[q];
	delete L_r_pn;
	for(q=0; q<N_fft; q++)
		for(i=0; i<Mc*T; i++)
   			delete L_q_np[q][i];
	for(q=0; q<N_fft; q++)
		delete L_q_np[q];
	delete L_q_np;
	for(q=0; q<N_fft; q++)
		delete L_q[q];
	delete L_q;
	delete err_count;
	delete err_rate;
	delete I;
	delete Q;
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

void AWGN_noise(float mu, long double variance, long double *noise)
{
//	const  float Pi = 3.14159265358979;
   long double u1, u2;
   do
   {
   	u1 = (long double)rand()/(long double)RAND_MAX;
      u2 = (long double)rand()/(long double)RAND_MAX;
   }
   while(u1 == 0.0 || u2 == 0.0);

   *(noise+0) = (sqrtl(-2.0*logl(u1))*cosl(2*Pi*u2))*sqrtl(variance/2.0)+mu/sqrtl(2.0);
   *(noise+1) = (sqrtl(-2.0*logl(u1))*sinl(2*Pi*u2))*sqrtl(variance/2.0)+mu/sqrtl(2.0);
}

int Error_count(int x, int y)
{
	if(x == y)
   	return 0;
   else
   	return 1;
}

void JakesFading(long double f_c/*Hz*/, long double v/*m/s*/, long double t/*s*/, int type, long double *fade)
{
	//const double C = 3.0e8;     // (m/s)
   //const float Pi = 3.14159265358979;
   int n, Np, N_o = 32;
   double lamda, w_m, beta_n, w_n, alpha, T_c2, T_s2, theta_n;

   lamda = C/f_c;     // wave length (meter)
   w_m = 2.0*Pi*v/lamda;    // maximum Doppler frequency
   Np = 2*(2*N_o+1);

   switch(type)
   {
   	case 1:
   		alpha = 0.0;
         T_c2 = (long double)N_o;
         T_s2 = (long double)N_o + 1.0;
         break;
      case 2:
      	alpha = 0.0;
         T_c2 = (long double)N_o + 1.0;
         T_s2 = (long double)N_o;
         break;
      case 3:
      	alpha = Pi/4.0;
         T_c2 = (long double)N_o + 0.5;
         T_s2 = (long double)N_o + 0.5;
         break;
      default:
      	printf("\nInvalid type selection for Jake's fading channel model.\n");
         break;
   }

   if(v == 0.0)
   {
   	*(fade+0) = 1.0;
      *(fade+1) = 0.0;
   }
   else
   {
   	*(fade+0) = sqrtl(1.0/T_c2)*cosl(alpha)*cosl(w_m*t);
      *(fade+1) = sqrtl(1.0/T_s2)*sinl(alpha)*cosl(w_m*t);

      for(n = 1; n <= N_o; n++)
      {
      	switch(type)
         {
         	case 1:
            	beta_n = (long double)n*Pi/((long double)N_o+1.0);
               break;
            case 2:
            	beta_n = (long double)n*Pi/(long double)N_o;
               break;
            case 3:
            	beta_n = (long double)n*Pi/(long double)N_o;
               break;
         	default:
            	break;
         }
         w_n = w_m*cosl(2.0*Pi*(long double)n/(long double)Np);
//            theta_n = 2.0*Pi*((long double)rand()/(long double)RAND_MAX);  // random phase
			theta_n = 0.0;
         *(fade+0) += sqrtl(2.0/T_c2)*cosl(beta_n)*cosl(w_n*t+theta_n);
         *(fade+1) += sqrtl(2.0/T_s2)*sinl(beta_n)*cosl(w_n*t+theta_n);
		}
	}
}

void Multipath_fade_pattern(long double ****fade, int path_num)
{
	int i, j, l;
	long double gain[2];

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

void dfour1(long double data[],unsigned long nn,int isign)
//long double data[];
//unsigned long nn;
//int isign;
{
	unsigned long n,mmax,m,j,istep,i;
   long double wtemp,wr,wpr,wpi,wi,theta;
   long double tempr,tempi;

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
      wtemp=sinl(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sinl(theta);
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
