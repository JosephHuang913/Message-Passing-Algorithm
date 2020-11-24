/*======================================================*/
/* Author: Chao-wang Huang                                                                  */
/* Date: Wednesday, September 10, 2008                                             */
/* A Bit-based Message Passing Algorithm for ICI Channel is simulated */
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
void JakesFading(double, double, double, int, double *);
void Multipath_fade_pattern(double *****, int);
void dfour1(double *, unsigned long, int);

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 		3.14159265358979
#define num_packet 500						// number of packets simulated
const int Mc = 1;									// modulation order (1 for BPSK)
const int T = 2;										// Number of transmit antenna
const int R = 2;										// Number of receive antenna
const int Iteration = 3;
const int N_fft =  1024;
const int CP_length = N_fft/8;				// Length of Cyclic Prefix
const double sample_rate = 11.2e6;		/* Sampling Frequency */
/* Power Delay Profile of ITU Vehicular A Channel */
const int path_num = 6;
const int path_delay[6] = {0,3,8,12,19,28};	// Delay (sample) {0,310,710,1090,1730,2510} (ns)
const double path_gain[6] = {0.485,0.3853,0.0611,0.0485,0.0153,0.0049};	// Power Gain (linear scale)
const int delay_spread = 28;				// Maximum path delay (samples)
const double vc = 350.0;						/* speed of vehicle in km/hr */
const double C = 3.0e8;						/* light speed */
const double fc = 2.5e9;						/* carrier frequency */
const double OFDM_sym_duration = 91.43e-6;		/* OFDMA Symbol Duration without Cyclic Prefix */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
const double fdxT = Doppler * OFDM_sym_duration;
const int D = 2;										// Total ICI Subcarriers = 2D
const int ICI = 2*D+1;
const int B = Mc*T*ICI;
const int X = (int)pow(2.0,B-1.0);
double ***ts;
int *delay;

int main(void)
{
	time_t  ti, start, end;
	int d, i, j, k, l, m, n, p, r, s, t, x, z, *data, *data_bit, ***Ak, *err_count;
	double snr, Eb_No, noise_pwr, noise[2], Y[2], ***Yk, *err_rate, Ps;
	double P, sum, RQ_tmp[2], pdf, ******R_pn, ******Q_np, **sym_I, **sym_Q;
	double ****Out, /*fft,*/ *gain, *****fade, /*****fade_path,*/ ***F_l, ****H_I, ****H_Q;
	//double ****ISI_I, ****ISI_Q, **TX_signal_I, **TX_signal_Q, **RX_signal_I, **RX_signal_Q;
	FILE *ber, *records;

	start = time(NULL);
	printf("BER Performance of Bit-based Message Passing in MIMO Channel\n");
	printf("Number of transmit antenna is %d\n", T);
	printf("Number of receive antenna is %d\n", R);
	printf("Maximum Doppler frequency: %f\n", Doppler);
	printf("Channel fading rate: %f\n", fdxT);
	printf("Maximum number of bits of simulation = %d\n\n", Mc*T*num_packet*N_fft);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Record_Bit-based_Message_Passing.log", "a");
	fprintf(records, "BER Performance of Bit-based Message Passing in MIMO Channel\n");
	fprintf(records, "Number of transmit antenna is %d\n", T);
	fprintf(records, "Number of receive antenna is %d\n", R);
	fprintf(records, "Maximum Doppler frequency: %f\n", Doppler);
	fprintf(records, "Channel fading rate: %f\n", fdxT);
	fprintf(records, "Maximum number of bits of simulation = %d\n\n", Mc*T*num_packet*N_fft);
	fprintf(records, "Eb/No     BER\n");
	fflush(records);
	
	data_bit = new int[N_fft*Mc*T];
	err_count = new int[Iteration];
	err_rate = new double[Iteration];
	ts = new double**[R];
	for(i=0; i<R; i++)
		ts[i] = new double*[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			ts[i][j] = new double[path_num];
	gain = new double[path_num];
	delay = new int[path_num];
/*
	ISI_I = new double***[R];
	for(i=0; i<R; i++)
		ISI_I[i] = new double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			ISI_I[i][j] = new double*[path_num];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				ISI_I[i][j][l] = new double[delay_spread];
	ISI_Q = new double***[R];
	for(i=0; i<R; i++)
		ISI_Q[i] = new double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			ISI_Q[i][j] = new double*[path_num];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				ISI_Q[i][j][l] = new double[delay_spread];
	fft = new double[2*N_fft+1];
*/
	sym_I = new double*[T];
	for(i=0; i<T; i++)
		sym_I[i] = new double[N_fft];
	sym_Q = new double*[T];
	for(i=0; i<T; i++)
		sym_Q[i] = new double[N_fft];
/*
	TX_signal_I = new double*[T];
	for(i=0; i<T; i++)
		TX_signal_I[i] = new double[N_fft+CP_length];
	TX_signal_Q = new double*[T];
	for(i=0; i<T; i++)
		TX_signal_Q[i] = new double[N_fft+CP_length];
*/
	fade = new double****[R];
	for(i=0; i<R; i++)
		fade[i] = new double***[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			fade[i][j] = new double**[path_num];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				fade[i][j][l] = new double*[N_fft+CP_length];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				for(t=0; t<N_fft+CP_length; t++)
					fade[i][j][l][t] = new double[2];
/*
	fade_path = new double****[R];
	for(i=0; i<R; i++)
		fade_path[i] = new double***[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			fade_path[i][j] = new double**[path_num];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				fade_path[i][j][l] = new double*[N_fft+CP_length+delay_spread];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				for(t=0; t<N_fft+CP_length+delay_spread; t++)
					fade_path[i][j][l][t] = new double[2];
	RX_signal_I = new double*[R];
	for(i=0; i<R; i++)
		RX_signal_I[i] = new double[N_fft+CP_length];
	RX_signal_Q = new double*[R];
	for(i=0; i<R; i++)
		RX_signal_Q[i] = new double[N_fft+CP_length];
*/
	Yk = new double**[R];
	for(i=0; i<R; i++)
   		Yk[i] = new double*[N_fft];
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			Yk[i][l] = new double[2];
	F_l = new double**[path_num];
	for(i=0; i<path_num; i++)
		F_l[i] = new double*[N_fft];
	for(i=0; i<path_num; i++)
		for(l=0; l<N_fft; l++)
			F_l[i][l] = new double[2];
	H_I = new double***[R];
	for(i=0; i<R; i++)
  		H_I[i] = new double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_I[i][j] = new double*[N_fft];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				H_I[i][j][l] = new double[N_fft];
	H_Q = new double***[R];
	for(i=0; i<R; i++)
   		H_Q[i] = new double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_Q[i][j] = new double*[N_fft];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				H_Q[i][j][l] = new double[N_fft];
	Q_np = new double*****[T];
	for(i=0; i<T; i++)
   		Q_np[i] = new double****[R];
	for(i=0; i<T; i++)
   		for(j=0; j<R; j++)
   			Q_np[i][j] = new double***[N_fft];
	for(i=0; i<T; i++)
   		for(j=0; j<R; j++)
			for(t=0; t<N_fft; t++)
   				Q_np[i][j][t] = new double**[N_fft];
	for(i=0; i<T; i++)
   		for(j=0; j<R; j++)
			for(t=0; t<N_fft; t++)
				for(d=0; d<N_fft; d++)
   					Q_np[i][j][t][d] = new double*[Mc];
	for(i=0; i<T; i++)
   		for(j=0; j<R; j++)
			for(t=0; t<N_fft; t++)
				for(d=0; d<N_fft; d++)
					for(m=0; m<Mc; m++)
   						Q_np[i][j][t][d][m] = new double[2];
	R_pn = new double*****[R];
	for(i=0; i<R; i++)
   		R_pn[i] = new double****[T];
	for(i=0; i<R; i++)
   		for(j=0; j<T; j++)
   			R_pn[i][j] = new double***[N_fft];
	for(i=0; i<R; i++)
   		for(j=0; j<T; j++)
			for(t=0; t<N_fft; t++)
   				R_pn[i][j][t] = new double**[N_fft];
	for(i=0; i<R; i++)
   		for(j=0; j<T; j++)
			for(t=0; t<N_fft; t++)
				for(d=0; d<N_fft; d++)
   					R_pn[i][j][t][d] = new double*[Mc];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(t=0; t<N_fft; t++)
				for(d=0; d<N_fft; d++)
					for(m=0; m<Mc; m++)
   						R_pn[i][j][t][d][m] = new double[2];
	data = new int[Mc*T*N_fft];
	Ak = new int**[T];
	for(i=0; i<T; i++)
	Ak[i] = new int*[N_fft];
	for(i=0; i<T; i++)
		for(t=0; t<N_fft; t++)
			Ak[i][t] = new int[Mc];
	Out = new double***[T];
	for(i=0; i<T; i++)
   		Out[i] = new double**[N_fft];
	for(i=0; i<T; i++)
		for(t=0; t<N_fft; t++)
   			Out[i][t] = new double*[Mc];
	for(i=0; i<T; i++)
		for(t=0; t<N_fft; t++)
			for(m=0; m<Mc; m++)
   				Out[i][t][m] = new double[2];
							
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
	fopen_s(&ber, "BER_Bit-based_Message_Passing.log", "w");
	for(snr=30; snr<=30; snr+=5)
	{
   		for(s=0; s<Iteration; s++)
   			err_count[s] = 0;

   		// noise power calculation
		Eb_No = (double)snr;
		Ps = 1 * T;
		//noise_pwr = 0.5*Ps*R/((float)N_fft*T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Time Domain)
		//noise_pwr = 0.5*Ps*R/((float)T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Freq. Domain)
		noise_pwr = 1.0*Ps*R/((float)T*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption (Freq. Domain)
		//noise_pwr = 0.5*Ps*R/((k/(float)n)*T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Freq. Domain)
		//noise_pwr = 0.5*Ps*R/(N_fft*(k/(float)n)*T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Time Domain)
		printf("Eb_No = %f\n", Eb_No);

		// Initial time of Rayleigh fading pattern
		for(i=0; i<R; i++)
			for(j=0; j<T; j++)
				for(l=0; l<path_num; l++)
					ts[i][j][l] = 1000.0 + i*1000.0 + j*500.0 + l*50.0;
/*
		// Initialize ISI buffer
		for(i=0; i<R; i++)
			for(j=0; j<T; j++)
				for(l=1; l<path_num; l++)		// path index
					for(t=0; t<delay[l]; t++)	// path delay
					{
						ISI_I[i][j][l][t] = 0.0;
						ISI_Q[i][j][l][t] = 0.0;
					}
*/
		p = 0;
		do
		{
			cout << "Packet: " << p << endl; 
			for(l=0; l<N_fft; l++)
				for(i=0; i<T; i++)
					for(m=0; m<Mc; m++)
						// Generate random information bit stream
						if(rand()/(float)RAND_MAX>=0.5)
							data_bit[l*Mc*T+i*Mc+m] = 1;
						else
							data_bit[l*Mc*T+i*Mc+m] = 0;

			for(i=0; i<T; i++)
				for(l=0; l<N_fft; l++)
				{
					// BPSK Mapping
					sym_I[i][l] = (2*data_bit[l*Mc*T+i*Mc]-1);
         			sym_Q[i][l] = 0.0;
					// QPSK Mapping
					//sym_I[i][l] = (2*data_bit[l*Mc*T+i*Mc]-1)/sqrt(2.0);
         			//sym_Q[i][l] = (2*data_bit[l*Mc*T+i*Mc+1]-1)/sqrt(2.0);
				}
/*
			// IFFT
			for(i=0; i<T; i++)
			{
				for(l=0; l<N_fft; l++)
				{
					fft[2*l+1] = sym_I[i][l];
					fft[2*l+2] = sym_Q[i][l];
				}

				dfour1(&fft[0], N_fft, 1);  // IFFT for transmitted signal must be multiplied by N_fft^(-1)

				for(l=0; l<N_fft; l++)
				{
					sym_I[i][l] = fft[2*l+1] / (double)N_fft;
					sym_Q[i][l] = fft[2*l+2] / (double)N_fft;
				}
			}

			// Cyclic Prefix
			for(i=0; i<T; i++)
				for(j=0; j<CP_length; j++)
				{
					TX_signal_I[i][j] = sym_I[i][j+N_fft-CP_length];
					TX_signal_Q[i][j] = sym_Q[i][j+N_fft-CP_length];
				}
			
			for(i=0; i<T; i++)
				for(j=CP_length; j<N_fft+CP_length; j++)
				{
					TX_signal_I[i][j] = sym_I[i][j-CP_length];
					TX_signal_Q[i][j] = sym_Q[i][j-CP_length];
				}
*/
			// Multi-path Rayleigh fading pattern
			Multipath_fade_pattern(fade, path_num);
/*
			// Multipath Signal
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
					for(l=0; l<path_num; l++)
						for(t=delay[l]; t<N_fft+CP_length+delay[l]; t++)
						{
							fade_path[i][j][l][t][0] = gain[l] * (TX_signal_I[j][t-delay[l]] * fade[i][j][l][t][0] - TX_signal_Q[j][t-delay[l]] * fade[i][j][l][t][1]);
							fade_path[i][j][l][t][1] = gain[l] * (TX_signal_I[j][t-delay[l]] * fade[i][j][l][t][1] + TX_signal_Q[j][t-delay[l]] * fade[i][j][l][t][0]);
						}

			// ISI Records
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
					for(l=1; l<path_num; l++)
						for(t=0; t<delay[l]; t++)
						{
							fade_path[i][j][l][t][0] = ISI_I[i][j][l][t];
							fade_path[i][j][l][t][1] = ISI_Q[i][j][l][t];
							ISI_I[i][j][l][t] = fade_path[i][j][l][t+N_fft+CP_length][0];
							ISI_Q[i][j][l][t] = fade_path[i][j][l][t+N_fft+CP_length][1];
						}

			// Multipath Channel
			for(i=0; i<R; i++)
				for(t=0; t<N_fft+CP_length; t++)
				{
					RX_signal_I[i][t] = 0.0;
					RX_signal_Q[i][t] = 0.0;
					for(j=0; j<T; j++)
						for(l=0; l<path_num; l++)
						{
							RX_signal_I[i][t] += fade_path[i][j][l][t][0];
							RX_signal_Q[i][t] += fade_path[i][j][l][t][1];
						}
				
					// AWGN noise
					AWGN_noise(0.0, noise_pwr, &noise[0]);
      	   			RX_signal_I[i][t] += noise[0];
		         	RX_signal_Q[i][t] += noise[1];
				}

			// Remove Cyclic Prefix
			for(i=0; i<R; i++)
				for(t=0; t<N_fft; t++)
				{
					Yk[i][t][0] = RX_signal_I[i][t+CP_length];
					Yk[i][t][1] = RX_signal_Q[i][t+CP_length];
				}

			// FFT
			for(i=0; i<R; i++)
			{
				for(t=0; t<N_fft; t++)
				{
					fft[2*t+1] = Yk[i][t][0];
					fft[2*t+2] = Yk[i][t][1];
				}

				dfour1(&fft[0], N_fft, -1);

				for(t=0; t<N_fft; t++)
				{
					Yk[i][t][0] = fft[2*t+1];
					Yk[i][t][1] = fft[2*t+2];
				}
			}
*/
			// Channel Frequency Response
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
				{
					for(l=0; l<path_num; l++)
						for(k=0; k<N_fft; k++)
						{
							F_l[l][k][0] = F_l[l][k][1] = 0.0;
							for(t=0; t<N_fft; t++)
							{
								F_l[l][k][0] += gain[l]*fade[i][j][l][t+CP_length-delay[l]][0] * cos(2*Pi*t*k/(float)N_fft) + gain[l]*fade[i][j][l][t+CP_length-delay[l]][1] * sin(2*Pi*t*k/(float)N_fft);
								F_l[l][k][1] += gain[l]*fade[i][j][l][t+CP_length-delay[l]][1] * cos(2*Pi*t*k/(float)N_fft) - gain[l]*fade[i][j][l][t+CP_length-delay[l]][0] * sin(2*Pi*t*k/(float)N_fft);
							}
						}

					// Channel Frequency Response & ICI Channel Coefficients
					for(t=0; t<N_fft; t++)
						for(d=0; d<N_fft; d++)
						{
							H_I[i][j][t][d] = H_Q[i][j][t][d] = 0.0;
							for(l=0; l<path_num; l++)
							{
								H_I[i][j][t][d] += F_l[l][d][0] * cos(2*Pi*delay[l]*(t-d)/(float)N_fft) + F_l[l][d][1] * sin(2*Pi*delay[l]*(t-d)/(float)N_fft);
								H_Q[i][j][t][d] += F_l[l][d][1] * cos(2*Pi*delay[l]*(t-d)/(float)N_fft) - F_l[l][d][0] * sin(2*Pi*delay[l]*(t-d)/(float)N_fft);
							}
		
							H_I[i][j][t][d] /= (float)N_fft;
							H_Q[i][j][t][d] /= (float)N_fft;
						}
				}

			// Freq. Domain Transmission
			for(i=0; i<R; i++)
				for(t=0; t<N_fft; t++)
					{
						Yk[i][t][0] = Yk[i][t][1] = 0.0;
						for(d=0; d<N_fft; d++)
							for(j=0; j<T; j++)
							{
								Yk[i][t][0] += (sym_I[j][(N_fft+t-d)%N_fft] * H_I[i][j][t][d] - sym_Q[j][(N_fft+t-d)%N_fft] * H_Q[i][j][t][d]);
								Yk[i][t][1] += (sym_I[j][(N_fft+t-d)%N_fft] * H_Q[i][j][t][d] + sym_Q[j][(N_fft+t-d)%N_fft] * H_I[i][j][t][d]);
							}

						// AWGN noise
						AWGN_noise(0.0, noise_pwr, &noise[0]);
		  	   			Yk[i][t][0] += noise[0];
				       	Yk[i][t][1] += noise[1];
					}

/*==================================================*/
/*  B i t - b a s e d    M e s s a g e    P a s s i n g    A l g o r i t h m  */
/*==================================================*/

			// Initilization of extrinsic(intrinsic) probability
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
					for(m=0; m<Mc; m++)
						for(t=0; t<N_fft; t++)
							for(d=0; d<N_fft; d++)
      							for(k=0; k<2; k++)
								{
   									Q_np[j][i][t][d][m][k] = 0.5;
									R_pn[i][j][t][d][m][k] = 0.5;
								}

			for(s=0; s<Iteration; s++)
			{
				// Extrinsic information of function nodes
            	for(i=0; i<R; i++)
					for(t=0; t<N_fft; t++)
						for(d=0; d<N_fft; d++)
							if(d <= D || d >= N_fft-D)
								for(j=0; j<T; j++)
									for(m=0; m<Mc; m++)
									{
										for(k=0; k<2; k++)		// Probability of 1 and 0
										{
											RQ_tmp[k] = 0.0;
											for(x=0; x<X; x++)
   											{
												z = 0;
												data[((N_fft+t-d)%N_fft)*Mc*T+j*Mc+m] = k;
												for(l=0; l<N_fft; l++)
													if(l <= D || l >= N_fft-D)
														for(r=0; r<T; r++)
															for(n=0; n<Mc; n++)
		           				     							if(((N_fft+t-l)%N_fft)*Mc*T+r*Mc+n != ((N_fft+t-d)%N_fft)*Mc*T+j*Mc+m)
					   											{
							  	      								data[((N_fft+t-l)%N_fft)*Mc*T+r*Mc+n] = (x>>z) & 1;
				                     								z++;
						                						}

												for(l=0; l<N_fft; l++)
													if(l <= D || l >= N_fft-D)
														for(r=0; r<T; r++)
														{
															// BPSK Mapping
         	            									sym_I[r][(N_fft+t-l)%N_fft] = (2*data[((N_fft+t-l)%N_fft)*Mc*T+r*Mc]-1);
				       										sym_Q[r][(N_fft+t-l)%N_fft] = 0.0;
														}
      												
												Y[0] = Y[1] = 0.0;
												for(l=0; l<N_fft; l++)
													if(l <= D || l >= N_fft-D)
														for(r=0; r<T; r++)
														{
															Y[0] += (sym_I[r][(N_fft+t-l)%N_fft] * H_I[i][r][t][l] - sym_Q[r][(N_fft+t-l)%N_fft] * H_Q[i][r][t][l]);
															Y[1] += (sym_I[r][(N_fft+t-l)%N_fft] * H_Q[i][r][t][l] + sym_Q[r][(N_fft+t-l)%N_fft] * H_I[i][r][t][l]);
														}

												P = 1.0;
												if(s != 0)
						                   			for(l=0; l<N_fft; l++)
														if(l <= D || l >= N_fft-D)
															for(r=0; r<T; r++)
																for(n=0; n<Mc; n++)
										     						if(((N_fft+t-l)%N_fft)*Mc*T+r*Mc+n != ((N_fft+t-d)%N_fft)*Mc*T+j*Mc+m)
																		P *= Q_np[r][i][(N_fft+t-l)%N_fft][t][n][data[((N_fft+t-l)%N_fft)*Mc*T+r*Mc+n]];

                     							sum = 0.0;
												for(z=0; z<2; z++)
                        							sum += pow(Yk[i][t][z]-Y[z],2.0);

												pdf = exp(-sum/noise_pwr)/(Pi*noise_pwr);
		                     					RQ_tmp[k] += pdf * P;
                  							}
										}

										R_pn[i][j][t][((N_fft+t-d)%N_fft)][m][1] = RQ_tmp[1] / (RQ_tmp[0] + RQ_tmp[1]);
										R_pn[i][j][t][((N_fft+t-d)%N_fft)][m][0] = RQ_tmp[0] / (RQ_tmp[0] + RQ_tmp[1]);
									}
            	
				// Extrinsic information of bit nodes
				for(t=0; t<N_fft; t++)
					for(j=0; j<T; j++)
						for(m=0; m<Mc; m++)
							for(d=0; d<N_fft; d++)
								if(d <= D || d >= N_fft-D)
									for(i=0; i<R; i++)
									{
										RQ_tmp[0] = RQ_tmp[1] = 1.0;
		           						for(l=0; l<N_fft; l++)
											if(l <= D || l >= N_fft-D)
												for(n=0; n<R; n++)
           											if(l != d || n != i)
													{
			                 							RQ_tmp[1] *= R_pn[n][j][((N_fft+t-l)%N_fft)][t][m][1];
														RQ_tmp[0] *= R_pn[n][j][((N_fft+t-l)%N_fft)][t][m][0];
													}
												
										Q_np[j][i][t][((N_fft+t-d)%N_fft)][m][1] = RQ_tmp[1] / (RQ_tmp[1] + RQ_tmp[0]);
										Q_np[j][i][t][((N_fft+t-d)%N_fft)][m][0] = RQ_tmp[0] / (RQ_tmp[1] + RQ_tmp[0]);
									}

				// Bit Nodes Decision Output
				for(t=0; t<N_fft; t++)
					for(j=0; j<T; j++)
						for(m=0; m<Mc; m++)
						{
							RQ_tmp[0] = RQ_tmp[1] = 1.0;
							for(d=0; d<N_fft; d++)
								if(d <= D || d>= N_fft-D)
									for(i=0; i<R; i++)
									{
										RQ_tmp[1] *= R_pn[i][j][((N_fft+t-d)%N_fft)][t][m][1];
										RQ_tmp[0] *= R_pn[i][j][((N_fft+t-d)%N_fft)][t][m][0];
									}

							Out[j][t][m][1] = RQ_tmp[1] / (RQ_tmp[1] + RQ_tmp[0]);
							Out[j][t][m][0] = RQ_tmp[0] / (RQ_tmp[1] + RQ_tmp[0]);

							if(Out[j][t][m][1] >= Out[j][t][m][0])
            					Ak[j][t][m] = 1;
       						else
		            			Ak[j][t][m] = 0;

            				err_count[s] += Error_count(data_bit[t*Mc*T+j*Mc+m], Ak[j][t][m]);
						}
			}
	
			p++;
		} while(err_count[Iteration-1] <= 1000 && p < num_packet);

		// Statistics and records
		cout << "Error Rate = ";
		for(s=0; s<Iteration; s++)
		{
      		err_rate[s] = err_count[s] / (double)(Mc*T*p*N_fft);
      		printf("%e, ", err_rate[s]);
		}
		cout << endl;

		fprintf(ber, "%f ", Eb_No);
		fprintf(records, "%f ", Eb_No);
		for(s=0; s<Iteration; s++)
		{
			fprintf(ber, "%e ", err_rate[s]);
			fprintf(records, "%e %d ", err_rate[s], err_count[s]);
		}
		fprintf(ber, "\n");
		fprintf(records, "\n");
		fflush(records);
		fflush(ber);
	}

	delete data_bit;
	delete err_count;
	delete err_rate;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete ts[i][j];
	for(i=0; i<R; i++)
		delete ts[i];
	delete ts;
	delete gain;
	delete delay;
/*
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				delete ISI_I[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete ISI_I[i][j];
	for(i=0; i<R; i++)
		delete ISI_I[i];
	delete ISI_I;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				delete ISI_Q[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete ISI_Q[i][j];
	for(i=0; i<R; i++)
		delete ISI_Q[i];
	delete ISI_Q;
	delete fft;
*/
	for(i=0; i<T; i++)
		delete sym_I[i];
	delete sym_I;
	for(i=0; i<T; i++)
		delete sym_Q[i];
	delete sym_Q;
/*
	for(i=0; i<T; i++)
		delete TX_signal_I[i];
	delete TX_signal_I;
	for(i=0; i<T; i++)
		delete TX_signal_Q[i];
	delete TX_signal_Q;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				for(t=0; t<N_fft+CP_length+delay_spread; t++)
					delete fade_path[i][j][l][t];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				delete fade_path[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete fade_path[i][j];
	for(i=0; i<R; i++)
		delete fade_path[i];
	delete fade_path;
	for(i=0; i<R; i++)
		delete RX_signal_I[i];
	delete RX_signal_I;
	for(i=0; i<R; i++)
		delete RX_signal_Q[i];
	delete RX_signal_Q;
*/
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				for(t=0; t<N_fft+CP_length; t++)
					delete fade[i][j][l][t];
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
	for(i=0; i<path_num; i++)
		for(l=0; l<N_fft; l++)
			delete F_l[i][l];
	for(i=0; i<path_num; i++)
		delete F_l[i];
	delete F_l;
	delete data;
	for(i=0; i<T; i++)
		for(t=0; t<N_fft; t++)
			delete Ak[i][t];
	for(i=0; i<T; i++)
		delete Ak[i];
	delete Ak;
	for(i=0; i<T; i++)
		for(t=0; t<N_fft; t++)
			for(m=0; m<Mc; m++)
   				delete Out[i][t][m];
	for(i=0; i<T; i++)
		for(t=0; t<N_fft; t++)
   			delete Out[i][t];
	for(i=0; i<T; i++)
		delete Out[i];
	delete Out;
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			delete Yk[i][l];
	for(i=0; i<R; i++)
   		delete Yk[i];
	delete Yk;
	for(i=0; i<T; i++)
   		for(j=0; j<R; j++)
			for(t=0; t<N_fft; t++)
				for(d=0; d<N_fft; d++)
					for(m=0; m<Mc; m++)
   						delete Q_np[i][j][t][d][m];
	for(i=0; i<T; i++)
   		for(j=0; j<R; j++)
			for(t=0; t<N_fft; t++)
				for(d=0; d<N_fft; d++)
   					delete Q_np[i][j][t][d];
	for(i=0; i<T; i++)
   		for(j=0; j<R; j++)
			for(t=0; t<N_fft; t++)
   				delete Q_np[i][j][t];
	for(i=0; i<T; i++)
   		for(j=0; j<R; j++)
   			delete Q_np[i][j];
	for(i=0; i<T; i++)
   		delete Q_np[i];
	delete Q_np;
	for(i=0; i<R; i++)
   		for(j=0; j<T; j++)
			for(t=0; t<N_fft; t++)
				for(d=0; d<N_fft; d++)
					for(m=0; m<Mc; m++)
   						delete R_pn[i][j][t][d][m];
	for(i=0; i<R; i++)
   		for(j=0; j<T; j++)
			for(t=0; t<N_fft; t++)
				for(d=0; d<N_fft; d++)
					delete R_pn[i][j][t][d];
	for(i=0; i<R; i++)
   		for(j=0; j<T; j++)
			for(t=0; t<N_fft; t++)
   				delete R_pn[i][j][t];
	for(i=0; i<R; i++)
 		for(j=0; j<T; j++)
		   	delete R_pn[i][j];
	for(i=0; i<R; i++)
   		delete R_pn[i];
	delete R_pn;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				delete H_I[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete H_I[i][j];
	for(i=0; i<R; i++)
   		delete H_I[i];
	delete H_I;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				delete H_Q[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete H_Q[i][j];
	for(i=0; i<R; i++)
		   delete H_Q[i];
	delete H_Q;

	end = time(NULL);
	printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
	fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
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
   int n, Np, N_o = 32;
   double lamda, w_m, beta_n, w_n, alpha, T_c2, T_s2, theta_n;

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

   if(v == 0.0)
   {
   	*(fade+0) = 1.0;
      *(fade+1) = 0.0;
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
//            theta_n = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
			theta_n = 0.0;
         *(fade+0) += sqrt(2.0/T_c2)*cos(beta_n)*cos(w_n*t+theta_n);
         *(fade+1) += sqrt(2.0/T_s2)*sin(beta_n)*cos(w_n*t+theta_n);
		}
	}
}

void Multipath_fade_pattern(double *****fade, int path_num)
{
	int i, j, l, t;
	double gain[2];

	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				for(t=0; t<N_fft+CP_length; t++)
				{
					JakesFading(fc, vc*1000/3600.0, ts[i][j][l], 2, &gain[0]);
					ts[i][j][l] += OFDM_sym_duration / (float)N_fft;
					fade[i][j][l][t][0] = gain[0];
					fade[i][j][l][t][1] = gain[1];
				}
}

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
