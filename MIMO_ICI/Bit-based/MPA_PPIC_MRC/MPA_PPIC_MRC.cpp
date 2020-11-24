/*========================================================*/
/* Author: Chao-wang Huang                                                                       */
/* Date: Wednesday, March 04, 2008                                                         */
/* A Bit-based Message Passing Algorithm for MIMO Channel is simulated */
/* 2 by 2 QPSK                                                                                           */
/* MIMO Detector: MPA with Sum-Product Rule                                          */
/* ICI Canceller: Progressive PIC (PPIC)                                                    */
/* Maximum Ratio Combining                                                                      */
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
void Multipath_fade_pattern(double *****, int);

#define Pi 	3.14159265358979
#define num_packet 125					// number of packets simulated
const int Mc = 4;								// modulation order (4 for 16QAM)
//const int Mc = 2;								// modulation order (2 for QPSK)
const int T = 2;									// Number of transmit antenna
const int R = 2;									// Number of receive antenna
const int Iterate2 = 6;
const int Iterate1 = 3;
const int N_fft =  1024;
const int CP_length = N_fft/8;				// Length of Cyclic Prefix
const double sample_rate = 11.2e6;		/* Sampling Frequency */
/* Power Delay Profile of ITU Vehicular A Channel */
const int path_num = 6;
const int path_delay[6] = {0,3,8,12,19,28};	// Delay (sample) {0,310,710,1090,1730,2510} (ns)
const double path_gain[6] = {0.485,0.3853,0.0611,0.0485,0.0153,0.0049};	// Power Gain (linear scale)
const int delay_spread = 28;				// Maximum path delay (samples)
const double vc = 500.0;						/* speed of vehicle in km/hr */
const double C = 3.0e8;						/* light speed */
const double fc = 2.5e9;						/* carrier frequency */
const double OFDM_sym_duration = 91.43e-6;		/* OFDMA Symbol Duration without Cyclic Prefix */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
const double fdxT = Doppler * OFDM_sym_duration;
double ***ts;

int main(void)
{
	time_t  ti, start, end;
	int i, j, k, l, m, n, p, r, s1, s2, s3, x, z, t, d, *data, **data_bit, *Ak, **err_count, *delay;
	double snr, Eb_No, noise_pwr, noise[2], Y[2], ***Yk, **err_rate, Ps, ***Out;
	double P, sum, RQ_tmp[2], pdf, ****R_pn, ****Q_np, **sym_I, **sym_Q, *I, *Q;
	double *gain, *****fade, ***F_l, ****H_I, ****H_Q, ***Est_sym_I, ***Est_sym_Q;
	double **Re_ICI_I, **Re_ICI_Q, ***Re_signal, **Sk;
	FILE *ber, *records;

	start = time(NULL);
	printf("BER Performance of Bit-based Message Passing in MIMO Channel\n");
	printf("Receiver Channel Selection\n");
	printf("Number of transmit antenna is %d\n", T);
	printf("Number of receive antenna is %d\n", R);
	printf("Maximum Doppler frequency: %f\n", Doppler);
	printf("Channel fading rate: %f\n", fdxT);
	printf("Maximum number of bits of simulation = %d\n\n", Mc*T*num_packet*N_fft);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Record_Bit-based_Message_Passing.log", "a");
	fprintf(records, "BER Performance of Bit-based Message Passing in MIMO Channel\n");
	fprintf(records, "Receiver Channel Selection\n");
	fprintf(records, "Number of transmit antenna is %d\n", T);
	fprintf(records, "Number of receive antenna is %d\n", R);
	fprintf(records, "Maximum Doppler frequency: %f\n", Doppler);
	fprintf(records, "Channel fading rate: %f\n", fdxT);
	fprintf(records, "Maximum number of bits of simulation = %d\n\n", Mc*T*num_packet*N_fft);
	fprintf(records, "Eb/No     BER\n");
	fflush(records);

	data_bit = new int*[N_fft];
	for(i=0; i<N_fft; i++)
		data_bit[i] = new int[Mc*T];
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
	Ak = new int[Mc*T];
	Sk = new double*[N_fft];
	for(i=0; i<N_fft; i++)
		Sk[i] = new double[Mc*T];
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
	Yk = new double**[R];
	for(i=0; i<R; i++)
   		Yk[i] = new double*[N_fft];
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			Yk[i][l] = new double[2];
	Re_signal = new double**[R];
	for(i=0; i<R; i++)
   		Re_signal[i] = new double*[N_fft];
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			Re_signal[i][l] = new double[2];
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
	R_pn = new double***[N_fft];
	for(k=0; k<N_fft; k++)
		R_pn[k] = new double**[R];
	for(k=0; k<N_fft; k++)
		for(i=0; i<R; i++)
   			R_pn[k][i] = new double*[Mc*T];
	for(k=0; k<N_fft; k++)
		for(i=0; i<R; i++)
   			for(j=0; j<Mc*T; j++)
   				R_pn[k][i][j] = new double[2];
	Q_np = new double***[N_fft];
	for(k=0; k<N_fft; k++)
		Q_np[k] = new double**[Mc*T];
	for(k=0; k<N_fft; k++)
		for(i=0; i<Mc*T; i++)
   			Q_np[k][i] = new double*[R];
	for(k=0; k<N_fft; k++)
		for(i=0; i<Mc*T; i++)
   			for(j=0; j<R; j++)
   				Q_np[k][i][j] = new double[2];
	Out = new double**[N_fft];
	for(k=0; k<N_fft; k++)
		Out[k] = new double*[Mc*T];
	for(k=0; k<N_fft; k++)
		for(i=0; i<Mc*T; i++)
   			Out[k][i] = new double[2];
	err_count = new int*[Iterate2];
	for(i=0; i<Iterate2; i++)
		err_count[i] = new int[Iterate1];	
	err_rate = new double*[Iterate2];
	for(i=0; i<Iterate2; i++)
		err_rate[i] = new double[Iterate1];
	Est_sym_I = new double**[Iterate2];
	for(i=0; i<Iterate2; i++)
		Est_sym_I[i] = new double*[T];
	for(i=0; i<Iterate2; i++)
		for(j=0; j<T; j++)
			Est_sym_I[i][j] = new double[N_fft];
	Est_sym_Q = new double**[Iterate2];
	for(i=0; i<Iterate2; i++)
		Est_sym_Q[i] = new double*[T];
	for(i=0; i<Iterate2; i++)
		for(j=0; j<T; j++)
			Est_sym_Q[i][j] = new double[N_fft];
	Re_ICI_I = new double*[R];
	for(i=0; i<R; i++)
		Re_ICI_I[i] = new double[N_fft];
	Re_ICI_Q = new double*[R];
	for(i=0; i<R; i++)
		Re_ICI_Q[i] = new double[N_fft];
	I = new double[T];
	Q = new double[T];

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
	fopen_s(&ber, "BER_Bit_based_Message_Passing.log", "w");
	for(snr=0; snr<=30; snr+=5)
	{
   		for(s2=0; s2<Iterate2; s2++)
			for(s1=0; s1<Iterate1; s1++)
   				err_count[s2][s1] = 0;

   		// noise power calculation
		Eb_No = (double)snr;
		Ps = 1 * T;
		//noise_pwr = 0.5*Ps*R/((float)N_fft*T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Time Domain)
		//noise_pwr = 0.5*Ps*R/((float)T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Freq. Domain)
		noise_pwr = 0.25*Ps*R/((float)T*pow(10.0, Eb_No/10.0));	// 16QAM, Nyquist filter assumption (Freq. Domain)
		printf("Eb_No = %f\n", Eb_No);

		// Initial time of Rayleigh fading pattern
		for(i=0; i<R; i++)
			for(j=0; j<T; j++)
				for(l=0; l<path_num; l++)
					ts[i][j][l] = 1000.0 + i*1000.0 + j*500.0 + l*50.0;

		p = 0;
		do
		{
			cout << "Packet: " << p << endl;
			for(l=0; l<N_fft; l++)
				for(i=0; i<T; i++)
					for(m=0; m<Mc; m++)
						// Generate random information bit stream
						if(rand()/(float)RAND_MAX>=0.5)
							data_bit[l][i*Mc+m] = 1;
						else
							data_bit[l][i*Mc+m] = 0;

			for(l=0; l<N_fft; l++)
   				for(i=0; i<T; i++)
	   			{
					// QPSK Mapping
					//sym_I[i][l] = (2*data_bit[l][i*Mc]-1)/sqrt(2.0);
         			//sym_Q[i][l] = (2*data_bit[l][i*Mc+1]-1)/sqrt(2.0);

					// 16-QAM Mapping	
					if(data_bit[l][Mc*i] == 0 && data_bit[l][Mc*i+2] == 0)
   		  				sym_I[i][l] = 1.0/sqrt(10.0);
					else if(data_bit[l][Mc*i] == 0 && data_bit[l][Mc*i+2] == 1)
   			  			sym_I[i][l] = 3.0/sqrt(10.0);
					else if(data_bit[l][Mc*i] == 1 && data_bit[l][Mc*i+2] == 1)
		  		 		sym_I[i][l] = -3.0/sqrt(10.0);
				  	else
  		   				sym_I[i][l] = -1.0/sqrt(10.0);

		   			if(data_bit[l][Mc*i+1] == 0 && data_bit[l][Mc*i+3] == 0)
  	   					sym_Q[i][l] = 1.0/sqrt(10.0);
					else if(data_bit[l][Mc*i+1] == 0 && data_bit[l][Mc*i+3] == 1)
   				  		sym_Q[i][l] = 3.0/sqrt(10.0);
				     else if(data_bit[l][Mc*i+1] == 1 && data_bit[l][Mc*i+3] == 1)
	   			 		sym_Q[i][l] = -3.0/sqrt(10.0);
					 else
      		 			sym_Q[i][l] = -1.0/sqrt(10.0);
				}

			// Multi-path Rayleigh fading pattern
			Multipath_fade_pattern(fade, path_num);

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

/*==============================================*/
/*  Bit-based Message Passing Algorithm & Progressive PIC */
/*==============================================*/

			for(s1=0; s1<Iterate1; s1++)
				for(t=0; t<N_fft; t++)
				{
					// Extrinsic information of function nodes
		    		for(i=0; i<R; i++)
						for(j=0; j<T; j++)
							for(m=0; m<Mc; m++)
           					{
								for(k=0; k<2; k++)		// Probability of 1 and 0
								{
               						RQ_tmp[k] = 0.0;
		       						for(x=0; x<pow(2.0,Mc*T-1.0); x++)
	   								{
		      	            			l = 0;
			   							for(n=0; n<Mc*T; n++)
				   		       				if(n != j*Mc+m)
           									{
					  			  				data[n] = (x>>l) & 1;
						 						l++;
											}

										data[j*Mc+m] = k;
   										for(r=0; r<T; r++)
		  								{
											// QPSK Mapping
											//I[r] = (2*data[r*Mc]-1)/sqrt(2.0);
         									//Q[r] = (2*data[r*Mc+1]-1)/sqrt(2.0);

											// 16-QAM Mapping	
											if(data[Mc*r] == 0 && data[Mc*r+2] == 0)
						   		  				I[r] = 1.0/sqrt(10.0);
											else if(data[Mc*r] == 0 && data[Mc*r+2] == 1)
   			  									I[r] = 3.0/sqrt(10.0);
											else if(data[Mc*r] == 1 && data[Mc*r+2] == 1)
   		 										I[r] = -3.0/sqrt(10.0);
			  								else
  		   										I[r] = -1.0/sqrt(10.0);
	
				   							if(data[Mc*r+1] == 0 && data[Mc*r+3] == 0)
  			   									Q[r] = 1.0/sqrt(10.0);
											else if(data[Mc*r+1] == 0 && data[Mc*r+3] == 1)
   					  							Q[r] = 3.0/sqrt(10.0);
										     else if(data[Mc*r+1] == 1 && data[Mc*r+3] == 1)
								   		 		Q[r] = -3.0/sqrt(10.0);
										     else
      	 										Q[r] = -1.0/sqrt(10.0);
           								}
	
		                 				Y[0] = Y[1] = 0.0;
										for(r=0; r<T; r++)
										{
							   				Y[0] += (I[r] * H_I[i][r][t][0] - Q[r] * H_Q[i][r][t][0]);
											Y[1] += (I[r] * H_Q[i][r][t][0] + Q[r] * H_I[i][r][t][0]);
										}
									
					 					P = 1.0;
										if(s1 !=0)
	                   						for(r=0; r<T; r++)
												for(n=0; n<Mc; n++)
		             								if(r*Mc+n != j*Mc+m)
				       									P *= Q_np[t][r*Mc+n][i][data[r*Mc+n]];

										sum = 0.0;
										for(z=0; z<2; z++)
				               					sum += pow(Yk[i][t][z]-Y[z],2.0);

										pdf = exp(-sum/noise_pwr)/(Pi*noise_pwr);
					     				RQ_tmp[k] += pdf * P;
               						}
           						}

								R_pn[t][i][j*Mc+m][1] = RQ_tmp[1] / (RQ_tmp[0] + RQ_tmp[1]);
          						R_pn[t][i][j*Mc+m][0] = RQ_tmp[0] / (RQ_tmp[0] + RQ_tmp[1]);
							}

					// Extrinsic information of bit nodes
					for(j=0; j<T; j++)
						for(m=0; m<Mc; m++)
							for(i=0; i<R; i++)
							{
	           					RQ_tmp[0] = RQ_tmp[1] = 1.0;
			   					for(n=0; n<R; n++)
           							if(n != i)
									{
						 				RQ_tmp[1] *= R_pn[t][n][j*Mc+m][1];
										RQ_tmp[0] *= R_pn[t][n][j*Mc+m][0];
									}

								Q_np[t][j*Mc+m][i][1] = RQ_tmp[1] / (RQ_tmp[1] + RQ_tmp[0]);
								Q_np[t][j*Mc+m][i][0] = RQ_tmp[0] / (RQ_tmp[1] + RQ_tmp[0]);
		 					}

					// Decision Output
					for(j=0; j<T; j++)
						for(m=0; m<Mc; m++)
						{
       						RQ_tmp[0] = RQ_tmp[1] = 1.0;
							for(i=0; i<R; i++)
							{	
								RQ_tmp[1] *= R_pn[t][i][j*Mc+m][1];
								RQ_tmp[0] *= R_pn[t][i][j*Mc+m][0];
							}

							Out[t][j*Mc+m][1] = RQ_tmp[1] / (RQ_tmp[1] + RQ_tmp[0]);
							Out[t][j*Mc+m][0] = RQ_tmp[0] / (RQ_tmp[1] + RQ_tmp[0]);
						}

					for(j=0; j<Mc*T; j++)
					{
						// Hard Decision
				 		if(Out[t][j][1] >= Out[t][j][0])
							Ak[j] = 1;
       					else
            				Ak[j] = 0;

	           			err_count[0][s1] += Error_count(data_bit[t][j], Ak[j]);
					}
				}

				for(t=0; t<N_fft; t++)
				{
					for(j=0; j<Mc*T; j++)
						// Bit Soft Decision
	            		Sk[t][j] = tanh(0.5*log(Out[t][j][1]/Out[t][j][0]));

					for(j=0; j<T; j++)
					{
						// QPSK Soft Mapping
						//Est_sym_I[0][j][t] = Sk[t][Mc*j] / sqrt(2.0);
						//Est_sym_Q[0][j][t] = Sk[t][Mc*j+1] / sqrt(2.0);

						// 16-QAM Soft Mapping
						Est_sym_I[0][j][t] = (-1.0/sqrt(10.0)) * Sk[t][Mc*j] * (2 + Sk[t][Mc*j+2]);
						Est_sym_Q[0][j][t] = (-1.0/sqrt(10.0)) * Sk[t][Mc*j+1] * (2 + Sk[t][Mc*j+3]);
					}
				}
			}

			for(s2=1; s2<Iterate2; s2++)
			{
				// PPIC ICI Cancellation
				for(t=0; t<N_fft; t++)
					for(i=0; i<R; i++)
					{
						// ICI Interference Reconstruction
						Re_ICI_I[i][t] = Re_ICI_Q[i][t] = 0.0;
						s3 = s2-1;
						for(d=1; d<=s2; d++)
						{
							for(j=0; j<T; j++)
							{
								Re_ICI_I[i][t] += (Est_sym_I[s3][j][(N_fft+t-d)%N_fft] * H_I[i][j][t][d] - Est_sym_Q[s3][j][(N_fft+t-d)%N_fft] * H_Q[i][j][t][d]);
								Re_ICI_Q[i][t] += (Est_sym_I[s3][j][(N_fft+t-d)%N_fft] * H_Q[i][j][t][d] + Est_sym_Q[s3][j][(N_fft+t-d)%N_fft] * H_I[i][j][t][d]);
							}
							s3--;
						}

						s3 = s2-1;
						for(d=N_fft-1; d>=N_fft-s2; d--)
						{
							for(j=0; j<T; j++)
							{
								Re_ICI_I[i][t] += (Est_sym_I[s3][j][(N_fft+t-d)%N_fft] * H_I[i][j][t][d] - Est_sym_Q[s3][j][(N_fft+t-d)%N_fft] * H_Q[i][j][t][d]);
								Re_ICI_Q[i][t] += (Est_sym_I[s3][j][(N_fft+t-d)%N_fft] * H_Q[i][j][t][d] + Est_sym_Q[s3][j][(N_fft+t-d)%N_fft] * H_I[i][j][t][d]);
							}
							s3--;
						}

						// ICI Cancellation
						Re_signal[i][t][0] = Yk[i][t][0] - Re_ICI_I[i][t];
						Re_signal[i][t][1] = Yk[i][t][1] - Re_ICI_Q[i][t];
					}

				// MRC
				for(t=0; t<N_fft; t++)
					for(i=0; i<R; i++)
					{
						MRC_out[i][t][0] = (Re_signal[i][t][0] * H_I[i][j][t][0] + Re_signal[i][t][1] * H_Q[i][j][t][0]);
						MRC_out[i][t][1] = (Re_signal[i][t][1] * H_I[i][j][t][0] - Re_signal[i][t][0] * H_Q[i][j][t][0]);
						for(d=1; d<=s2; d++)
							for(j=0; j<T; j++)
							{
								MRC_out[i][t][0] += (Est_sym_I[s2-1][j][t] * (pow(H_I[i][j][(N_fft+t-d)%N_fft][d],2) +  pow(H_Q[i][j][(N_fft+t-d)%N_fft][d], 2)));
								MRC_out[i][t][1] += (Est_sym_Q[s2-1][j][t] * (pow(H_I[i][j][(N_fft+t-d)%N_fft][d],2) +  pow(H_Q[i][j][(N_fft+t-d)%N_fft][d], 2)));
								
								MRC_out[i][t][0] += (Est_sym_I[s2-1][j][t] * (pow(H_I[i][j][(t+d)%N_fft][N_fft-d],2) +  pow(H_Q[i][j][(t+d)%N_fft][N_fft-d], 2)));
								MRC_out[i][t][1] += (Est_sym_Q[s2-1][j][t] * (pow(H_I[i][j][(t+d)%N_fft][N_fft-d],2) +  pow(H_Q[i][j][(t+d)%N_fft][N_fft-d], 2)));
							}
					}

				// Soft
			}

			p++;
		} while(err_count[Iterate2-1][0] <= 1000 && p < num_packet);

		// Statistics and records
		cout << "Error Rate = " << endl;
		for(s2=0; s2<Iterate2; s2++)
		{
			for(s1=0; s1<Iterate1; s1++)
			{
      			err_rate[s2][s1] = err_count[s2][s1] / (double)(Mc*T*p*N_fft);
      			printf("%e, ", err_rate[s2][s1]);
			}
			cout << endl;
		}
		cout << endl;

		fprintf(ber, "%f ", Eb_No);
		fprintf(records, "%f:\n", Eb_No);
		for(s2=0; s2<Iterate2; s2++)
		{
			for(s1=0; s1<Iterate1; s1++)
			{
				fprintf(ber, "%e ", err_rate[s2][s1]);
				fprintf(records, "%e ", err_rate[s2][s1]);
			}
			fprintf(records, "\n");
		}
		fprintf(ber, "\n");
		fprintf(records, "\n");
		fflush(records);
		fflush(ber);
	}

	for(i=0; i<N_fft; i++)
		delete data_bit[i];
	delete data_bit;
	delete data;
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
	delete Ak;
	for(t=0; t<N_fft; t++)
		delete Sk[t];
	delete Sk;
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			delete Yk[i][l];
	for(i=0; i<R; i++)
   		delete Yk[i];
	delete Yk;
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			delete Re_signal[i][l];
	for(i=0; i<R; i++)
   		delete Re_signal[i];
	delete Re_signal;
	for(k=0; k<N_fft; k++)
		for(i=0; i<R; i++)
   			for(j=0; j<Mc*T; j++)
   				delete R_pn[k][i][j];
	for(k=0; k<N_fft; k++)
		for(i=0; i<R; i++)
   			delete R_pn[k][i];
	for(k=0; k<N_fft; k++)
   		delete R_pn[k];
	delete R_pn;
	for(k=0; k<N_fft; k++)
		for(i=0; i<Mc*T; i++)
   			for(j=0; j<R; j++)
   				delete Q_np[k][i][j];
	for(k=0; k<N_fft; k++)
		for(i=0; i<Mc*T; i++)
   			delete Q_np[k][i];
	for(k=0; k<N_fft; k++)
		delete Q_np[k];
	delete Q_np;
	for(k=0; k<N_fft; k++)
		for(i=0; i<Mc*T; i++)
   			delete Out[k][i];
	for(k=0; k<N_fft; k++)
   		delete Out[k];
	delete Out;
	for(i=0; i<Iterate2; i++)
		delete err_count[i];
	delete err_count;
	for(i=0; i<Iterate2; i++)
		delete err_rate[i];
	delete err_rate;
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
	for(i=0; i<Iterate2; i++)
		for(j=0; j<T; j++)
			delete Est_sym_I[i][j];
	for(i=0; i<Iterate2; i++)
		delete Est_sym_I[i];
	delete Est_sym_I;
	for(i=0; i<Iterate2; i++)
		for(j=0; j<T; j++)
			delete Est_sym_Q[i][j];
	for(i=0; i<Iterate2; i++)
		delete Est_sym_Q[i];
	delete Est_sym_Q;
	for(i=0; i<R; i++)
		delete Re_ICI_I[i];
	delete Re_ICI_I;
	for(i=0; i<R; i++)
		delete Re_ICI_Q[i];
	delete Re_ICI_Q;
	delete I;
	delete Q;

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
