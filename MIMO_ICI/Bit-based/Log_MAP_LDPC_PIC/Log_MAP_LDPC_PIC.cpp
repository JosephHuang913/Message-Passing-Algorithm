/*========================================================*/
/* Author: Chao-wang Huang                                                                       */
/* Date: Wednesday, April 7, 2009                                                              */
/* A Bit-based Message Passinlg Algorithm for MIMO Channel is simulated */
/* 2 by 2 QPSK                                                                                           */
/* MIMO Detector: Log-domain MPA with Sum-Product Rule                       */
/* LDPC code: Gallager (20, 4, 3)                                                               */
/* Parallel Interference Cancellation (PIC)                                                  */
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
void Multipath_fade_pattern(long double *****, int);
//void Multipath_fade_pattern(long double ****, int);
void LDPC_Encode(int *, int *);
int Belief_Propagation(long double *, long double *, long double, int, int);
long double F(long double);
void dfour1(long double *, unsigned long, int);

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 	3.14159265358979
#define Num_packet 2500					// number of packets simulated
//#define Num_codeword 409			// number of codewords per packet (16-QAM)
#define Num_codeword 204				// number of codewords per packet (QPSK)
//const int Mc = 4;								// modulation order (4 for 16QAM)
const int Mc = 2;								// modulation order (2 for QPSK)
const int T = 2;									// Number of transmit antenna
const int R = 2;									// Number of receive antenna
const int Iterate2 = 10;						// For PPIC
const int Iterate1 = 1;						// For Message Passing MIMO Detection
const int Iterate3 = 1;						// For LDPC Decoding
const int N_fft =  1024;
const int CP_length = N_fft/8;				// Length of Cyclic Prefix
const long double sample_rate = 11.2e6;		/* Sampling Frequency */
/* Power Delay Profile of ITU Vehicular A Channel */
const int path_num = 6;
const int path_delay[6] = {0, 3, 8, 12, 19, 28};	// Delay (sample) {0,310,710,1090,1730,2510} (ns)
const long double path_gain[6] = {0.485, 0.3853, 0.0611, 0.0485, 0.0153, 0.0049};	// Power Gain (linear scale)
const int delay_spread = 28;				// Maximum path delay (samples)
const long double vc = 350.0;						/* speed of vehicle in km/hr */
const long double C = 3.0e8;						/* light speed */
const long double fc = 2.5e9;						/* carrier frequency */
//const long double OFDM_sym_duration = 91.43e-6;		/* OFDMA Symbol Duration without Cyclic Prefix */
const long double OFDM_sym_duration = 102.86e-6;		/* OFDMA Symbol Duration with Cyclic Prefix */
const long double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
const long double fdxT = Doppler * OFDM_sym_duration;
#define NODES          20							// Maximum of number of code/check nodes
#define MAX_CHK      3   						// Maximum number of checks per code bit
#define MAX_BIT        4   						// Maximum number of code bits per check
int Max_col_weight;                 					// Maximum column weight
int Max_row_weight;                 				// Maximum row weight
const int Cn = 20;               						// Codeword length of a linear block code
const int Ck = 7;                						// Information bits of a linear block code
int M, N;                          							// Size of parity-check matrix (Row, Column) = (M, N)
int **G, **H;												// Generator matrix & parity check matrix
int s1, s2;
long double ***ts;
#define N_o		8
double theta_n[N_o];

// ---------------
// NODE STRUCTURES
// ---------------
struct bit_node {
	int Num_chks_bit;									// Number of check nodes of a bit node
    int index[MAX_CHK];							// Check nodes set of a bit node
    long double pi1[MAX_CHK];							// messages "pi" to check node
    };

struct check_node {
    int Num_bits_chk;									// Number of bit nodes of a check node		
    int index[MAX_BIT];								// Bit nodes set of a check node
    long double lambda1[MAX_BIT];					// messages "lambda" to bit node
    };

struct bit_node			Bit_node[Num_codeword][NODES];
struct check_node	CHK_node[Num_codeword][NODES];

int main(void)
{
	time_t  ti, start, end;
	int i, j, k, q, l, m, p, r, x, z, t, d, *data, **data_bit, **Ak, *S_random;
	int **err_count, *delay, **coded_bit, flag, *interleaved, *inter, *coded, *bit, *de_inter;
	long double Eb_No, snr, noise_pwr, noise[2], Y[2], ***Yk, **err_rate, Ps;
	long double sum, pdf, ***L_r_pn, ***L_q_np, **sym_I, **sym_Q, *I, *Q;
	long double *gain, *****fade, ***F_l, ****H_I, ****H_Q, ***Est_sym_I, ***Est_sym_Q;
	long double **Re_ICI_I, **Re_ICI_Q, ***Re_signal, **Sk, *fft, *Out;
	long double **L_q, L_sum, L_tmp[2], *L_a_u, *L_a_d, *L_in, *L_out, *L_o;
	long double ****ICI_gain, sigma_e, ****H_est_I, ****H_est_Q;
	FILE *ber, *records, *s_inter, *parity_chk;

	start = time(NULL);
	printf("BER Performance of Bit-based Message Passinlg in MIMO Channel\n");
	printf("Gallager (20, 7) LDPC code\n");
	printf("Dimension of Parity Check matrix: (15,20)\n");
	printf("Column Weight = 3, Row Weight = 4\n");
	printf("Number of transmit antenna is %d\n", T);
	printf("Number of receive antenna is %d\n", R);
	printf("Vehicular speed: %f km/h\n", vc);
	printf("Maximum Doppler frequency: %f\n", Doppler);
	printf("Channel fading rate: %f\n", fdxT);
	printf("Maximum number of bits of simulation = %d\n\n", Cn*Num_codeword*Num_packet);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Record_Bit-based_Message_Passinlg.log", "a");
	fprintf(records, "BER Performance of Bit-based Message Passinlg in MIMO Channel\n");
	fprintf(records, "Gallager (20, 7) LDPC code\n");
	fprintf(records, "Dimension of Parity Check matrix: (15,20)\n");
	fprintf(records, "Column Weight = 3, Row Weight = 4\n");
	fprintf(records, "Number of transmit antenna is %d\n", T);
	fprintf(records, "Number of receive antenna is %d\n", R);
	fprintf(records, "Vehicular speed: %f km/h\n", vc);
	fprintf(records, "Maximum Doppler frequency: %f\n", Doppler);
	fprintf(records, "Channel fading rate: %f\n", fdxT);
	fprintf(records, "Maximum number of bits of simulation = %d\n\n", Cn*Num_codeword*Num_packet);
	fprintf(records, "Eb/No     BER\n");
	fflush(records);

	data_bit = new int*[Num_codeword];
	for(i=0; i<Num_codeword; i++)
		data_bit[i] = new int[Ck];
	data = new int[Mc*T];
	coded_bit = new int*[Num_codeword];
	for(i=0; i<Num_codeword; i++)
		coded_bit[i] = new int[Cn];
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
	Ak = new int*[Num_codeword];
	for(i=0; i<Num_codeword; i++)
		Ak[i] = new int[Cn];
	Sk = new long double*[N_fft];
	for(i=0; i<N_fft; i++)
		Sk[i] = new long double[Mc*T];
	fade = new long double****[R];
	for(i=0; i<R; i++)
		fade[i] = new long double***[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			fade[i][j] = new long double**[path_num];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				fade[i][j][l] = new long double*[N_fft+CP_length];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				for(t=0; t<N_fft+CP_length; t++)
					fade[i][j][l][t] = new long double[2];
	Yk = new long double**[R];
	for(i=0; i<R; i++)
   		Yk[i] = new long double*[N_fft];
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			Yk[i][l] = new long double[2];
	Re_signal = new long double**[R];
	for(i=0; i<R; i++)
   		Re_signal[i] = new long double*[N_fft];
	for(i=0; i<R; i++)
		for(l=0; l<N_fft; l++)
			Re_signal[i][l] = new long double[2];
	F_l = new long double**[path_num];
	for(i=0; i<path_num; i++)
		F_l[i] = new long double*[N_fft];
	for(i=0; i<path_num; i++)
		for(l=0; l<N_fft; l++)
			F_l[i][l] = new long double[2];
	H_I = new long double***[R];
	for(i=0; i<R; i++)
  		H_I[i] = new long double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_I[i][j] = new long double*[N_fft];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				H_I[i][j][l] = new long double[N_fft];
	H_Q = new long double***[R];
	for(i=0; i<R; i++)
   		H_Q[i] = new long double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_Q[i][j] = new long double*[N_fft];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				H_Q[i][j][l] = new long double[N_fft];
	H_est_I = new long double***[R];
	for(i=0; i<R; i++)
  		H_est_I[i] = new long double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_est_I[i][j] = new long double*[N_fft];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				H_est_I[i][j][l] = new long double[N_fft];
	H_est_Q = new long double***[R];
	for(i=0; i<R; i++)
   		H_est_Q[i] = new long double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_est_Q[i][j] = new long double*[N_fft];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				H_est_Q[i][j][l] = new long double[N_fft];
	ICI_gain = new long double***[R];
	for(i=0; i<R; i++)
   		ICI_gain[i] = new long double**[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			ICI_gain[i][j] = new long double*[N_fft];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				ICI_gain[i][j][l] = new long double[N_fft];
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
	err_count = new int*[Iterate2];
	for(i=0; i<Iterate2; i++)
		err_count[i] = new int[Iterate1];	
	err_rate = new long double*[Iterate2];
	for(i=0; i<Iterate2; i++)
		err_rate[i] = new long double[Iterate1];
	Est_sym_I = new long double**[Iterate2];
	for(i=0; i<Iterate2; i++)
		Est_sym_I[i] = new long double*[T];
	for(i=0; i<Iterate2; i++)
		for(j=0; j<T; j++)
			Est_sym_I[i][j] = new long double[N_fft];
	Est_sym_Q = new long double**[Iterate2];
	for(i=0; i<Iterate2; i++)
		Est_sym_Q[i] = new long double*[T];
	for(i=0; i<Iterate2; i++)
		for(j=0; j<T; j++)
			Est_sym_Q[i][j] = new long double[N_fft];
	Re_ICI_I = new long double*[R];
	for(i=0; i<R; i++)
		Re_ICI_I[i] = new long double[N_fft];
	Re_ICI_Q = new long double*[R];
	for(i=0; i<R; i++)
		Re_ICI_Q[i] = new long double[N_fft];
	I = new long double[T];
	Q = new long double[T];
	H = new int*[15];
	for(i=0; i<15; i++)
		H[i] = new int[20];
	G = new int*[Ck];
	for(i=0; i<Ck; i++)
		G[i] = new int[Cn];
	L_a_u = new long double[N_fft*T*Mc];
	L_a_d = new long double[N_fft*T*Mc];
	L_o = new long double[N_fft*T*Mc];
	L_in = new long double[Cn];
	L_out = new long double[Cn];
	fft = new long double[2*N_fft+1];
	interleaved = new int[N_fft*T*Mc];
	inter = new int[N_fft*T*Mc];
	coded = new int[Cn];
	bit = new int[Ck];
	S_random = new int[N_fft*T*Mc];
	de_inter = new int[N_fft*T*Mc];
	Out = new long double[N_fft*T*Mc];

	// Channel Weighting Gain and Channel path delay
	for(i=0; i<path_num; i++)
	{
		delay[i] = path_delay[i];
		gain[i] = sqrtl(path_gain[i]);
	}

	srand((unsigned) time(&ti));

	// Define Parity Check Matrix
	fopen_s(&parity_chk, "gallager_20_4_3.log", "r");
	fscanf_s(parity_chk, "%d %d", &N, &M);
	fscanf_s(parity_chk, "%d %d", &Max_row_weight, &Max_col_weight);

	for (i=0; i<M; i++)
	    fscanf_s(parity_chk, "%d", &CHK_node[0][i].Num_bits_chk);
	for (i=0; i<N; i++)
	    fscanf_s(parity_chk, "%d", &Bit_node[0][i].Num_chks_bit);
	
		// Read index sets for check nodes
	for (i=0; i<M; i++)
		for (j=0; j<CHK_node[0][i].Num_bits_chk; j++)
			fscanf_s(parity_chk, "%d", &CHK_node[0][i].index[j]);
    
	   // Read index sets for bit nodes
	for (i=0; i<N; i++)
		for (j=0; j<Bit_node[0][i].Num_chks_bit; j++)
			fscanf_s(parity_chk, "%d", &Bit_node[0][i].index[j]);

	fclose(parity_chk);

	for(i=0; i<M; i++)
		for(j=0; j<N; j++)
			H[i][j] = 0;

   for(i=0; i<M; i++)
	   for(j=0; j<CHK_node[0][i].Num_bits_chk; j++)
		   H[i][CHK_node[0][i].index[j]-1] = 1;

   fopen_s(&parity_chk, "G.log", "r");
	for(i=0; i<Ck; i++)
		for(j=0; j<Cn; j++)
			fscanf_s(parity_chk, "%d", &G[i][j]);
	fclose(parity_chk);

/*===========================================*/
/* S-random interleaver and de-interleaver (S = 64 or 45) */
/*===========================================*/
	fopen_s(&s_inter, "s_random.log", "r");
	for(i=0; i<N_fft*T*Mc; i++)		// Interleaver
   		fscanf_s(s_inter, "%d %d", &de_inter[i], &S_random[i]);
	for(i=0; i<N_fft*T*Mc; i++)		// De-interleaver
   		de_inter[S_random[i]] = i;
	fclose(s_inter);

/*======================*/
/* main simulation loop            */
/*======================*/
	fopen_s(&ber, "BER_MPA_LDPC.log", "w");
	for(snr=0; snr<=16; snr+=2)
	{
   		for(s2=0; s2<Iterate2; s2++)
			for(s1=0; s1<Iterate1; s1++)
   				err_count[s2][s1] = 0;

   		// noise powler calculation
		Eb_No = (double)snr;
		Ps = 1 * T;
		//noise_pwr = 0.5*Ps*R/((float)N_fft*T*((double)Ck/(double)Cn)*powl(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Time Domain)
		noise_pwr = 0.5*Ps*R/((float)T*((double)Ck/(double)Cn)*powl(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Freq. Domain)
		//noise_pwr = 0.25*Ps*R/((float)T*((double)Ck/(double)Cn)*powl(10.0, Eb_No/10.0));	// 16QAM, Nyquist filter assumption (Freq. Domain)
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
			// Generate random information bit 
			for(i=0; i<Num_codeword; i++)
			{
				for(j=0; j<Ck; j++)
					if(rand()/(long double)RAND_MAX>=0.5)
					{
						bit[j] = 1;
						data_bit[i][j] = bit[j];
					}
					else
					{
						bit[j] = 0;
						data_bit[i][j] = bit[j];
					}

				// LDPC encoding: Gallager (20, 4, 3)
				LDPC_Encode(bit, coded);
				for(l=0; l<Cn; l++)
				{
					coded_bit[i][l] = coded[l];
					inter[i*Cn+l] = coded[l];
				}
			}

			// Zero padding
			for(l=Num_codeword*Cn; l<N_fft*T*Mc; l++)
				inter[l] = 0;

			// Interleaving: (8192, 64) S-random interleaver
			// Interleaving: (4096, 45) S-random interleaver
			for(i=0; i<N_fft*T*Mc; i++)
      			interleaved[S_random[i]] = inter[i];
				//interleaved[i] = inter[i]; //		No interleaving

			for(l=0; l<N_fft; l++)
   				for(i=0; i<T; i++)
	   			{
					// QPSK Mapping
					sym_I[i][l] = (2*interleaved[l*Mc*T+i*Mc]-1)/sqrtl(2.0);
         			sym_Q[i][l] = (2*interleaved[l*Mc*T+i*Mc+1]-1)/sqrtl(2.0);
/*
					// 16-QAM Mapping	
					if(interleaved[l*Mc*T+i*Mc] == 0 && interleaved[l*Mc*T+i*Mc+2] == 0)
   		  				sym_I[i][l] = 1.0/sqrtl(10.0);
					else if(interleaved[l*Mc*T+i*Mc] == 0 && interleaved[l*Mc*T+i*Mc+2] == 1)
   			  			sym_I[i][l] = 3.0/sqrtl(10.0);
					else if(interleaved[l*Mc*T+i*Mc] == 1 && interleaved[l*Mc*T+i*Mc+2] == 1)
		  		 		sym_I[i][l] = -3.0/sqrtl(10.0);
				  	else
  		   				sym_I[i][l] = -1.0/sqrtl(10.0);

		   			if(interleaved[l*Mc*T+i*Mc+1] == 0 && interleaved[l*Mc*T+i*Mc+3] == 0)
  	   					sym_Q[i][l] = 1.0/sqrtl(10.0);
					else if(interleaved[l*Mc*T+i*Mc+1] == 0 && interleaved[l*Mc*T+i*Mc+3] == 1)
   				  		sym_Q[i][l] = 3.0/sqrtl(10.0);
				     else if(interleaved[l*Mc*T+i*Mc+1] == 1 && interleaved[l*Mc*T+i*Mc+3] == 1)
	   			 		sym_Q[i][l] = -3.0/sqrtl(10.0);
					 else
      		 			sym_Q[i][l] = -1.0/sqrtl(10.0);*/
				}
				
			// Random Phase of Rayleigh fading pattern
			if(p%10 == 0)
				for(i=0; i<N_o; i++)
					theta_n[i] = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase

			// Multi-path Rayleigh fading pattern with ICI
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
								F_l[l][k][0] += gain[l]*fade[i][j][l][t+CP_length-delay[l]][0] * cosl(2*Pi*t*k/(float)N_fft) + gain[l]*fade[i][j][l][t+CP_length-delay[l]][1] * sinl(2*Pi*t*k/(float)N_fft);
								F_l[l][k][1] += gain[l]*fade[i][j][l][t+CP_length-delay[l]][1] * cosl(2*Pi*t*k/(float)N_fft) - gain[l]*fade[i][j][l][t+CP_length-delay[l]][0] * sinl(2*Pi*t*k/(float)N_fft);
							}
						}

					// Channel Frequency Response & ICI Channel Coefficients
					for(t=0; t<N_fft; t++)
						for(d=0; d<N_fft; d++)
						{
							H_I[i][j][t][d] = H_Q[i][j][t][d] = 0.0;
							for(l=0; l<path_num; l++)
							{
								H_I[i][j][t][d] += F_l[l][d][0] * cosl(2*Pi*delay[l]*(t-d)/(float)N_fft) + F_l[l][d][1] * sinl(2*Pi*delay[l]*(t-d)/(float)N_fft);
								H_Q[i][j][t][d] += F_l[l][d][1] * cosl(2*Pi*delay[l]*(t-d)/(float)N_fft) - F_l[l][d][0] * sinl(2*Pi*delay[l]*(t-d)/(float)N_fft);
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

			// Channel Estimation Error
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
					for(t=0; t<N_fft; t++)
					{
						sum = 0.0;
						for(d=0; d<N_fft; d++)
							sum += powl(H_I[i][j][t][d], 2.0) + powl(H_Q[i][j][t][d], 2.0);

						for(d=0; d<N_fft; d++)
							ICI_gain[i][j][t][d] = sqrtl(powl(H_I[i][j][t][d], 2.0)+powl(H_Q[i][j][t][d], 2.0) / sum);
					}

/*========================================================*/
/*  Log-Domain Bit-based Message Passinlg Algorithm & PIC                     */
/* Channel Decoder: Gallager (20, 4, 3) Regular LDPC code                      */
/*========================================================*/
			
			// Variance of Channel Estimation Error
			sigma_e = 0.03 + 0.8 / powl(10.0, Eb_No/10.0);
			//sigma_e = 0.1;
			//sigma_e = 0.01;

			/*
			// Channel Estimation Error
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
					for(t=0; t<N_fft; t++)
					{
						AWGN_noise(0.0, sigma_e, &noise[0]);

						for(d=0; d<N_fft; d++)
						{
							H_est_I[i][j][t][d] = H_I[i][j][t][d] + noise[0] * ICI_gain[i][j][t][d];
							H_est_Q[i][j][t][d] = H_Q[i][j][t][d] + noise[1] * ICI_gain[i][j][t][d];
						}
					}
*/
			// MIMO Detector: LLR Initialization
			for(t=0; t<N_fft; t++)
				for(j=0; j<T; j++)
					for(m=0; m<Mc; m++)
						for(i=0; i<R; i++)
							L_q_np[t][j*Mc+m][i] = 0.0;

			for(s2=0; s2<Iterate2; s2++)
			{
				// Channel Estimation Error
				for(i=0; i<R; i++)
					for(j=0; j<T; j++)
						for(t=0; t<N_fft; t++)
						{
							AWGN_noise(0.0, sigma_e, &noise[0]);

							for(d=0; d<N_fft; d++)
							{
								H_est_I[i][j][t][d] = H_I[i][j][t][d] + noise[0] * ICI_gain[i][j][t][d];
								H_est_Q[i][j][t][d] = H_Q[i][j][t][d] + noise[1] * ICI_gain[i][j][t][d];
							}
						}

				sigma_e /= 2.0;

				// PIC ICI Cancellation
				if(s2 != 0)
				{
					for(t=0; t<N_fft; t++)
						for(i=0; i<R; i++)
						{
							// ICI Interference Reconstruction
							Re_ICI_I[i][t] = Re_ICI_Q[i][t] = 0.0;
							for(d=1; d<=6; d++)
								for(j=0; j<T; j++)
								{
									Re_ICI_I[i][t] += (Est_sym_I[s2-1][j][(N_fft+t-d)%N_fft] * H_est_I[i][j][t][d] - Est_sym_Q[s2-1][j][(N_fft+t-d)%N_fft] * H_est_Q[i][j][t][d]);
									Re_ICI_Q[i][t] += (Est_sym_I[s2-1][j][(N_fft+t-d)%N_fft] * H_est_Q[i][j][t][d] + Est_sym_Q[s2-1][j][(N_fft+t-d)%N_fft] * H_est_I[i][j][t][d]);
								}

							for(d=N_fft-1; d>=N_fft-6; d--)
								for(j=0; j<T; j++)
								{
									Re_ICI_I[i][t] += (Est_sym_I[s2-1][j][(N_fft+t-d)%N_fft] * H_est_I[i][j][t][d] - Est_sym_Q[s2-1][j][(N_fft+t-d)%N_fft] * H_est_Q[i][j][t][d]);
									Re_ICI_Q[i][t] += (Est_sym_I[s2-1][j][(N_fft+t-d)%N_fft] * H_est_Q[i][j][t][d] + Est_sym_Q[s2-1][j][(N_fft+t-d)%N_fft] * H_est_I[i][j][t][d]);
								}

							// ICI Cancellation
							Re_signal[i][t][0] = Yk[i][t][0] - Re_ICI_I[i][t];
							Re_signal[i][t][1] = Yk[i][t][1] - Re_ICI_Q[i][t];

							// Obtain frequency diversity
							for(d=1; d<=6; d++)
								for(j=0; j<T; j++)
								{
									Re_signal[i][t][0] += (Est_sym_I[s2-1][j][t] * H_est_I[i][j][(N_fft+t-d)%N_fft][d] - Est_sym_Q[s2-1][j][t] * H_est_Q[i][j][(N_fft+t-d)%N_fft][d]);
									Re_signal[i][t][1] += (Est_sym_I[s2-1][j][t] * H_est_Q[i][j][(N_fft+t-d)%N_fft][d] + Est_sym_Q[s2-1][j][t] * H_est_I[i][j][(N_fft+t-d)%N_fft][d]);
									
									Re_signal[i][t][0] += (Est_sym_I[s2-1][j][t] * H_est_I[i][j][(t+d)%N_fft][N_fft-d] - Est_sym_Q[s2-1][j][t] * H_est_Q[i][j][(t+d)%N_fft][N_fft-d]);
									Re_signal[i][t][1] += (Est_sym_I[s2-1][j][t] * H_est_Q[i][j][(t+d)%N_fft][N_fft-d] + Est_sym_Q[s2-1][j][t] * H_est_I[i][j][(t+d)%N_fft][N_fft-d]);
								}
						}
				}

				for(s1=0; s1<Iterate1; s1++)
				{
					for(t=0; t<N_fft; t++)
					{
						// Extrinsic information of function nodes
			    		for(i=0; i<R; i++)
							for(j=0; j<T; j++)
								for(m=0; m<Mc; m++)
           						{
									flag = 0;
									for(r=Num_codeword*Cn; r<N_fft*Mc*T; r++)
										if(t*Mc*T+j*Mc+m == S_random[r])			
										{
											L_r_pn[t][i][j*Mc+m] = -50.0;
											flag = 1;
										}
									
									if(flag == 0)
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
														for(d=Num_codeword*Cn; d<N_fft*Mc*T; d++)
															if(t*Mc*T+r == S_random[d])
																data[r] = 0;
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
							       					Y[0] += (I[r] * H_est_I[i][r][t][0] - Q[r] * H_est_Q[i][r][t][0]);
													Y[1] += (I[r] * H_est_Q[i][r][t][0] + Q[r] * H_est_I[i][r][t][0]);
												}
										
								 				L_sum = 0.0;
												if(s2 != 0 || s1 != 0)
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
													if(s2 == 0)
							        					sum += powl(Yk[i][t][z]-Y[z],2.0);
													else
														sum += powl(Re_signal[i][t][z]-Y[z],2.0);

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
					}

					// MIMO Detector to LDPC Decoder
					for(t=0; t<N_fft; t++)
						for(j=0; j<T; j++)
							for(m=0; m<Mc; m++)
							{	
								L_o[t*Mc*T+j*Mc+m] = 0.0;
								for(i=0; i<R; i++)
									L_o[t*Mc*T+j*Mc+m] += L_r_pn[t][i][j*Mc+m];

								if (L_o[t*Mc*T+j*Mc+m] < -50.0)
									L_o[t*Mc*T+j*Mc+m] = -50.0;
								if (L_o[t*Mc*T+j*Mc+m] > 50.0)
									L_o[t*Mc*T+j*Mc+m] = 50.0;
							}

					// De-interleaving
       				for(i=0; i<N_fft*T*Mc; i++)
         				L_a_u[de_inter[i]] = L_o[i];
						//L_a_u[i] = L_o[i];		// No de-interleaving
				
					// LDPC Decoding
					for(i=0; i<Num_codeword; i++)
					{
						for(l=0; l<Cn; l++)
							L_in[l] = L_a_u[i*Cn+l];

						Belief_Propagation(L_in,  L_out, snr, Iterate3, i);
						for(l=0; l<Cn; l++)
							L_o[i*Cn+l] = L_out[l];

						// Hard Decision
						for(l=0; l<Cn; l++)
						{
					 		if( (L_o[i*Cn+l]+L_a_u[i*Cn+l]) >= 0.0)
								Ak[i][l] = 1;
       						else
       							Ak[i][l] = 0;

	            			err_count[s2][s1] += Error_count(coded_bit[i][l], Ak[i][l]);
						}
					}

					for(l=Num_codeword*Cn; l<N_fft*T*Mc; l++)
						L_o[l] = -50.0;

					// Interleaving
					for(i=0; i<N_fft*T*Mc; i++)
					{
						Out[S_random[i]] = L_o[i] + L_a_u[i];
						//Out[S_random[i]] = L_o[i];
      					L_a_d[S_random[i]] = L_o[i];
						//L_a_d[i] = L_o[i];		// No interleaving
					}

					// LDPC Decoder to MIMO Detector
					// Extrinsic information of bit nodes
					for(t=0; t<N_fft; t++)
						for(j=0; j<T; j++)
							for(m=0; m<Mc; m++)
							{
								flag = 0;
								for(i=Num_codeword*Cn; i<N_fft*Mc*T; i++)
									if(t*Mc*T+j*Mc+m == S_random[i])
										for(l=0; l<R; l++)
										{
											L_q_np[t][j*Mc+m][l] = -50.0;
											flag = 1;
										}
								
								if(flag == 0)
									for(i=0; i<R; i++)
									{
										L_q_np[t][j*Mc+m][i] = 0.0;
				   						for(r=0; r<R; r++)
		       								if(r != i)
												L_q_np[t][j*Mc+m][i] += L_r_pn[t][r][j*Mc+m];
	
										L_q_np[t][j*Mc+m][i] += L_a_d[t*Mc*T+j*Mc+m];
										
										if (L_q_np[t][j*Mc+m][i] < -50.0)
											L_q_np[t][j*Mc+m][i] = -50.0;
										if (L_q_np[t][j*Mc+m][i] > 50.0)
											L_q_np[t][j*Mc+m][i] = 50.0;
				 					}
							}
				}

				// For ICI Cancellation
				for(t=0; t<N_fft; t++)
				{
					// Bit Soft Decision
					for(j=0; j<Mc*T; j++)
            			Sk[t][j] = tanhl(0.5*Out[t*Mc*T+j]);

					for(j=0; j<T; j++)
					{
						// QPSK Soft Mapping
						Est_sym_I[s2][j][t] = Sk[t][Mc*j] / sqrtl(2.0);
						Est_sym_Q[s2][j][t] = Sk[t][Mc*j+1] / sqrtl(2.0);

						// 16-QAM Soft Mapping
						//Est_sym_I[s2][j][t] = (-1.0/sqrtl(10.0)) * Sk[t][Mc*j] * (2 + Sk[t][Mc*j+2]);
						//Est_sym_Q[s2][j][t] = (-1.0/sqrtl(10.0)) * Sk[t][Mc*j+1] * (2 + Sk[t][Mc*j+3]);
					}
				}
			}

			p++;
		} while(err_count[Iterate2-1][Iterate1-1] <= 10000 && p < Num_packet);

		// Statistics and records
		cout << "Error Rate = " << endl;
		for(s2=0; s2<Iterate2; s2++)
		{
			for(s1=0; s1<Iterate1; s1++)
			{
      			err_rate[s2][s1] = err_count[s2][s1] / (long double)(Cn*Num_codeword*p);
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

	for(i=0; i<Num_codeword; i++)
		delete data_bit[i];
	delete data_bit;
	delete data;
	for(i=0; i<Num_codeword; i++)
		delete coded_bit[i];
	delete coded_bit;
	for(i=0; i<Num_codeword; i++)
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
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				delete H_est_I[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete H_est_I[i][j];
	for(i=0; i<R; i++)
   		delete H_est_I[i];
	delete H_est_I;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				delete H_est_Q[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete H_est_Q[i][j];
	for(i=0; i<R; i++)
		   delete H_est_Q[i];
	delete H_est_Q;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<N_fft; l++)
				delete ICI_gain[i][j][l];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete ICI_gain[i][j];
	for(i=0; i<R; i++)
		   delete ICI_gain[i];
	delete ICI_gain;
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
	for(i=0; i<M; i++)
		delete H[i];
	delete H;
	for(i=0; i<Ck; i++)
		delete G[i];
	delete G;
	delete L_a_u;
	delete L_a_d;
	delete L_o;
	delete L_in;
	delete L_out;
	delete fft;
	delete interleaved;
	delete inter;
	delete coded;
	delete bit;
	delete S_random;
	delete de_inter;
	delete Out;

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
   int n, Np/*, N_o = 8*/;
   long double lamda, w_m, beta_n, w_n, alpha, T_c2, T_s2/*, theta_n*/;

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
            //theta_n = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
			//theta_n = 0.0;
			*(fade+0) += sqrt(2.0/T_c2)*cos(beta_n)*cos(w_n*t+theta_n[n-1]);
			*(fade+1) += sqrt(2.0/T_s2)*sin(beta_n)*cos(w_n*t+theta_n[n-1]);
		}
	}
}

void Multipath_fade_pattern(long double *****fade, int path_num)
{
	int i, j, l, t;
	long double gain[2];

	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			for(l=0; l<path_num; l++)
				for(t=0; t<N_fft+CP_length; t++)
				{
					JakesFading(fc, vc*1000/3600.0, ts[i][j][l], 2, &gain[0]);
					ts[i][j][l] += OFDM_sym_duration / float(N_fft+CP_length);
					fade[i][j][l][t][0] = gain[0];
					fade[i][j][l][t][1] = gain[1];
				}
}
/*
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
*/
int Belief_Propagation(long double *received, long double *decoded, long double snr, int max_iter, int ind)
{
// Iterative decoding by belief propagation in code's Bayesian network
// Based on Pearl's book and MacKay's paper

	int i, j, l, iter, m, aux, sign;
	long double delt;
	long double llrp1[NODES];									// Prior probabilities (channel)
	long double q1[NODES];										// Pseudo-posterior probabilities
	//long double snr_rms = 2.0 * sqrtl(2.0*(long double)Ck/(long double)Cn*(powl(10.0,(snr/10.0))));

	// -------------------
	// ***** STEP 0 *****
	// INITIALIZATION STEP
	// -------------------

	// Prior logl-likelihood ratios (channel metrics)
	for (i=0; i<N; i++)
    {
		// LOOK-UP TABLE (LUT)
		//llrp1[i] = received[i]*snr_rms;	// ln{P(1)/P(-1)}
		// From Detector
		llrp1[i] = received[i];					// ln{P(1)/P(-1)}
    }

	iter = 0;                  // Counter of iterations
	do
	{

		// For every (m,l) such that there is a link between parents and
		// children, qm0[i][j] and qm1[i][j] are initialized to pl[j].
		// Notation: pi (Pearl) = q (MacKay)

		//if( (s1 == 0) && (iter == 0) )
		if( (s2 == 0) && (iter == 0) )
			for (i=0; i<N; i++)                         // run over code nodes
			{
				//printf("received = %10.7lf\n", received[i]);
				for (j=0; j<Bit_node[0][i].Num_chks_bit; j++)       // run over check nodes
				{
					Bit_node[ind][i].pi1[j] = llrp1[i];
					//printf("i, j,  pi1 = %d,%d   %10.7lf\n", i, j, Bit_node[ind][i].pi1[j]);
				}
			}
		else
		{

		// ------------------------------------
		//         ***** STEP 1 *****
		// VERTICAL STEP = TOP-DOWN PROPAGATION
		// ------------------------------------
		//
		// MacKay:
		// Take the computed values of rm0, rm1 and update the values of
		// the probabilities qm0, qm1
		//
		// Pearl:
		// Each node u_l computes new "pi" messages to be send to its
		// children x_1, x_2, ..., x_J

			for (i=0; i<N; i++)
				for (j=0; j<Bit_node[0][i].Num_chks_bit; j++)
				{
					Bit_node[ind][i].pi1[j] = 0.0;
					for (l=0; l<Bit_node[0][i].Num_chks_bit; l++)
					{
						aux = Bit_node[0][i].index[l] - 1; 
						if ( aux != (Bit_node[0][i].index[j]-1) )
						{
							// Compute index "m" of message from children
							m = 0;
							while (  ( (CHK_node[0][aux].index[m]-1) != i ) && ( m < CHK_node[0][aux].Num_bits_chk )  ) 
								m++;

							Bit_node[ind][i].pi1[j] += CHK_node[ind][aux].lambda1[m];
						}
					}

					Bit_node[ind][i].pi1[j] += llrp1[i];
					//printf("---->  i, j,  pi1 = %d,%d    %10.7lf\n", i, j, Bit_node[ind][i].pi1[j]);

					if (Bit_node[ind][i].pi1[j] < -50.0)
						Bit_node[ind][i].pi1[j] = -50.0;
					if (Bit_node[ind][i].pi1[j] > 50.0)
						Bit_node[ind][i].pi1[j] = 50.0;
				}
		}

		// ---------------------------------------
		//         ***** STEP 2 *****
		// HORIZONTAL STEP = BOTTOM-UP PROPAGATION
		// ---------------------------------------
		//
		// MacKay:
		// Run through the checks m and compute, for each n in N(m) the
		// probabilitiy of a check symbol when code symbol is 0 (or 1)
		// given that the other code symbols have distribution qm0, qm1
		//
		// Pearl:
		// Node x_m computes new "lambda" messages to be sent to its parents
		// u_1, u_2, ..., u_K

		for (i=0; i<M; i++)
			for (j=0; j<CHK_node[0][i].Num_bits_chk; j++)
			{	
				delt = 0.0;
				sign = 0;                           // Keep track of sign of delt, 0: positive, 1: negative

				for (l=0; l<CHK_node[0][i].Num_bits_chk; l++)
				{
					aux = CHK_node[0][i].index[l];
					if (aux != CHK_node[0][i].index[j])
					{
						// --------------------------------------------------------
						//  Compute the index "m" of the message from parent node
						// --------------------------------------------------------
						m = 0;
						while (  ( (Bit_node[0][aux-1].index[m]-1) != i ) && ( m < Bit_node[0][aux-1].Num_chks_bit)  ) 
							m++;

						if (Bit_node[ind][aux-1].pi1[m] < 0.0) 
							sign ^= 1;
					
						delt += F(fabs(Bit_node[ind][aux-1].pi1[m]));
						//printf("pi1, delt =  %lf, %lf \n", Bit_node[ind][aux-1].pi1[m], delt);
					}
				}
      
				if (sign == 0)
					CHK_node[ind][i].lambda1[j] = F(delt);
				else
					CHK_node[ind][i].lambda1[j] = -F(delt);

				// Prevent numerical underflow
				if (CHK_node[ind][i].lambda1[j] < -50.0)
					CHK_node[ind][i].lambda1[j] = -50.0;
				if (CHK_node[ind][i].lambda1[j] > 50.0)
					CHK_node[ind][i].lambda1[j] = 50.0;

				//printf("i, j, lambda1 = %d,%d  %10.7lf\n", i, j, CHK_node[ind][i].lambda1[j]);
			}
		
		// Increment the number of iterations, and check if maximum reached
		iter++;
	} while (iter < max_iter);

	// DECODING:
	// MacKay: At this step we also compute the (unconditional) pseudo-
	// posterior probalilities "q0, q1" to make tentative decisions

	for (i=0; i<N; i++)
    {
		q1[i] = 0.0;
		for (j=0; j<Bit_node[0][i].Num_chks_bit; j++)
		{
			aux = Bit_node[0][i].index[j] - 1; 

			// Compute index "m" of message from children
			m = 0;
			while (  ( (CHK_node[0][aux].index[m]-1) != i ) && ( m < CHK_node[0][aux].Num_bits_chk )  ) 
				m++;

			q1[i] += CHK_node[ind][aux].lambda1[m];
		}

		if (q1[i] < -50.0)
			q1[i] = -50.0;
		if (q1[i] > 50.0)
			q1[i] = 50.0;

		// To the MIMO Detector
		decoded[i] = q1[i];
	}

	return 0;
}

long double F(long double x)
{
	long double interm;
	
	if (x == 0.0) 
		return(1.0e+30);
	if (fabs(x) > 50.0)
		interm = 0.0;
	else
		interm = logl ( (expl(x)+1.0)/(expl(x)-1.0) );
	
	return(interm);
}

void LDPC_Encode(int *data, int *codeword)
// Systematic encoding 
{
	int i, j;
	
	for (j=0; j<Cn; j++)
	{
		if (j>=Cn-Ck)                     // information bits
			codeword[j] = data[j-(Cn-Ck)];
		else								// parity bits
		{
			codeword[j] = 0;
			for (i=0; i<Ck; i++)
				codeword[j] ^= ( data[i] * G[i][j] ) & 0x01;
		}
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
