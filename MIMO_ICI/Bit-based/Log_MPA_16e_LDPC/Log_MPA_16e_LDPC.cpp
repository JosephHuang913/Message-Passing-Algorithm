/*========================================================*/
/* Author: Chao-wang Huang                                                                       */
/* Date: Thursday, December 16, 2010                                                       */
/* A Bit-based Message Passing Algorithm for MIMO Channel is simulated */
/* 4 by 4 QPSK                                                                                           */
/* MIMO Detector: Log-domain MPA with Sum-Product Rule                       */
/* LDPC code: 802.16e LDPC code                                                           */
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
void Encoder_H_Generator(int, int, int **, int **);
void MatrixXORMultiply (int **, int **, int **, int, int, int);
void SPA_decoder(int, int, int, long double *, long double *, long double, int);
long double Phi(long double);

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 	3.14159265358979
#define Num_packet 500					// number of packets simulated
#define Num_codeword 1 				// number of codewords per packet
//const int Mc = 4;								// modulation order (4 for 16QAM)
const int Mc = 2;								// modulation order (2 for QPSK)
const int T = 1;									// Number of transmit antenna
const int R = 1;									// Number of receive antenna
const int Iterate2 = 1;
const int Iterate1 = 1;
const int Iterate3 = 1;
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
int s1;
long double ***ts;
long double *L_Q, ***L_q, ***L_r, *Lpi;
int **mark_v, **mark_h, *wc_per_column,*wr_per_row;

int main(void)
{
	time_t  ti, start, end;
	int p, i, j, q, l/*, s1*/, s2/*, s3*/, x, z, t, d, *data, **data_bit, **Ak, *S_random, Cn;
	int **err_count, *delay, **coded_bit, flag, *interleaved, *inter, *coded, *bit, *de_inter;
	int CodeType, z_factor, max_itr, k, m, UncBitSize, v, phi, count;
	int **uT, **p1T, **AU, **Bp1, **AUBp1, **p2T, H_row, H_column;
	int **inv_T, **H, **A, **B, **C, **D, **E, **ET, **ETA, **ETAC, **inv_phi;
	long double r, Eb_No, snr, noise_pwr, noise[2], Y[2], ***Yk, **err_rate, Ps;
	long double sum, pdf, ***L_r_pn, ***L_q_np, **sym_I, **sym_Q, *I, *Q, *gain, ****fade, ***H_f, *fft;
	long double L_sum, L_tmp[2], *L_a_u, *L_a_d, *L_in, *L_out, *L_o, snr_rms, sigma_N;
	FILE *ber, *records, *s_inter;

	// CodeType: 1/2=>1, 2/3A=>2, 2/3B=>3, 3/4A=>4, 3/4B=>5, 5/6=>6
	// Z factor: 24, 28, 32, 36, 40, 44,  48,  52,  56,  60,  64,  68,  72,  76,  80,  84,  88,  92,  96
	// n bits:  576,672,768,864,960,1056,1152,1248,1344,1440,1536,1632,1728,1824,1920,2016,2112,2208,2304

	CodeType = 1;
	z_factor = 84;
	max_itr = 50;
	Cn = z_factor * 24;

	switch(CodeType)
    {
		case 1 : k = 12; m = 12; r = 1/2.0; break;
		case 2 : k = 16; m = 8; r = 2/3.0; break;
		case 3 : k = 16; m = 8; r = 2/3.0; break;
		case 4 : k = 18; m = 6; r = 3/4.0; break;
		case 5 : k = 18; m = 6; r = 3/4.0; break;
		case 6 : k = 20; m = 4; r = 5/6.0; break;
		default:   cout << "Undefined CodeType\n";  break;
    }

	UncBitSize = (int)(z_factor * 24 * r);
	H_row = m * z_factor;
	H_column = 24 * z_factor;
	phi = (int)(80 * z_factor / 96.0);

	start = time(NULL);
	printf("BER Performance of Bit-based Message Passing in MIMO Channel\n");
	printf("802.16e LDPC code\n");
	printf("Number of transmit antenna is %d\n", T);
	printf("Number of receive antenna is %d\n", R);
	printf("Vehicular speed: %f km/h\n", vc);
	printf("Maximum Doppler frequency: %f\n", Doppler);
	printf("Channel fading rate: %f\n", fdxT);
	printf("Maximum number of bits of simulation = %d\n\n", Cn*Num_codeword*Num_packet);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Record_Bit-based_Message_Passing.log", "a");
	fprintf(records, "BER Performance of Bit-based Message Passing in MIMO Channel\n");
	fprintf(records, "802.16e LDPC code\n");
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
		data_bit[i] = new int[UncBitSize];
	data = new int[Mc*T];
	coded_bit = new int*[Num_codeword];
	for(i=0; i<Num_codeword; i++)
		coded_bit[i] = new int[Cn];
	interleaved = new int[N_fft*T*Mc];
	inter = new int[N_fft*T*Mc];
	coded = new int[Cn];
	bit = new int[UncBitSize];
	S_random = new int[N_fft*T*Mc];
	de_inter = new int[N_fft*T*Mc];
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
	err_count = new int*[Iterate2];
	for(i=0; i<Iterate2; i++)
		err_count[i] = new int[Iterate1];	
	err_rate = new long double*[Iterate2];
	for(i=0; i<Iterate2; i++)
		err_rate[i] = new long double[Iterate1];
	I = new long double[T];
	Q = new long double[T];
	L_a_u = new long double[N_fft*T*Mc];
	L_a_d = new long double[N_fft*T*Mc];
	L_o = new long double[N_fft*T*Mc];
	L_in = new long double[Cn];
	L_out = new long double[Cn];
	H_f = new long double**[R];
	for(i=0; i<R; i++)
   		H_f[i] = new long double*[T];
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			H_f[i][j] = new long double[2*N_fft+1];
	fft = new long double[2*N_fft+1];
	H = new int*[H_row];
	for (i=0; i<H_row; i++)
		H[i] = new int[H_column];
   	inv_T = new int*[H_row - z_factor];
	for (i=0; i<(H_row - z_factor); i++)
		inv_T[i] = new int[H_row - z_factor];
	A = new int*[H_row - z_factor];
	for (i=0; i<(H_row - z_factor); i++)
		A[i] = new int[H_column - H_row];
   	B = new int*[H_row - z_factor];
	for (i=0; i<(H_row - z_factor); i++)
		B[i] = new int[z_factor];
	C = new int*[z_factor];
	for (i=0; i<z_factor; i++)
		C[i] = new int[H_column - H_row];
   	D = new int*[z_factor];
	for (i=0; i<z_factor; i++)
		D[i] = new int[z_factor];
   	E = new int*[z_factor];
	for (i=0; i<(H_row - z_factor); i++)
		E[i] = new int[H_row - z_factor];
	ET = new int*[z_factor];
	for (i=0; i<z_factor; i++)
		ET[i] = new int[H_row - z_factor];
	ETA = new int*[z_factor];
	for (i=0; i<z_factor; i++)
		ETA[i] = new int[H_column - H_row];
	ETAC = new int*[z_factor];
	for (i=0; i<z_factor; i++)
		ETAC[i] = new int[H_column - H_row];
	inv_phi = new int*[z_factor];
	for (i=0; i<z_factor; i++)
		inv_phi[i] = new int[z_factor];
	uT = new int*[UncBitSize];
	for (i=0; i<UncBitSize; i++)
		uT[i] = new int[1];
	p1T = new int*[z_factor];
	for (i=0; i<z_factor; i++)
		p1T[i] = new int[1];
	AU = new int*[H_row - z_factor];
	for (i=0; i<(H_row - z_factor); i++)
		AU[i] = new int[1];
	Bp1 = new int*[H_row - z_factor];
	for (i=0; i<(H_row - z_factor); i++)
		Bp1[i] = new int[1];
	AUBp1 = new int*[H_row - z_factor];
	for (i=0; i<(H_row - z_factor); i++)
		AUBp1[i] = new int[1];
	p2T = new int*[H_row - z_factor];
	for (i=0; i<(H_row - z_factor); i++)
		p2T[i] = new int[1];
	L_Q = new long double[H_column];
	L_q = new long double**[Num_codeword];
	for(i=0; i<Num_codeword; i++)
		L_q[i] = new long double*[H_column];
	for(i=0; i<Num_codeword; i++)
		for(j=0; j<H_column; j++)
			L_q[i][j] = new long double[H_row];
	L_r = new long double**[Num_codeword];
	for(i=0; i<Num_codeword; i++)
		L_r[i] = new long double*[H_row];
	for(i=0; i<Num_codeword; i++)
		for(j=0; j<H_row; j++)
			L_r[i][j] = new long double[H_column];
	//Lpi = new long double[H_column];
	Lpi = new long double[2*N_fft];
	mark_v = new int*[H_column];
	for(i=0; i<H_column; i++)
		mark_v[i] = new int[6];			// maximum wc is 6
	mark_h = new int*[H_row];
	for(i=0; i<H_row; i++)
		mark_h[i] = new int[24];		// maximum wr is 20
	wc_per_column = new int[H_column];
	wr_per_row = new int[H_row];

	// Channel Weighting Gain and Channel path delay
	for(i=0; i<path_num; i++)
	{
		delay[i] = path_delay[i];
		gain[i] = sqrtl(path_gain[i]);
	}

	srand((unsigned) time(&ti));

	// Generate H matrix
	Encoder_H_Generator(CodeType, z_factor, H, inv_T);

	// Generate R&U algorithm submatrix
	for (i=0; i<(m-1)*z_factor; i++)
		for (j=0; j<k*z_factor; j++)
			A[i][j] = H[i][j];

	for (i=0; i<(m-1)*z_factor; i++)
		for (j=0; j<z_factor; j++)
			B[i][j] = H[i][j + k*z_factor];

	for (i=0; i<z_factor; i++)
		for (j=0; j<k*z_factor; j++)
			C[i][j] = H[i + (m-1)*z_factor][j];

	for (i=0; i<z_factor; i++)
		for (j=0; j<z_factor; j++)
			D[i][j] = H[i + (m-1)*z_factor][j + k*z_factor];
	
	for (i=0; i<z_factor; i++)
		for (j=0; j<(24-k-1)*z_factor; j++)
			E[i][j] = H[i + (m-1)*z_factor][j + (k+1)*z_factor];

	// ET=E*inv_T
	MatrixXORMultiply(E, inv_T, ET, z_factor, (H_row-z_factor), (H_row-z_factor));

	// ETA=ET*A
	MatrixXORMultiply(ET, A, ETA, z_factor, (H_row-z_factor), (H_column-H_row));

	// ETAC=ETA+C
	for (i=0; i<z_factor; i++)
		for (j=0; j<(H_column - H_row); j++)
			ETAC[i][j] = (int)((ETA[i][j] + C[i][j])%2);

	// coding_rate=3/4B  (E*inv_T*B+D != I), Adjust ETAC
	if  (CodeType == 5)
	{
		for (i=0; i<z_factor; i++)
			for (j=0; j<z_factor; j++)
				if ((int)((i+z_factor-phi)%z_factor) == j)
					inv_phi[i][j] = 1;				
				else
					inv_phi[i][j] = 0;

		MatrixXORMultiply(inv_phi, ETAC, ETAC, z_factor, z_factor, (H_column-H_row));
	}

	for(i=0; i<H_column; i++)
		for(j=0; j<6; j++)
			mark_v[i][j] = -1;
	
	for(i=0; i<H_row; i++)
		for(j=0; j<24; j++)
			mark_h[i][j] = -1;

	for(i=0; i<H_column; i++)
	{
		count = 0;
		for(j=0; j<H_row; j++) 
			if(H[j][i] == 1)
			{
				mark_v[i][count] = j;
				count++;
			}

		wc_per_column[i] = count;
	}

	for(j=0; j<H_row; j++)
	{
		count = 0;
		for(i=0; i<H_column; i++) 
			if(H[j][i] == 1)
			{	
				mark_h[j][count] = i;
				count++;
			}

		wr_per_row[j] = count;
	}

/*======================================*/
/* S-random interleaver and de-interleaver (S=64) */
/*======================================*/
	fopen_s(&s_inter, "s_random.log", "r");
	for(i=0; i<N_fft*T*Mc; i++)		// Interleaver
   		fscanf_s(s_inter, "%d %d", &de_inter[i], &S_random[i]);
	for(i=0; i<N_fft*T*Mc; i++)		// De-interleaver
   		de_inter[S_random[i]] = i;
	fclose(s_inter);

/*======================*/
/* main simulation loop            */
/*======================*/
	fopen_s(&ber, "BER_MPA_16e_LDPC.log", "w");
	for(snr=1; snr<=4; snr+=1)
	{
   		for(s2=0; s2<Iterate2; s2++)
			for(s1=0; s1<Iterate1; s1++)
   				err_count[s2][s1] = 0;

   		// noise power calculation
		Eb_No = (double)snr;
		Ps = 1 * T;
		//noise_pwr = 0.5*Ps*R/((float)N_fft*T*((double)UncBitSize/(double)Cn)*powl(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Time Domain)
		noise_pwr = 0.5*Ps*R/((float)T*r*powl(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption (Freq. Domain)
		//noise_pwr = 0.25*Ps*R/((float)T*((double)UncBitSize/(double)Cn)*powl(10.0, Eb_No/10.0));	// 16QAM, Nyquist filter assumption (Freq. Domain)
		printf("Eb_No = %f\n", Eb_No);
		snr_rms = 2.0 * sqrtl(1.0*r*(powl(10.0, (snr/10.0))));
		sigma_N = powl(0.5/powl(10, (double)snr/10.0), 0.5);

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
				for(j=0; j<UncBitSize; j++)
					if(rand()/(double)RAND_MAX >= 0.5)
					{
						bit[j] = 1;
						data_bit[i][j] = bit[j];
					}
					else
					{
						bit[j] = 0;
						data_bit[i][j] = bit[j];
					}

				// Generate coded bits [uncoded_bit p1 p2]
				// uT=u'
				for (j=0; j<UncBitSize; j++)
					uT[j][0] = data_bit[i][j];

				// p1'=((E*inv(T)*A+C)*u');
				MatrixXORMultiply(ETAC, uT, p1T, z_factor, UncBitSize, 1);
				MatrixXORMultiply(A, uT, AU, H_row-z_factor, UncBitSize, 1);
				MatrixXORMultiply(B, p1T, Bp1, H_row-z_factor, z_factor, 1);

				for (j=0; j<(H_row-z_factor); j++)
					AUBp1[j][0] = (int)((AU[j][0] + Bp1[j][0])%2);

				// p2=(inv(T)*(A*u'+B*p1'))';
				MatrixXORMultiply(inv_T, AUBp1, p2T, (H_row-z_factor), (H_row-z_factor), 1);

				//coded_bit=[coded_bit u p1' p2'];
				v = (int)(z_factor*24*r);

				for (j=0; j<z_factor*24; j++)
				{
					if (j < v)
						coded[j] = data_bit[i][j];
					else if( (j>=v) & (j<(v+z_factor)) )
						coded[j] = p1T[j - v][0];
					else
						coded[j] = p2T[j - (v+z_factor)][0];
				}

				for(l=0; l<z_factor*24; l++)
				{
					coded_bit[i][l] = coded[l];
					inter[i*z_factor*24 + l] = coded[l];
				}
			}

			// Zero padding
			for(l=Num_codeword*z_factor*24; l<N_fft*T*Mc; l++)
				inter[l] = 0;

			// Interleaving: (8192, 64) S-random interleaver
			for(i=0; i<N_fft*T*Mc; i++)
      			//interleaved[S_random[i]] = inter[i];
				interleaved[i] = inter[i];		// No interleaving

			for(l=0; l<N_fft; l++)
   				for(i=0; i<T; i++)
	   			{
					// QPSK Mapping
					//sym_I[i][l] = (2*interleaved[l*Mc*T+i*Mc]-1)/sqrtl(2.0);
         			//sym_Q[i][l] = (2*interleaved[l*Mc*T+i*Mc+1]-1)/sqrtl(2.0);
					sym_I[i][l] = (-2*interleaved[l*Mc*T+i*Mc]+1)/sqrtl(2.0);
         			sym_Q[i][l] = (-2*interleaved[l*Mc*T+i*Mc+1]+1)/sqrtl(2.0);
				}				
/*
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
*/			
         	//AWGN channel
         	for(l=0; l<N_fft; l++)
				for(i=0; i<R; i++)
				{
   					AWGN_noise(0, noise_pwr, &noise[0]);
	      	   		//Yk[i][l][0] += noise[0];
		     		//Yk[i][l][1] += noise[1];
					Yk[i][l][0] = sym_I[i][l] + noise[0];
         			Yk[i][l][1] = sym_Q[i][l] + noise[1];
					//Lpi[2*l] = Yk[i][l][0] * snr_rms * sqrtl(2.0);
					//Lpi[2*l+1] = Yk[i][l][1] * snr_rms * sqrtl(2.0);
					Lpi[2*l] = 2.0 * Yk[i][l][0] / noise_pwr;
					Lpi[2*l+1] = 2.0 * Yk[i][l][1] / noise_pwr;
				}

/*========================================================*/
/*  Log-Domain Bit-based Message Passing Algorithm & Progressive PIC */
/* Channel Decoder: 802.16e LDPC code                                                  */
/*========================================================*/
/*			
			// MIMO Detector: LLR Initialization
			for(t=0; t<N_fft; t++)
				for(j=0; j<T; j++)
					for(m=0; m<Mc; m++)
						for(i=0; i<R; i++)
							L_q_np[t][j*Mc+m][i] = 0.0;
	
			for(j=0; j<H_column; j++)
				for(l=0; l<H_row; l++)
					L_q[i][j][l] = 0.0;

			for(j=0; j<H_row; j++)
				for(l=0; l<H_column; l++)
					L_r[i][j][l] = 0.0;
*/
//			for(s2=0; s2<Iterate2; s2++)
//			{
//				for(s1=0; s1<Iterate1; s1++)
//				{
				s2 = s1 = 0;
/*					for(t=0; t<N_fft; t++)
					{
						// Extrinsic information of function nodes
			    		for(i=0; i<R; i++)
							for(j=0; j<T; j++)
								for(m=0; m<Mc; m++)
           						{
									flag = 0;
									for(z=Num_codeword*Cn; z<N_fft*Mc*T; z++)
										if(t*Mc*T+j*Mc+m == S_random[z])
										{
											L_r_pn[t][i][j*Mc+m] = -50.0;
											flag = 1;
										}
									
									//if(t*Mc*T+j*Mc+m >= Num_codeword*z_factor*24)
									//	L_r_pn[t][i][j*Mc+m] = -50.0;
									//else
									if(flag == 0)
									{
										for(q=0; q<2; q++)		// Probability of 1 and 0
										{
               								L_tmp[q] = 0.0;
		           							for(x=0; x<powl(2.0,Mc*T-1.0); x++)
	   										{
		                  						l = 0;
				   								for(z=0; z<Mc*T; z++)
					   		         				if(z != j*Mc+m)
       												{
								  			  			data[z] = (x>>l) & 1;
							 							l++;
														for(d=Num_codeword*Cn; d<N_fft*Mc*T; d++)
															if(t*Mc*T+z == S_random[d])
														//if( t*Mc*T+z >= Num_codeword*z_factor*24)
																data[z] = 0;
													}

												data[j*Mc+m] = q;
			   									for(z=0; z<T; z++)
				  								{
													// QPSK Mapping
													I[z] = (2*data[z*Mc]-1)/sqrtl(2.0);
         											Q[z] = (2*data[z*Mc+1]-1)/sqrtl(2.0);
			       								}
	
					         					Y[0] = Y[1] = 0.0;
												for(z=0; z<T; z++)
												{
													Y[0] += (I[z] * H_f[i][z][2*t] - Q[z] * H_f[i][z][2*t+1]);
               										Y[1] += (I[z] * H_f[i][z][2*t+1] + Q[z] * H_f[i][z][2*t]);
												}
										
								 				L_sum = 0.0;
												if(s2 != 0 || s1 != 0)
	                     							for(z=0; z<T; z++)
														for(d=0; d<Mc; d++)
			             									if( (z*Mc+d != j*Mc+m) && (data[z*Mc+d] == 1) )
				       											L_sum += L_q_np[t][z*Mc+d][i];
																
												if (L_sum < -50.0)
													L_sum = -50.0;
												if (L_sum > 50.0)
													L_sum = 50.0;

												sum = 0.0;
												for(z=0; z<2; z++)
							        				sum += powl(Yk[i][t][z] - Y[z], 2.0);

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
*/				
					// LDPC Decoding
					for(i=0; i<Num_codeword; i++)
					{
						//for(l=0; l<Cn; l++)
						//	L_in[l] = L_a_u[i*Cn+l];

						//SPA_decoder(max_itr, H_column, H_row, L_in, L_out, snr_rms, i);
						SPA_decoder(max_itr, H_column, H_row, Lpi, L_out, sigma_N, i);

						for(l=0; l<Cn; l++)
							L_o[i*Cn+l] = L_out[l];

						// Hard Decision
						for(l=0; l<Cn; l++)
						{
					 		//if( (L_o[i*Cn+l]+L_a_u[i*Cn+l]) >= 0.0)
							if( L_o[i*Cn+l] < 0.0)
								Ak[i][l] = 1;
       						else
       							Ak[i][l] = 0;

	            			//err_count[s2][s1] += Error_count(coded_bit[i][l], Ak[i][l]);
							err_count[s2][s1] += Error_count(interleaved[i*Cn+l], Ak[i][l]);
						}
					}
/*
					for(l=Num_codeword*Cn; l<N_fft*T*Mc; l++)
						L_o[l] = -50.0;

					// Interleaving
					for(i=0; i<N_fft*T*Mc; i++)
      					L_a_d[S_random[i]] = L_o[i];
						//L_a_d[i] = L_o[i];		// No interleaving

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

								//if(t*Mc*T+j*Mc+m >= Num_codeword*Cn)
								//	L_q_np[t][j*Mc+m][l] = -50.0;
								//else
								if(flag == 0)
									for(i=0; i<R; i++)
									{
										L_q_np[t][j*Mc+m][i] = 0.0;
				   						for(z=0; z<R; z++)
		       								if(z != i)
												L_q_np[t][j*Mc+m][i] += L_r_pn[t][z][j*Mc+m];
	
										L_q_np[t][j*Mc+m][i] += L_a_d[t*Mc*T+j*Mc+m];
										
										if (L_q_np[t][j*Mc+m][i] < -50.0)
											L_q_np[t][j*Mc+m][i] = -50.0;
										if (L_q_np[t][j*Mc+m][i] > 50.0)
											L_q_np[t][j*Mc+m][i] = 50.0;
				 					}
							}
				}
			}
*/
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
	for(i=0; i<Iterate2; i++)
		delete err_count[i];
	delete err_count;
	for(i=0; i<Iterate2; i++)
		delete err_rate[i];
	delete err_rate;
	delete I;
	delete Q;
	for(i=0; i<H_row; i++)
		delete H[i];
	delete H;
	delete L_a_u;
	delete L_a_d;
	delete L_o;
	delete L_in;
	delete L_out;
	for(i=0; i<R; i++)
		for(j=0; j<T; j++)
			delete H_f[i][j];
	for(i=0; i<R; i++)
		   delete H_f[i];
	delete H_f;
	delete fft;
	delete interleaved;
	delete inter;
	delete coded;
	delete bit;
	delete S_random;
	delete de_inter;
	for (i=0;i<H_row;i++)
		delete H[i];
	delete H;
	for (i=0;i<(H_row-z_factor);i++)
		delete inv_T[i];
   	delete inv_T;
	for (i=0;i<(H_row-z_factor);i++)
		delete A[i];
	delete A;
	for (i=0;i<(H_row-z_factor);i++)
		delete B[i];
   	delete B;
	for (i=0;i<z_factor;i++)
		delete C[i];
	delete C;
	for (i=0;i<z_factor;i++)
		delete D[i];
   	delete D;
	for (i=0;i<(H_row-z_factor);i++)
		delete E[i];
   	delete E;
	for (i=0;i<z_factor;i++)
		delete ET[i];
	delete ET;
	for (i=0;i<z_factor;i++)
		delete ETA[i];
	delete ETA;
	for (i=0;i<z_factor;i++)
		delete ETAC[i];
	delete ETAC;
	for (i=0;i<z_factor;i++)
		delete inv_phi[i];
	delete inv_phi;
	delete Yk;
	for (i=0;i<UncBitSize;i++)
		delete uT[i];
	delete uT;
	for (i=0;i<z_factor;i++)
		delete p1T[i];
	delete p1T;
	for (i=0;i<(H_row-z_factor);i++)
		delete AU[i];
	delete AU;
	for (i=0;i<(H_row-z_factor);i++)
		delete Bp1[i];
	delete Bp1;
	for (i=0;i<(H_row-z_factor);i++)
		delete AUBp1[i];
	delete AUBp1;
	for (i=0;i<(H_row-z_factor);i++)
		delete p2T[i];
	delete p2T;
	delete L_Q;
	for(i=0; i<Num_codeword; i++)
		for(j=0; j<H_column; j++)
			delete L_q[i][j];
	for(i=0; i<Num_codeword; i++)
		delete L_q[i];
	delete L_q;
	for(i=0; i<Num_codeword; i++)
		for(j=0; j<H_row; j++)
			delete L_r[i][j];
	for(i=0; i<Num_codeword; i++)
		delete L_r[i];
	delete L_r;
	delete Lpi;
	for(i=0; i<H_column; i++)
		delete mark_v[i];
	delete mark_v;
	for(i=0; i<H_row; i++)
		delete mark_h[i];
	delete mark_h;
	delete wc_per_column;
	delete wr_per_row;

	end = time(NULL);
	printf("Total elapsed time: %.0f(sec)\n", difftime(end, start));
	fprintf(records, "Total elapsed time: %.0f(sec)\n\n", difftime(end, start));
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

void Encoder_H_Generator(int CodeType, int z_factor, int **H, int **inv_T)
{
	int i, j, k, m, HbmRow, HbmCol;
	int Hbm[12][24];

	// CodeType=1/2
	int A1[12][24] = {{-1, 94, 73, -1, -1, -1, -1, -1, 55, 83, -1, -1,  7,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
					   {-1, 27, -1, -1, -1, 22, 79,  9, -1, -1, -1, 12, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1},
					   {-1, -1, -1, 24, 22, 81, -1, 33, -1, -1, -1,  0, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1},
					   {61, -1, 47, -1, -1, -1, -1, -1, 65, 25, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1, -1},
					   {-1, -1, 39, -1, -1, -1, 84, -1, -1, 41, 72, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1, -1},
					   {-1, -1, -1, -1, 46, 40, -1, 82, -1, -1, -1, 79,  0, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1, -1},
					   {-1, -1, 95, 53, -1, -1, -1, -1, -1, 14, 18, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1, -1},
					   {-1, 11, 73, -1, -1, -1,  2, -1, -1, 47, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1},
					   {12, -1, -1, -1, 83, 24, -1, 43, -1, -1, -1, 51, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1},
					   {-1, -1, -1, -1, -1, 94, -1, 59, -1, -1, 70, 72, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1},
					   {-1, -1,  7, 65, -1, -1, -1, -1, 39, 49, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0},
					   {43, -1, -1, -1, -1, 66, -1, 41, -1, -1, -1, 26,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0}};

	// CodeType=2/3A
	int A2[8][24]	= {{ 3,  0, -1, -1,  2,  0, -1,  3,  7, -1,  1,  1, -1, -1, -1, -1,  1,  0, -1, -1, -1, -1, -1, -1},
					   {-1, -1,  1, -1, 36, -1, -1, 34, 10, -1, -1, 18,  2, -1,  3,  0, -1,  0,  0, -1, -1, -1, -1, -1},
					   {-1, -1, 12,  2, -1, 15, -1, 40, -1,  3, -1, 15, -1,  2, 13, -1, -1, -1,  0,  0, -1, -1, -1, -1},
					   {-1, -1, 19, 24, -1,  3,  0, -1,  6, -1, 17, -1, -1, -1,  8, 39, -1, -1, -1,  0,  0, -1, -1, -1},
					   {20, -1,  6, -1, -1, 10, 29, -1, -1, 28, -1, 14, -1, 38, -1, -1,  0, -1, -1, -1,  0,  0, -1, -1},
					   {-1, -1, 10, -1, 28, 20, -1, -1,  8, -1, 36, -1,  9, -1, 21, 45, -1, -1, -1, -1, -1,  0,  0, -1},
					   {35, 25, -1, 37, -1, 21, -1, -1,  5, -1, -1,  0, -1,  4, 20, -1, -1, -1, -1, -1, -1, -1,  0,  0},
					   {-1,  6,  6, -1, -1, -1,  4, -1, 14, 30, -1,  3, 36, -1, 14, -1,  1, -1, -1, -1, -1, -1, -1,  0}};

	// CodeType=2/3B
	int A3[8][24]	= {{ 2, -1, 19, -1, 47, -1, 48, -1, 36, -1, 82, -1, 47, -1, 15, -1, 95,  0, -1, -1, -1, -1, -1, -1},
					   {-1, 69, -1, 88, -1, 33, -1,  3, -1, 16, -1, 37, -1, 40, -1, 48, -1,  0,  0, -1, -1, -1, -1, -1},
					   {10, -1, 86, -1, 62, -1, 28, -1, 85, -1, 16, -1, 34, -1, 73, -1, -1, -1,  0,  0, -1, -1, -1, -1},
		   			   {-1, 28, -1, 32, -1, 81, -1, 27, -1, 88, -1,  5, -1, 56, -1, 37, -1, -1, -1,  0,  0, -1, -1, -1},
					   {23, -1, 29, -1, 15, -1, 30, -1, 66, -1, 24, -1, 50, -1, 62, -1, -1, -1, -1, -1,  0,  0, -1, -1},
					   {-1, 30, -1, 65, -1, 54, -1, 14, -1,  0, -1, 30, -1, 74, -1,  0, -1, -1, -1, -1, -1,  0,  0, -1},
					   {32, -1,  0, -1, 15, -1, 56, -1, 85, -1,  5, -1,  6, -1, 52, -1,  0, -1, -1, -1, -1, -1,  0,  0},
					   {-1,  0, -1, 47, -1, 13, -1, 61, -1, 84, -1, 55, -1, 78, -1, 41, 95, -1, -1, -1, -1, -1, -1,  0}};

	// CodeType=3/4A
	int A4[6][24]	= {{ 6, 38,  3, 93, -1, -1, -1, 30, 70, -1, 86, -1, 37, 38,  4, 11, -1, 46, 48,  0, -1, -1, -1, -1},
					   {62, 94, 19, 84, -1, 92, 78, -1, 15, -1, -1, 92, -1, 45, 24, 32, 30, -1, -1,  0,  0, -1, -1, -1},
					   {71, -1, 55, -1, 12, 66, 45, 79, -1, 78, -1, -1, 10, -1, 22, 55, 70, 82, -1, -1,  0,  0, -1, -1},
					   {38, 61, -1, 66,  9, 73, 47, 64, -1, 39, 61, 43, -1, -1, -1, -1, 95, 32,  0, -1, -1,  0,  0, -1},
					   {-1, -1, -1, -1, 32, 52, 55, 80, 95, 22,  6, 51, 24, 90, 44, 20, -1, -1, -1, -1, -1, -1,  0,  0},
					   {-1, 63, 31, 88, 20, -1, -1, -1,  6, 40, 56, 16, 71, 53, -1, -1, 27, 26, 48, -1, -1, -1, -1,  0}};

	// CodeType=3/4B
	int A5[6][24]	= {{-1, 81, -1, 28, -1, -1, 14, 25, 17, -1, -1, 85, 29, 52, 78, 95, 22, 92,  0,  0, -1, -1, -1, -1},
					   {42, -1, 14, 68, 32, -1, -1, -1, -1, 70, 43, 11, 36, 40, 33, 57, 38, 24, -1,  0,  0, -1, -1, -1},
					   {-1, -1, 20, -1, -1, 63, 39, -1, 70, 67, -1, 38,  4, 72, 47, 29, 60,  5, 80, -1,  0,  0, -1, -1},
					   {64,  2, -1, -1, 63, -1, -1,  3, 51, -1, 81, 15, 94,  9, 85, 36, 14, 19, -1, -1, -1,  0,  0, -1},
					   {-1, 53, 60, 80, -1, 26, 75, -1, -1, -1, -1, 86, 77,  1,  3, 72, 60, 25, -1, -1, -1, -1,  0,  0},
					   {77, -1, -1, -1, 15, 28, -1, 35, -1, 72, 30, 68, 85, 84, 26, 64, 11, 89,  0, -1, -1, -1, -1,  0}};

	// CodeType=5/6
	int A6[4][24]	= {{ 1, 25, 55, -1, 47,  4, -1, 91, 84,  8, 86, 52, 82, 33,  5,  0, 36, 20,  4, 77, 80,  0, -1, -1},
					   {-1,  6, -1, 36, 40, 47, 12, 79, 47, -1, 41, 21, 12, 71, 14, 72,  0, 44, 49,  0,  0,  0,  0, -1},
					   {51, 81, 83,  4, 67, -1, 21, -1, 31, 24, 91, 61, 81,  9, 86, 78, 60, 88, 67, 15, -1, -1,  0,  0},
					   {68, -1, 50, 15, -1, 36, 13, 10, 11, 20, 53, 90, 29, 92, 57, 30, 84, 92, 11, 66, 80, -1, -1,  0}};

	// Define Hbm matrix
	switch (CodeType)
	{
		case 1:   // CodeType=1/2

			HbmRow = 12;
			HbmCol = 24;
			
			for (i=0;i<HbmRow;i++)
				for (j=0;j<HbmCol;j++)
					Hbm[i][j]=A1[i][j];
			break;

		case 2:   // CodeType=2/3A

			HbmRow = 8;
			HbmCol = 24;
			
			for (i=0;i<HbmRow;i++)
				for (j=0;j<HbmCol;j++)
					Hbm[i][j]=A2[i][j];
			break;

		case 3:   // CodeType=2/3B

			HbmRow = 8;
			HbmCol = 24;
			
			for (i=0;i<HbmRow;i++)
				for (j=0;j<HbmCol;j++)
					Hbm[i][j]=A3[i][j];
			break;

		case 4:   // CodeType=3/4A

			HbmRow = 6;
			HbmCol = 24;
			
			for (i=0;i<HbmRow;i++)
				for (j=0;j<HbmCol;j++)
					Hbm[i][j]=A4[i][j];
			break;

		case 5:   // CodeType=3/4B

			HbmRow = 6;
			HbmCol = 24;
			
			for (i=0;i<HbmRow;i++)
				for (j=0;j<HbmCol;j++)
					Hbm[i][j]=A5[i][j];
			break;

		case 6:   // CodeType=5/6

			HbmRow = 4;
			HbmCol = 24;
			
			for (i=0;i<HbmRow;i++)
				for (j=0;j<HbmCol;j++)
					Hbm[i][j]=A6[i][j];
			break;

		default:
			cout << "Undefined Coding Type!!!" << endl;
			break;
	}

	// Adjust shift parameter
	switch (CodeType)
	{
		case 2:   // coding_rate=2/3A
			for (i=0;i<HbmRow;i++)
				for (j=0;j<HbmCol;j++)
					if (Hbm[i][j]>0)
						Hbm[i][j]=(int)(Hbm[i][j]%z_factor);
			break;

		default:
			for (i=0;i<HbmRow;i++)
				for (j=0;j<HbmCol;j++)
					if (Hbm[i][j]>0)
						Hbm[i][j]=(int)(Hbm[i][j]*z_factor/96);
			break;
	}

	// Create Zero matrix and Identity matrix
	int **hhZeros = new int *[z_factor];
	int **hhEyes = new int *[z_factor];

	for (i=0;i<z_factor;i++)
	{
		hhZeros[i] = new int[z_factor];
		hhEyes[i] = new int[z_factor];
	}

	for (i=0;i<z_factor;i++)
		for (j=0;j<z_factor;j++)
			if (i==j)
			{
				hhEyes[i][j]=1;
				hhZeros[i][j]=0;
			}
			else
			{
				hhEyes[i][j]=0;
				hhZeros[i][j]=0;
			}

	// Generate H according to the Hbm
	for (i=0;i<HbmRow;i++)
		for (j=0;j<HbmCol;j++)
			if(Hbm[i][j]==-1)
				for (k=0;k<z_factor;k++)
					for (m=0;m<z_factor;m++)
						H[i*z_factor+k][j*z_factor+m]=hhZeros[k][m];
			else if(Hbm[i][j]==0)
				for (k=0;k<z_factor;k++)
					for (m=0;m<z_factor;m++)
						H[i*z_factor+k][j*z_factor+m]=hhEyes[k][m];
			else
				for (k=0;k<z_factor;k++)
					for (m=0;m<z_factor;m++)
						H[i*z_factor+k][j*z_factor+m]=hhEyes[k][(int)((m+z_factor-(Hbm[i][j]))%z_factor)];

	// Generate inv_T
	for (i=0;i<(HbmRow-1);i++)
		for (j=0;j<(HbmRow-1);j++)
			if (i>=j)
				for (k=0;k<z_factor;k++)
					for (m=0;m<z_factor;m++)
						inv_T[i*z_factor+k][j*z_factor+m]=hhEyes[k][m];
			else
				for (k=0;k<z_factor;k++)
					for (m=0;m<z_factor;m++)
						inv_T[i*z_factor+k][j*z_factor+m]=hhZeros[k][m];

	for (i=0;i<z_factor;i++)
	{
		delete hhZeros[i];
		delete hhEyes[i];
	}
	delete hhZeros;
	delete hhEyes;
}

void MatrixXORMultiply (int **A ,int **B,int **C,int rA,int cA,int cB)
{
	// Claculate A*B=C
	// rA is the row of matrix A
	// cA is the colum of matrix A  cA=rB
	// cB is the column of matrix B

	int i, j, k, sum;
	
	for (i=0;i<rA;i++)
		for (j=0;j<cB;j++)
		{
			sum = 0;
			for (k=0; k<cA; k++)
			{
				sum = sum + (A[i][k] * B[k][j]);
				sum = (int)sum%2;
			}

			C[i][j] = sum;
		}
}

void SPA_decoder(int max_itr, int column, int row, long double *receive_y, long double *L_Q, long double snr_rms, int idx)
{
	int i, j, k, itr, temp, count;
	double sum, sum_Q, temp_afa, temp_beta, sign;

	/***************  Start Decoding  ***********************/
	//for (i=0; i<column; i++)
	//{
		// From Detector
		//Lpi[i] = receive_y[i];					// ln{P(1)/P(0)}

		//Lpi[i] = 2.0 * receive_y[i] / powl(sigma_N, 2.0);
		//Lpi[i] = receive_y[i] * snr_rms;
	//}

	for( itr=0; itr<max_itr; itr++)
	{	
		if( (s1 == 0) && (itr == 0) )
		{
			for(i=0; i<column; i++)
				for(k=0; k<wc_per_column[i]; k++)
				{	
					temp = mark_v[i][k];
					L_q[idx][i][temp] = Lpi[i];
				}
		}
		else
		{
			for(i=0; i<column; i++)
				for(k=0; k<wc_per_column[i]; k++)
				{
					sum = 0.0;
					for(count=0; count<wc_per_column[i]; count++)
					{
						temp = mark_v[i][count];
						if(count != k)
							sum += L_r[idx][temp][i];
					}

			        L_q[idx][i][temp] =  sum + Lpi[i];
				}
		}

		for(j=0; j<row; j++)                                    
			for(k=0; k<wr_per_row[j]; k++)
			{
				temp_afa = 1.0;
				temp_beta = 0.0;
				for(count=0; count<wr_per_row[j]; count++)
				{
					temp = mark_h[j][count];
		
					if(count != k)
					{
						if(L_q[idx][temp][j] >= 0)	
							sign = 1.0;
						else
							sign = -1.0;		

						temp_afa = temp_afa * sign;
						temp_beta += Phi(sign * L_q[idx][temp][j]);
					}
				}

				temp = mark_h[j][k];
				L_r[idx][j][temp] = temp_afa * Phi(temp_beta);
			}
	}

	for(i=0; i<column; i++)                                     
	{
		sum_Q = 0.0;
		for(k=0; k<wc_per_column[i]; k++)
		{
			temp = mark_v[i][k];
			sum_Q += L_r[idx][temp][i];
		}

		L_Q[i] = sum_Q + Lpi[i];
	}
}

long double Phi(long double beta)
{
	if(beta == 0.0)
		return logl((expl(1E-12)+1)/(expl(1E-12)-1));
	else
		return logl((expl(beta)+1)/(expl(beta)-1));
}
