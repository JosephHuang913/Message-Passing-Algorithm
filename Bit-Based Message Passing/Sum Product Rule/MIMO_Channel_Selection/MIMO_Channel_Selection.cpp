/*========================================================*/
/* Author: Chao-wang Huang                                                                       */
/* Date: Tuesday, September 09, 2008                                                       */
/* A Bit-based Message Passing Algorithm for MIMO Channel is simulated */
/* Receiver Channel Selection                                                                    */
/* Sum-product Rule                                                                                    */
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

#define Pi 	3.14159265358979
#define num_packet 125000			// number of packets simulated
const int Mc = 2;								// modulation order (2 for QPSK)
const int T = 4;									// Number of transmit antenna
const int R = 4;									// Number of receive antenna
const int Iteration = 10;
double H_I[4][4], H_Q[4][4];
int RCS[4][4];

int main(void)
{
	time_t  t, start, end;
	int i, j, d, k, l, m, n, p, r, s, x, z, *data, *data_bit, *Ak, *err_count, min1[12], min2[12], flag;
	double snr, Eb_No, noise_pwr, noise[2], **Y, **Yk, *err_rate, Ps;
	double P, sum, sum0, sum1, pdf, ***R_pn, ***Q_np, *sym_I, *sym_Q, *I, *Q, **Rk;
	double ***R_tmp, **Out, AVE_CH;
	FILE *ber, *records;

	start = time(NULL);
	printf("BER Performance of Bit-based Message Passing in MIMO Channel\n");
	printf("Receiver Channel Selection\n");
	printf("Number of transmit antenna is %d\n", T);
	printf("Number of receive antenna is %d\n", R);
	printf("Maximum number of bits of simulation = %d\n\n", Mc*T*num_packet);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Record_Bit-based_Message_Passing.log", "a");
	fprintf(records, "BER Performance of Bit-based Message Passing in MIMO Channel\n");
	fprintf(records, "Receiver Channel Selection\n");
	fprintf(records, "Number of transmit antenna is %d\n", T);
	fprintf(records, "Number of receive antenna is %d\n", R);
	fprintf(records, "Maximum number of bits of simulation = %d\n\n", Mc*T*num_packet);
	fprintf(records, "Eb/No     BER\n");
	fflush(records);

	data_bit = new int[Mc*T];
	data = new int[Mc*T];
	sym_I = new double[T];
	sym_Q = new double[T];
	I = new double[T];
	Q = new double[T];
	Ak = new int[Mc*T];
	Yk = new double*[R];
	for(i=0; i<R; i++)
   		Yk[i] = new double[2];
	Y = new double*[R];
	for(i=0; i<R; i++)
   		Y[i] = new double[2];
	Rk = new double*[R];
	for(i=0; i<R; i++)
   		Rk[i] = new double[2];
	R_pn = new double**[R];
	for(i=0; i<R; i++)
   		R_pn[i] = new double*[Mc*T];
	for(i=0; i<R; i++)
   		for(j=0; j<Mc*T; j++)
   			R_pn[i][j] = new double[2];
	R_tmp = new double**[R];
	for(i=0; i<R; i++)
   		R_tmp[i] = new double*[Mc*T];
	for(i=0; i<R; i++)
   		for(j=0; j<Mc*T; j++)
   			R_tmp[i][j] = new double[2];
	Q_np = new double**[Mc*T];
	for(i=0; i<Mc*T; i++)
   		Q_np[i] = new double*[R];
	for(i=0; i<Mc*T; i++)
   		for(j=0; j<R; j++)
   			Q_np[i][j] = new double[2];
	Out = new double*[Mc*T];
	for(i=0; i<Mc*T; i++)
   		Out[i] = new double[2];
	err_count = new int[Iteration];
	err_rate = new double[Iteration];

	srand((unsigned) time(&t));

/*======================*/
/* main simulation loop            */
/*======================*/
	fopen_s(&ber, "BER_Bit_based_Message_Passing.log", "w");
	for(snr=0; snr<=16; snr+=2)
	{
   		for(s=0; s<Iteration; s++)
   			err_count[s] = 0;

   		// noise power calculation
		Eb_No = (double)snr;
		Ps = 1 * T;
		noise_pwr = 0.5*Ps*R/(T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption
		printf("Eb_No = %f\n", Eb_No);

		p = 0;
		do
		{
   			for(i=0; i<T; i++)
   			{
				// Generate random information bit stream
				if(rand()/(float)RAND_MAX>=0.5)
					data_bit[2*i] = 1;
				else
					data_bit[2*i] = 0;

				if(rand()/(float)RAND_MAX>=0.5)
					data_bit[2*i+1] = 1;
				else
					data_bit[2*i+1] = 0;

				sym_I[i] = (2*data_bit[2*i]-1)/sqrt(2.0);
         		sym_Q[i] = (2*data_bit[2*i+1]-1)/sqrt(2.0);
			}

			// Generate iid. Complex Gaussian MIMO Channel Coefficients
			for(i=0; i<R; i++)
         		for(j=0; j<T; j++)
				{
            		AWGN_noise(0, 1.0, &noise[0]);
            		H_I[i][j] = noise[0];
					H_Q[i][j] = noise[1];
				}

			// MIMO Fading Channel
			for(i=0; i<R; i++)
			{
         		Yk[i][0] = Yk[i][1] = 0.0;
         		for(j=0; j<T; j++)
				{
         			Yk[i][0] += (sym_I[j] * H_I[i][j] - sym_Q[j] * H_Q[i][j]);
					Yk[i][1] += (sym_I[j] * H_Q[i][j] + sym_Q[j] * H_I[i][j]);
				}
			}

			/* AWGN channel */
			for(i=0; i<R; i++)
			{
				AWGN_noise(0, noise_pwr, &noise[0]);
         		Rk[i][0] = Yk[i][0] + noise[0];
				Rk[i][1] = Yk[i][1] + noise[1];
			}

/*==================================================*/
/*  B i t - b a s e d    M e s s a g e    P a s s i n g    A l g o r i t h m  */
/*==================================================*/

			// Initilization of extrinsic probability
			for(i=0; i<Mc*T; i++)
   				for(j=0; j<R; j++)
      				for(k=0; k<2; k++)
   						Q_np[i][j][k] = 0.5;

			// Receiver Channel Selection (Criterion: Below Average)
			AVE_CH = 0.0;
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
				{
					RCS[i][j] = 1;
					AVE_CH += (pow(H_I[i][j],2.0) + pow(H_Q[i][j],2.0));
				}

			AVE_CH /= float(R*T);
			
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
					if((pow(H_I[i][j],2.0)+pow(H_Q[i][j],2.0)) < AVE_CH/5.0)
						RCS[i][j] = 0;
/*
			// Receiver Channel Selection (Criterion: delete the weakest channels)
			for(i=0; i<R; i++)
				for(j=0; j<T; j++)
					RCS[i][j] = 1;

			for(d=0; d<4; d++)
			{
				min1[d]= 0;
				min2[d] = 0;
				for(z=0; z<d; z++)
					if( (min1[d] == min1[z]) && (min2[d] == min2[z]) )
					{
						if(min2[d]++ < T)
							min2[d]++;
						else
						{
							min1[d]++;
							min2[d] = 0;
						}
					}

				for(i=0; i<R; i++)
					for(j=0; j<T; j++)
					{
						flag = 1;
						for(z=0; z<d; z++)
							if( (i == min1[z]) && (j == min2[z]) )
							{
								flag = 0;
								break;
							}

						if(flag == 1)
						{
							if((pow(H_I[i][j],2.0)+pow(H_Q[i][j],2.0)) <= (pow(H_I[min1[d]][min2[d]],2.0)+pow(H_Q[min1[d]][min2[d]],2.0)))
							{
								min1[d] = i;
								min2[d] = j;
							}
						}
					}

				RCS[min1[d]][min2[d]] = 0;
			}
*/
			for(s=0; s<Iteration; s++)
			{
         		for(k=0; k<2; k++)		// Probability of 1 and 0
            		for(i=0; i<R; i++)
            			for(j=0; j<T; j++)
							if(RCS[i][j] == 1)
								for(m=0; m<Mc; m++)
               					{
               						R_tmp[i][j*Mc+m][k] = 0.0;
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
			 	            				I[r] = (2*data[2*r]-1)/sqrt(2.0);
											Q[r] = (2*data[2*r+1]-1)/sqrt(2.0);
               							}
	
		                 				Y[i][0] = Y[i][1] = 0.0;
										for(r=0; r<T; r++)
										{
					        				Y[i][0] += (I[r] * H_I[i][r] - Q[r] * H_Q[i][r]);
											Y[i][1] += (I[r] * H_Q[i][r] + Q[r] * H_I[i][r]);
										}
									
                     					P = 1.0;
                     					for(r=0; r<T; r++)
											for(n=0; n<Mc; n++)
                     							if((RCS[i][r] == 1) && (r*Mc+n != j*Mc+m))
                     								P *= Q_np[r*Mc+n][i][data[r*Mc+n]];

										sum = 0.0;
										for(z=0; z<2; z++)
			                				sum += pow(Rk[i][z]-Y[i][z],2.0);

										pdf = exp(-sum/noise_pwr)/(Pi*noise_pwr);

					     				R_tmp[i][j*Mc+m][k] += pdf * P;
                  					}
               					}	
						
				for(i=0; i<R; i++)
            		for(j=0; j<T; j++)
						if(RCS[i][j] == 1)
							for(m=0; m<Mc; m++)
							{
								R_pn[i][j*Mc+m][1] = R_tmp[i][j*Mc+m][1]/(R_tmp[i][j*Mc+m][0] + R_tmp[i][j*Mc+m][1]);
          						R_pn[i][j*Mc+m][0] = R_tmp[i][j*Mc+m][0]/(R_tmp[i][j*Mc+m][0] + R_tmp[i][j*Mc+m][1]);
							}

				for(j=0; j<T; j++)
					for(i=0; i<R; i++)
						if(RCS[i][j] == 1)
							for(m=0; m<Mc; m++)           				
							{
		           				sum0 = sum1 = 1.0;
			       				for(n=0; n<R; n++)
           							if((n != i) && (RCS[n][j] == 1))
									{
						 				sum1 *= R_pn[n][j*Mc+m][1];
										sum0 *= R_pn[n][j*Mc+m][0];
									}

								Q_np[j*Mc+m][i][1] = sum1/(sum1 + sum0);
								Q_np[j*Mc+m][i][0] = sum0/(sum1 + sum0);
			 				}

				// Decision Output
				for(j=0; j<T; j++)
					for(m=0; m<Mc; m++)
					{
       					sum0 = sum1 = 1.0;
						for(i=0; i<R; i++)
							if(RCS[i][j] == 1)
							{	
								sum1 *= R_pn[i][j*Mc+m][1];
								sum0 *= R_pn[i][j*Mc+m][0];
							}

						Out[j*Mc+m][1] = sum1/(sum1 + sum0);
						Out[j*Mc+m][0] = sum0/(sum1 + sum0);
					}

				for(j=0; j<Mc*T; j++)
				{
         			if(Out[j][1] >= Out[j][0])
            			Ak[j] = 1;
       				else
            			Ak[j] = 0;

            		err_count[s] += Error_count(data_bit[j], Ak[j]);
				}
			}

			p++;
		} while(err_count[Iteration-1] <= 1000 && p < num_packet);

		// Statistics and records
		cout << "Error Rate = ";
		for(s=0; s<Iteration; s++)
		{
      		err_rate[s] = err_count[s] / (double)(Mc*T*p);
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

	delete data_bit;
	delete data;
	delete sym_I;
	delete sym_Q;
	delete I;
	delete Q;
	delete Ak;
	for(i=0; i<R; i++)
   		delete Yk[i];
	delete Yk;
	for(i=0; i<R; i++)
   		delete Y[i];
	delete Y;
	for(i=0; i<R; i++)
   		delete Rk[i];
	delete Rk;
	for(i=0; i<R; i++)
   		for(j=0; j<Mc*T; j++)
   			delete R_pn[i][j];
	for(i=0; i<R; i++)
   		delete R_pn[i];
	delete R_pn;
	for(i=0; i<R; i++)
   		for(j=0; j<Mc*T; j++)
   			delete R_tmp[i][j];
	for(i=0; i<R; i++)
   		delete R_tmp[i];
	delete R_tmp;
	for(i=0; i<Mc*T; i++)
   		for(j=0; j<R; j++)
   			delete Q_np[i][j];
	for(i=0; i<Mc*T; i++)
   		delete Q_np[i];
	delete Q_np;
	for(i=0; i<Mc*T; i++)
   		delete Out[i];
	delete Out;
	delete err_count;
	delete err_rate;

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
