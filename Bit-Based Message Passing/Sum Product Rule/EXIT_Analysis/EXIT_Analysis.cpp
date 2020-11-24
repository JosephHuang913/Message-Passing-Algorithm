/*========================================================*/
/* Author: Chao-wang Huang                                                                       */
/* Date: Monday, September 22, 2008                                                        */
/* A Bit-based Message Passing Algorithm for MIMO Channel is simulated */
/* EXIT Chart Analysis                                                                                */
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

void AWGN_noise(double, double, double *);
int Error_count(int, int);
void Histgram(double, double, double, int, int *);

#define Pi 	3.14159265358979
#define num_packet 100000			// number of packets simulated
const int Mc = 2;								// modulation order (2 for QPSK)
const int T = 4;									// Number of transmit antenna
const int R = 4;									// Number of receive antenna
double H_I[4][4], H_Q[4][4];
int RCS[4][4];

int main(void)
{
	time_t  t, start, end;
	int i, j, k, l, m, n, p, r, x, z, *data, *data_bit, min1, min2, num1, num0;
	int *Hist_E1, *Hist_E2, *Hist_E3, *Hist_E4, *Hist_E5, *Hist_E6;
	double snr, Eb_No, noise_pwr, noise[2], **Y, **Yk, Ps;
	double P, sum, sum0, sum1, pdf, ***R_pn, ***Q_np, *sym_I, *sym_Q, *I, *Q, **Rk;
	double ***R_tmp, AVE_CH, MI, A1, A2, B1, B2, C1, C2;
	double Sigma_a, Var_a, Mean_a, LLR_a[2], ***P_a, LLR_e;
	double *pdf_E1, *pdf_E2, *pdf_E3, *pdf_E4, *pdf_E5, *pdf_E6, MI_e1, MI_e2, MI_d;
	double **Out;
	FILE *exit1, *exit2, *pdf_E, *records;

	start = time(NULL);
	printf("EXIT Chart of Bit-based Message Passing in MIMO Channel\n");
	//printf("Receiver Channel Selection\n");
	printf("Number of transmit antenna is %d\n", T);
	printf("Number of receive antenna is %d\n", R);
	printf("Maximum number of bits of simulation = %d\n\n", Mc*T*num_packet);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Record_Bit-based_Message_Passing.log", "a");
	fprintf(records, "EXIT Chart of Bit-based Message Passing in MIMO Channel\n");
	//fprintf(records, "Receiver Channel Selection\n");
	fprintf(records, "Number of transmit antenna is %d\n", T);
	fprintf(records, "Number of receive antenna is %d\n", R);
	fprintf(records, "Maximum number of bits of simulation = %d\n\n", Mc*T*num_packet);
	fprintf(records, "Eb/No     I_a      I_e,fnd  I,e_bnd\n");
	fflush(records);

	data_bit = new int[Mc*T];
	data = new int[Mc*T];
	sym_I = new double[T];
	sym_Q = new double[T];
	I = new double[T];
	Q = new double[T];
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
	P_a = new double**[Mc*T];
	for(i=0; i<Mc*T; i++)
		P_a[i] = new double*[R];
	for(i=0; i<Mc*T; i++)
		for(j=0; j<R; j++)
			P_a[i][j] = new double[2];
	Out = new double*[Mc*T];
	for(i=0; i<Mc*T; i++)
   		Out[i] = new double[2];
	Hist_E1 = new int[10000];
	Hist_E2 = new int[10000];
	Hist_E3 = new int[10000];
	Hist_E4 = new int[10000];
	Hist_E5 = new int[10000];
	Hist_E6 = new int[10000];
	pdf_E1 = new double[10000];
	pdf_E2 = new double[10000];
	pdf_E3 = new double[10000];
	pdf_E4 = new double[10000];
	pdf_E5 = new double[10000];
	pdf_E6 = new double[10000];

	srand((unsigned) time(&t));

/*======================*/
/* main simulation loop            */
/*======================*/
	fopen_s(&exit1, "Func_exit.log", "w");
	fopen_s(&exit2, "Bit_exit.log", "w");
	fopen_s(&pdf_E, "LLR_pdf.log", "w");
	for(snr=10; snr<=10; snr+=3)
	{
   		// noise power calculation
		Eb_No = (double)snr;
		Ps = 1 * T;
		noise_pwr = 0.5*Ps*R/(T*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption
		
		for(MI=0.0; MI<1.0; MI+=0.099)
		{
			// Inverse J-function
			A1 = 1.09542; B1 = 0.214217; C1 = 2.33727;
			A2 = 0.706692; B2 = 0.386013; C2 = -1.75017;
			if(MI <= 0.3646)
				Sigma_a = A1*pow(MI,2.0) + B1*MI + C1*sqrt(MI);
			else
				Sigma_a = -A2 * log(-B2 * (MI - 1.0)) - C2 * MI;

			Var_a = pow(Sigma_a, 2.0);
			Mean_a = Var_a / 2.0;
			printf("Ia = %f, Eb_No = %f\n", MI, Eb_No);
			p = num1 = num0 = 0;
			for(i=0; i<10000; i++)
				Hist_E1[i] = Hist_E2[i] = Hist_E3[i] = Hist_E4[i] = Hist_E5[i] = Hist_E6[i] = 0;

			do
			{
   				for(i=0; i<T; i++)
   				{
					// Generate random information bit stream (QPSK)
					if(rand()/(float)RAND_MAX>=0.5)
					{
						data_bit[2*i] = 1;
						num1++;
					}
					else
					{
						data_bit[2*i] = 0;
						num0++;
					}

					if(rand()/(float)RAND_MAX>=0.5)
					{
						data_bit[2*i+1] = 1;
						num1++;
					}
					else
					{
						data_bit[2*i+1] = 0;
						num0++;
					}

					sym_I[i] = (2*data_bit[2*i]-1)/sqrt(2.0);
         			sym_Q[i] = (2*data_bit[2*i+1]-1)/sqrt(2.0);
				}

				// Initilization of intrinsic probability
				for(i=0; i<Mc*T; i++)
   					for(j=0; j<R; j++)
					{
						AWGN_noise(Mean_a*(2.0*data_bit[i]-1.0), 2.0*Var_a, &LLR_a[0]);
      					for(k=0; k<2; k++)		
   							P_a[i][j][k] = exp(k*LLR_a[0]) / (1 + exp(LLR_a[0]));
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
/*  Bit-based Message Passing Algorithm & EXIT Chart Analysis    */
/*==================================================*/
/*
				// Receiver Channel Selection (Criterion: Above Average)
				for(i=0; i<R; i++)
					for(j=0; j<T; j++)
						RCS[i][j] = 1;

				for(j=0; j<T; j++)
				{
					AVE_CH = 0.0;
					for(i=0; i<R; i++)
						AVE_CH += (1.0/(float)R) * (pow(H_I[i][j],2.0) + pow(H_Q[i][j],2.0));
					
					for(i=0; i<R; i++)
						if((pow(H_I[i][j],2.0)+pow(H_Q[i][j],2.0)) < AVE_CH)
							RCS[i][j] = 0;
				}
*/
				// Receiver Channel Selection (Criterion: Strongest)
				for(i=0; i<R; i++)
					for(j=0; j<T; j++)
						RCS[i][j] = 1;
/*
				for(j=0; j<T; j++)
				{
					min1= 0;
					for(i=1; i<R; i++)
						if((pow(H_I[i][j],2.0)+pow(H_Q[i][j],2.0)) <= (pow(H_I[min1][j],2.0)+pow(H_Q[min1][j],2.0)))
							min1 = i;

					if(min1 == 0)
						min2 = 1;
					else
						min2 = 0;
					for(i=min2+1; i<R; i++)
					{
						if(i == min1)
							;
						else if((pow(H_I[i][j],2.0)+pow(H_Q[i][j],2.0)) <= (pow(H_I[min2][j],2.0)+pow(H_Q[min2][j],2.0)))
							min2 = i;
					}
	
					//RCS[min1][j] = RCS[min2][j] = 0;
					RCS[min1][j] = 0;
				}
*/
				// Extrinsic probability calculation of function nodes
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
                     								P *= P_a[r*Mc+n][i][data[r*Mc+n]];

										sum = 0.0;
										for(z=0; z<2; z++)
			                				sum += pow(Rk[i][z]-Y[i][z],2.0);

										pdf = exp(-sum/noise_pwr)/(Pi*noise_pwr);

					     				R_tmp[i][j*Mc+m][k] += pdf * P;
                  					}
               					}	
						
				// Normalization of extrinsic probability
				for(i=0; i<R; i++)
            		for(j=0; j<T; j++)
						if(RCS[i][j] == 1)
							for(m=0; m<Mc; m++)
							{
								R_pn[i][j*Mc+m][1] = R_tmp[i][j*Mc+m][1]/(R_tmp[i][j*Mc+m][0] + R_tmp[i][j*Mc+m][1]);
          						R_pn[i][j*Mc+m][0] = R_tmp[i][j*Mc+m][0]/(R_tmp[i][j*Mc+m][0] + R_tmp[i][j*Mc+m][1]);
								LLR_e = log(R_pn[i][j*Mc+m][1] / R_pn[i][j*Mc+m][0]);
								
								// Histgram of extrinsic information
								if(data_bit[j*Mc+m] == 1)
									Histgram(LLR_e, -500.0, 500.0, 10000, Hist_E1);
								else
									Histgram(LLR_e, -500.0, 500.0, 10000, Hist_E2);
							}

				// Extrinsic probability calculation of bit nodes
				for(j=0; j<T; j++)
					for(i=0; i<R; i++)
						if(RCS[i][j] == 1)
							for(m=0; m<Mc; m++)           				
							{
		           				sum0 = sum1 = 1.0;
			       				for(n=0; n<R; n++)
           							if((n != i) && (RCS[n][j] == 1))
									{
						 				//sum1 *= P_a[j*Mc+m][n][1];
										//sum0 *= P_a[j*Mc+m][n][0];
										sum1 *= R_pn[n][j*Mc+m][1];
										sum0 *= R_pn[n][j*Mc+m][0];
									}

								// Normalization of extrinsic probability
								Q_np[j*Mc+m][i][1] = sum1/(sum1 + sum0);
								Q_np[j*Mc+m][i][0] = sum0/(sum1 + sum0);
								LLR_e = log(Q_np[j*Mc+m][i][1] / Q_np[j*Mc+m][i][0]);
								
								// Histgram of extrinsic information
								if(data_bit[j*Mc+m] == 1)
									Histgram(LLR_e, -500.0, 500.0, 10000, Hist_E3);
								else
									Histgram(LLR_e, -500.0, 500.0, 10000, Hist_E4);
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
								//sum1 *= P_a[j*Mc+m][i][1];
								//sum0 *= P_a[j*Mc+m][i][0];
							}

						Out[j*Mc+m][1] = sum1/(sum1 + sum0);
						Out[j*Mc+m][0] = sum0/(sum1 + sum0);
						LLR_e = log(Out[j*Mc+m][1] / Out[j*Mc+m][0]);
								
						// Histgram of decision information
						if(data_bit[j*Mc+m] == 1)
							Histgram(LLR_e, -500.0, 500.0, 10000, Hist_E5);
						else
							Histgram(LLR_e, -500.0, 500.0, 10000, Hist_E6);
					}

				p++;
			} while(p < num_packet);

			// Statistics (pdf) of extrinsic information
			for(i=0; i<10000; i++)
			{
				pdf_E1[i] = (Hist_E1[i] / (float)(num1*R)) / 0.1;
				pdf_E2[i] = (Hist_E2[i] / (float)(num0*R)) / 0.1;
				pdf_E3[i] = (Hist_E3[i] / (float)(num1*R)) / 0.1;
				pdf_E4[i] = (Hist_E4[i] / (float)(num0*R)) / 0.1;
				pdf_E5[i] = (Hist_E5[i] / (float)(num1)) / 0.1;
				pdf_E6[i] = (Hist_E6[i] / (float)(num0)) / 0.1;
				
				if(MI == 0.495)
					fprintf(pdf_E, "%f %f %f %f %f %f %f %f %f\n", MI, Eb_No, (i-5000)/10.0, pdf_E1[i], pdf_E2[i], pdf_E3[i], pdf_E4[i], pdf_E5[i], pdf_E6[i]);
			}
			fflush(pdf_E);

			// Mutual information of the Extrinsic output
			MI_e1 = MI_e2 = MI_d = 0.0;
			for(i=0; i<10000; i++)
			{
				if(pdf_E1[i] != 0.0)
					MI_e1 += 0.5 * (pdf_E1[i]*(log(2.0*pdf_E1[i]/(pdf_E1[i]+pdf_E2[i]))/log(2.0))*0.1);
				if(pdf_E2[i] != 0.0)
					MI_e1 += 0.5 * (pdf_E2[i]*(log(2.0*pdf_E2[i]/(pdf_E1[i]+pdf_E2[i]))/log(2.0))*0.1);
				if(pdf_E3[i] != 0.0)
					MI_e2 += 0.5 * (pdf_E3[i]*(log(2.0*pdf_E3[i]/(pdf_E3[i]+pdf_E4[i]))/log(2.0))*0.1);
				if(pdf_E4[i] != 0.0)
					MI_e2 += 0.5 * (pdf_E4[i]*(log(2.0*pdf_E4[i]/(pdf_E3[i]+pdf_E4[i]))/log(2.0))*0.1);
				if(pdf_E5[i] != 0.0)
					MI_d += 0.5 * (pdf_E5[i]*(log(2.0*pdf_E5[i]/(pdf_E5[i]+pdf_E6[i]))/log(2.0))*0.1);
				if(pdf_E6[i] != 0.0)
					MI_d += 0.5 * (pdf_E6[i]*(log(2.0*pdf_E6[i]/(pdf_E5[i]+pdf_E6[i]))/log(2.0))*0.1);
			}

			fprintf(exit1, "%f %f %f\n", Eb_No, MI, MI_e1);
			fprintf(exit2, "%f %f %f\n", MI_e1, MI_e2, MI_d);
			fprintf(records, "%f %f %f %f %f\n", Eb_No, MI, MI_e1, MI_e2, MI_d);
			fflush(exit1);
			fflush(exit2);
			fflush(records);
		}
	}

	delete data_bit;
	delete data;
	delete sym_I;
	delete sym_Q;
	delete I;
	delete Q;
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
		for(j=0; j<R; j++)
			delete P_a[i][j];
	for(i=0; i<Mc*T; i++)
		delete P_a[i];
	delete P_a;
	for(i=0; i<Mc*T; i++)
   		delete Out[i];
	delete Out;
	delete Hist_E1;
	delete Hist_E2;
	delete Hist_E3;
	delete Hist_E4;
	delete Hist_E5;
	delete Hist_E6;
	delete pdf_E1;
	delete pdf_E2;
	delete pdf_E3;
	delete pdf_E4;
	delete pdf_E5;
	delete pdf_E6;
	
	end = time(NULL);
	printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
	fprintf(records, "Total elapsed time: %.0f(sec)\n\n", difftime(end,start));
	fclose(exit1);
	fclose(exit2);
	fclose(pdf_E);
	fclose(records);
	printf("This program is ended. Press any key to continue.\n");
	getchar();

	return 0;
}

void AWGN_noise(double mu, double variance, double *noise)
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

void Histgram(double in, double left, double right, int section, int *statistic)
{
	int i;

   if(in<left || in>=right)
   {
   		printf("Warning!! The input number is out of range.\n");
		cout << in << endl;
		//getchar();
   }
   else
   {
   	for(i=0; i<section; i++)
   		if(in>=left+((right-left)/(float)section)*i && in<left+((right-left)/(float)section)*(i+1))
         {
      		*(statistic+i) = *(statistic+i) + 1;
            break;
         }
   }
}
