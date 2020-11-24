/*=====================================================================*/
/* Author: Chao-wang Huang                                             */
/* Date: Saturday, December 08, 2006                                   */
/* A Bit-based Message Passing Algorithm for MIMO Channel is simulated */
/*=====================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>

void AWGN_noise(float, double, double *);
int Error_count(int, int);

#define Pi 	3.14159265358979
#define num_packet 1000000							// number of packets simulated
const int m = 2;										// modulation order (2 for QPSK)
const int T = 4;										// Number of transmit antenna
const int R = 4;                             // Number of receive antenna
const int Iteration = 10;

const float H_I[4][4] = {{0.0, 1.0, 0.0, 1.0},	// Channel matrix
								 {1.0, 0.0, 1.0, 0.0},
                         {0.0, 0.0, 1.0, 1.0},
                         {1.0, 1.0, 0.0, 0.0}
								};
const float H_Q[4][4] = {{0.0, 1.0, 0.0, 1.0},
								 {1.0, 0.0, 1.0, 0.0},
                         {0.0, 0.0, 1.0, 1.0},
                         {1.0, 1.0, 0.0, 0.0}
								};

//double H_I[4][4], H_Q[4][4];

int main(void)
{
	time_t  t, start, end;
	int i, j, k, l, n, p, r, s, x, z, *data, *data_bit, *Ak, *err_count;
   double snr, Eb_No, noise_pwr, noise[2], **Y, **Yk, *err_rate, Ps;
   double P, sum, sum0, sum1, pdf, ***R_pn, ***Q_np, *sym_I, *sym_Q, *I, *Q, **Rk;
   double ***R_tmp, **Out;
   FILE *ber, *records;

   start = time(NULL);
   printf("BER Performance of Bit-based Message Passing in MIMO Channel\n");
	printf("Number of transmit antenna is %d\n", T);
   printf("Number of receive antenna is %d\n", R);
   printf("Maximum number of bits of simulation = %d\n\n", m*T*num_packet);
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("record_Bit-based_Message_Passing.log", "a");
   fprintf(records, "BER Performance of Bit-based Message Passing in MIMO Channel\n");
	fprintf(records, "Number of transmit antenna is %d\n", T);
   fprintf(records, "Number of receive antenna is %d\n", R);
   fprintf(records, "Maximum number of bits of simulation = %d\n\n", m*T*num_packet);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[m*T];
   data = new int[m*T];
   sym_I = new double[T];
   sym_Q = new double[T];
   I = new double[T];
   Q = new double[T];
   Ak = new int[m*T];
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
   	R_pn[i] = new double*[m*T];
   for(i=0; i<R; i++)
   	for(j=0; j<m*T; j++)
   		R_pn[i][j] = new double[2];
   R_tmp = new double**[R];
   for(i=0; i<R; i++)
   	R_tmp[i] = new double*[m*T];
   for(i=0; i<R; i++)
   	for(j=0; j<m*T; j++)
   		R_tmp[i][j] = new double[2];
   Q_np = new double**[m*T];
   for(i=0; i<m*T; i++)
   	Q_np[i] = new double*[R];
   for(i=0; i<m*T; i++)
   	for(j=0; j<R; j++)
   		Q_np[i][j] = new double[2];
   Out = new double*[m*T];
   for(i=0; i<m*T; i++)
   	Out[i] = new double[2];
   err_count = new int[Iteration];
   err_rate = new double[Iteration];

	srand((unsigned) time(&t));

/*======================*/
/* main simulation loop */
/*======================*/
	ber=fopen("ber_Bit-based_Message_Passing.log", "w");
   for(snr=0; snr<=16; snr++)
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
      //for(p=0; p<num_packet; p++)
      {
   		for(i=0; i<T; i++)
   		{
         	data_bit[2*i] = random(2);		// Generate random information bit stream
            data_bit[2*i+1] = random(2);

            sym_I[i] = (2*data_bit[2*i]-1)/sqrt(2.0);
         	sym_Q[i] = (2*data_bit[2*i+1]-1)/sqrt(2.0);
         }
         /*
         // Generate iid. Complex Gaussian MIMO Channel Coefficients
         for(i=0; i<R; i++)
         	for(j=0; j<T; j++)
            {
            	AWGN_noise(0, 1.0, &noise[0]);
            	H_I[i][j] = noise[0];
               H_Q[i][j] = noise[1];
            }
         
         // Record for y_p
         for(x=0; x<pow(2,m*T); x++)
         {
            for(n=0; n<m*T; n++)
            	data[n] = (x>>n) & 1;

            for(i=0; i<T; i++)
            {
            	I[i] = (2*data[2*i]-1)/sqrt(2.0);
         		Q[i] = (2*data[2*i+1]-1)/sqrt(2.0);
            }

            for(i=0; i<R; i++)
            {
         		Yp[i][x][0] = Yp[i][x][1] = 0.0;
	         	for(j=0; j<T; j++)
   	         {
      	   		Yp[i][x][0] += (I[j] * H_I[i][j] - Q[j] * H_Q[i][j]);
         	      Yp[i][x][1] += (I[j] * H_Q[i][j] + Q[j] * H_I[i][j]);
            	}
         	}
         }
           */
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

/*============================================================================*/
/*  B i t - b a s e d    M e s s a g e    P a s s i n g    A l g o r i t h m  */
/*============================================================================*/

			// Initilization of a-priori probability
			for(i=0; i<m*T; i++)
   			for(j=0; j<R; j++)
      			for(k=0; k<2; k++)
   					Q_np[i][j][k] = 0.5;

			for(s=0; s<Iteration; s++)
         {
         	for(k=0; k<2; k++)		// Probability of 1 and 0
            	for(i=0; i<R; i++)
            		for(j=0; j<m*T; j++)
               	{
               		R_tmp[i][j][k] = 0.0;
	               	for(x=0; x<pow(2,m*T-1); x++)
   	               {
      	            	l = 0;
         	            for(n=0; n<m*T; n++)
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
                        	Y[i][0] += (I[r] * H_I[i][r] - Q[r] * H_Q[i][r]);
                           Y[i][1] += (I[r] * H_Q[i][r] + Q[r] * H_I[i][r]);
                        }

                     	P = 1.0;
                     	for(n=0; n<m*T; n++)
                     		if(n!=j)
                     			P *= Q_np[n][i][data[n]];
                        /*
                     	for(q=0; q<pow(2,m*T); q++)
   	                     if(Yp[i][q][0] == Y[i][0] && Yp[i][q][1] == Y[i][1])
                           {
                           	P_xn = 1.0;
                              P_yp = 0.5;

                              sum = 0.0;
					               for(z=0; z<2; z++)
               						sum += pow(Rk[i][z]-Y[i][z],2.0);

							         pdf = exp(-sum/noise_pwr)/(Pi*noise_pwr);

                              R[i][j][k] += P_xn * pdf * P_yp * P;
                           }
                        */
                        sum = 0.0;
                        for(z=0; z<2; z++)
                        	sum += pow(Rk[i][z]-Y[i][z],2.0);

                        pdf = exp(-sum/noise_pwr)/(Pi*noise_pwr);

                     	R_tmp[i][j][k] += pdf * P;
                  	}
               	}

            for(i=0; i<R; i++)
            	for(j=0; j<m*T; j++)
               {
         			R_pn[i][j][1] = R_tmp[i][j][1]/(R_tmp[i][j][0] + R_tmp[i][j][1]);
            		R_pn[i][j][0] = R_tmp[i][j][0]/(R_tmp[i][j][0] + R_tmp[i][j][1]);
               }

            for(j=0; j<m*T; j++)
            	for(i=0; i<R; i++)
               {
               	sum0 = sum1 = 1.0;
               	for(n=0; n<R; n++)
                  	if(n != i)
                     {
                     	sum1 *= R_pn[n][j][1];
                        sum0 *= R_pn[n][j][0];
                     }

                  Q_np[j][i][1] = sum1/(sum1 + sum0);
                  Q_np[j][i][0] = sum0/(sum1 + sum0);
         		}

            // Decision Output
            for(j=0; j<m*T; j++)
            {
            	sum0 = sum1 = 1.0;
               for(i=0; i<R; i++)
               {
	               sum1 *= R_pn[i][j][1];
                  sum0 *= R_pn[i][j][0];
               }

               Out[j][1] = sum1/(sum1 + sum0);
               Out[j][0] = sum0/(sum1 + sum0);
            }

            for(j=0; j<m*T; j++)
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
      	err_rate[s] = err_count[s] / (double)(m*T*p);
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
   	for(j=0; j<m*T; j++)
   		delete R_pn[i][j];
   for(i=0; i<R; i++)
   	delete R_pn[i];
   delete R_pn;
   for(i=0; i<R; i++)
   	for(j=0; j<m*T; j++)
   		delete R_tmp[i][j];
   for(i=0; i<R; i++)
   	delete R_tmp[i];
   delete R_tmp;
   for(i=0; i<m*T; i++)
   	for(j=0; j<R; j++)
   		delete Q_np[i][j];
   for(i=0; i<m*T; i++)
   	delete Q_np[i];
   delete Q_np;
   for(i=0; i<m*T; i++)
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
   getch();

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

