#include "stdafx.h"
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <iostream>
using namespace std;

double Phi(double );

void SPA_decoder(int max_itr, int column, int row, double *receive_y, int **h, int *decoded_bit, double snr_rms)
{

	int i,j,k,count,FLAG_STOP=0,itr=1;
	int temp;
	double sum,sum_Q,temp_afa,temp_beta,sign;
	double *Q,**q,**r,*Lpi;
	int **mark_v,**mark_h;
	int *wc_per_column,*wr_per_row;
	
	
	Q= new double[column];
	q = new double*[column];
	for(i=0;i<column;i++)
		q[i] = new double[row];
	r = new double*[row];
	for(i=0;i<row;i++)
		r[i] = new double[column];
	Lpi = new double[column];

	/***************  find_ones  ***********************/

	mark_v = new int*[column];
	for(i=0;i<column;i++)
		mark_v[i] = new int[6];   // maximum wc is 6
	mark_h = new int*[row];
	for(i=0;i<row;i++)
		mark_h[i] = new int[24];   // maximum wr is 20
	wc_per_column = new int[column];
	wr_per_row = new int[row];

	for(i=0;i<column;i++)
		for(j=0;j<6;j++)
			mark_v[i][j]=-1;
	
	for(i=0;i<row;i++)
		for(j=0;j<24;j++)
			mark_h[i][j]=-1;

	for(i=0;i<column;i++)
	{
		count=0;
		for(j=0;j<row;j++) 
		{
			if(h[j][i]==1)
			{
				mark_v[i][count]=j;
				count++;
			}
		}
		wc_per_column[i]=count;
	}

	for(j=0;j<row;j++)
	{
		count=0;
		for(i=0;i<column;i++) 
		{
			if(h[j][i]==1)
			{	mark_h[j][count]=i;
				count++;
			}
		}
		wr_per_row[j]=count;
	}

	/***************  Start Decoding  ***********************/


	
	for (i=0;i<column;i++)
	{
		for (j=0;j<row;j++)
		{
			q[i][j]=0;
			r[j][i]=0;
		}
		//Lpi[i]=2*receive_y[i]/(pow(sigma_N,2));
		Lpi[i] = receive_y[i]*snr_rms;
	}

	for(i=0;i<column;i++)                                        //step 1 : initial q
	{
		for(k=0;k<wc_per_column[i];k++)
		{	
			temp=mark_v[i][k];
			q[i][temp]=Lpi[i];
		}
	}
	
	while(itr<=max_itr /*&& FLAG_STOP==0*/)
	{			
		for(j=0;j<row;j++)                                      //step 2: update r                                   
		{
			for(k=0;k<wr_per_row[j];k++)
			{
				temp_afa=1.0;
				temp_beta=0.0;
				for(count=0;count<wr_per_row[j];count++)
				{
					temp=mark_h[j][count];
		
					if(count!=k)
					{
						if(q[temp][j]>=0)	
							sign=1.0;
						else
							sign=-1.0;		
						temp_afa = temp_afa*sign;
						temp_beta += Phi(sign*q[temp][j]);
					}
				}

				temp=mark_h[j][k];
				r[j][temp]=temp_afa*Phi(temp_beta);
			}
		}														                       
	
		for(i=0;i<column;i++)                                          //step 3 : update q & Q
		{
			sum_Q=0.0;
			for(k=0;k<wc_per_column[i];k++)
			{
				sum=0.0;
				for(count=0;count<wc_per_column[i];count++)
				{
					temp=mark_v[i][count];
					if(count!=k)
						sum+=r[temp][i];
				}
				temp=mark_v[i][k];
				sum_Q+=r[temp][i];
                q[i][temp] =  sum + Lpi[i];
			}
			Q[i] = sum_Q + Lpi[i];
			
			if(Q[i]<0)            //decision
				decoded_bit[i]=1;
			else
				decoded_bit[i]=0;

			err_count[iter] += Error_count(coded_bit[i], decoded_bit[i]);
		}

		for(j=0;j<row;j++)
		{
			temp=0;
			for(i=0;i<column;i++)
				temp+=(decoded_bit[i]*h[j][i]);              // c*h^t
			if((temp%2)!=0)
				break;
			if(j==row-1)
				FLAG_STOP=1;
		}
		itr++;
	}
	
	delete Q;
	for(i=0;i<column;i++)
		delete q[i];
	delete q;
	for(i=0;i<row;i++)
		delete r[i];
	delete r;
	delete Lpi;
	for(i=0;i<column;i++)
		delete mark_v[i];
	delete mark_v;
	for(i=0;i<row;i++)
		delete mark_h[i];
	delete mark_h;
	delete wc_per_column;
	delete wr_per_row;
}

double Phi(double beta)
{
	if(beta==0.0)
		return log((exp(1E-12)+1)/(exp(1E-12)-1));
	else
		return log((exp(beta)+1)/(exp(beta)-1));
}
		



