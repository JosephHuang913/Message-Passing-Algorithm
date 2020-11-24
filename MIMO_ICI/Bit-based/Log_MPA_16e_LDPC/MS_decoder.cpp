#include "stdafx.h"
# include <stdio.h>
# include <stdlib.h>
# include <iostream>
using namespace std;


void MS_decoder(int max_itr, int column, int row,double *receive_y,int **h,int *decoded_bit)
{

	int i,j,k,count,FLAG_STOP=0,itr=1;
	int temp;
	double sum,sum_Q,temp_afa,temp_beta,sign,temp_beta_1;
	double *Q,**q,**r,*Lpi;
	int **mark_v,**mark_h;
	int *wc_per_column,*wr_per_row;
	
	
	Q= (double*) calloc(column,sizeof(double));

	q = (double**) calloc(column,sizeof(double *));
	for(i=0;i<column;i++)
		q[i]= (double*) calloc(row,sizeof(double));

	r = (double**) calloc(row,sizeof(double *));
	for(i=0;i<row;i++)
		r[i]= (double*) calloc(column,sizeof(double));
	
	Lpi = (double *) calloc(column,sizeof(double));


	/***************  find_ones  ***********************/

	mark_v = (int **) calloc(column,sizeof(int *));
	for(i=0;i<column;i++)
		mark_v[i]= (int *) calloc(6,sizeof(int));   // maximum wc is 6

	mark_h = (int **) calloc(row,sizeof(int *));
	for(i=0;i<row;i++)
		mark_h[i]= (int *) calloc(24,sizeof(int));   // maximum wr is 20
	
	wc_per_column = (int *) calloc(column,sizeof(int));
	wr_per_row = (int *) calloc(row,sizeof(int));


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
		Lpi[i]=receive_y[i];

	}

	for(i=0;i<column;i++)                                        //step 1 : initial q
	{
		for(k=0;k<wc_per_column[i];k++)
		{	
			temp=mark_v[i][k];
			q[i][temp]=Lpi[i];


		}
	}
	
	while(itr<=max_itr && FLAG_STOP==0)
	{			
		for(j=0;j<row;j++)                                      //step 2: update r                                   
		{
			for(k=0;k<wr_per_row[j];k++)
			{
				temp_afa=1.0;
				temp_beta=0.0;
				temp_beta_1=10000;
				
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
						temp_beta = sign*q[temp][j];
                        
						if(temp_beta<temp_beta_1)
							temp_beta_1=temp_beta;
                        
					}
				}

				temp=mark_h[j][k];
				r[j][temp]=temp_afa*temp_beta_1;


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
				
				sum_Q=sum+r[temp][i];
                q[i][temp] =  sum + Lpi[i];
                
			}
		
			Q[i] = sum_Q + Lpi[i];
			
			if(Q[i]<0)            //decision
				decoded_bit[i]=1;
			else
				decoded_bit[i]=0;

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
	

	free(Q);
	free(q);
	free(r);
	free(Lpi);
	free(mark_v);
	free(mark_h);
	free(wc_per_column);
	free(wr_per_row);
	
	
}

