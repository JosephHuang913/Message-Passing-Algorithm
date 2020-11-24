// WiMAX 802.16e LDPC Matrix XOR Multiply 
// Editor: Shih Yun Yi
// Date: 2008/10
#include "stdafx.h"
# include <stdio.h>
# include <stdlib.h>
# include <iostream>
using namespace std;


void MatrixXORMultiply (int **A ,int **B,int **C,int rA,int cA,int cB)
{

	// Claculate A*B=C
	// rA is the row of matrix A
	// cA is the colum of matrix A  cA=rB
	// cB is the column of matrix B

	int i,j,k,sum;
	for (i=0;i<rA;i++)
	{
		for (j=0;j<cB;j++)
		{
			sum=0;
			for (k=0;k<cA;k++)
			{
				sum=sum+(A[i][k]*B[k][j]);
				sum=(int)sum%2;
			}
			C[i][j]=sum;
		}
	}
}