#ifndef __MATRIX_H__
#define __MATRIX_h__
#include <assert.h>  // Defines the assert function.

class Matrix {

public:

// Default Constructor. Creates a 1 by 1 matrix; sets value to zero. 
Matrix () {
  nRow_ = 1; nCol_ = 1;
  data_ = new double [1];  // Allocate memory
  set(0.0);                // Set value of data_[0] to 0.0
}

// Regular Constructor. Creates an nR by nC matrix; sets values to zero.
// If number of columns is not specified, it is set to 1.
Matrix(int nR, int nC = 1) {
  assert(nR > 0 && nC > 0);    // Check that nC and nR both > 0.
  nRow_ = nR; nCol_ = nC;
  data_ = new double [nR*nC];  // Allocate memory
  assert(data_ != 0);          // Check that memory was allocated
  set(0.0);                    // Set values of data_[] to 0.0
}

// Copy Constructor.
// Used when a copy of an object is produced 
// (e.g., passing to a function by value)
Matrix(const Matrix& mat) {
  this->copy(mat);   // Call private copy function.
}

// Destructor. Called when a Matrix object goes out of scope.
~Matrix() {
  delete [] data_;   // Release allocated memory
}

// Assignment operator function.
// Overloads the equal sign operator to work with
// Matrix objects.
Matrix& operator=(const Matrix& mat) {
  if( this == &mat ) return *this;  // If two sides equal, do nothing.
  delete [] data_;                  // Delete data on left hand side
  this->copy(mat);                  // Copy right hand side to l.h.s.
  return *this;
}

// Simple "get" functions. Return number of rows or columns.
int nRow() const { return nRow_; }
int nCol() const { return nCol_; }

// Parenthesis operator function.
// Allows access to values of Matrix via (i,j) pair.
// Example: a(1,1) = 2*b(2,3); 
// If column is unspecified, take as 1.
double& operator() (int i, int j = 1) {
  assert(i > 0 && i <= nRow_);          // Bounds checking for rows
  assert(j > 0 && j <= nCol_);          // Bounds checking for columns
  return data_[ nCol_*(i-1) + (j-1) ];  // Access appropriate value
}

// Parenthesis operator function (const version).
const double& operator() (int i, int j = 1) const{
  assert(i > 0 && i <= nRow_);          // Bounds checking for rows
  assert(j > 0 && j <= nCol_);          // Bounds checking for columns
  return data_[ nCol_*(i-1) + (j-1) ];  // Access appropriate value
}

// Set function. Sets all elements of a matrix to a given value.
void set(double value) {
  int i, iData = nRow_*nCol_;
  for( i=0; i<iData; i++ )
    data_[i] = value;
}

//*********************************************************************
private:

// Matrix data.
int nRow_, nCol_;  // Number of rows, columns
double* data_;     // Pointer used to allocate memory for data.

// Private copy function.
// Copies values from one Matrix object to another.
void copy(const Matrix& mat) {
  nRow_ = mat.nRow_;
  nCol_ = mat.nCol_;
  int i, iData = nRow_*nCol_;
  data_ = new double [iData];
  for(i = 0; i<iData; i++ )
    data_[i] = mat.data_[i];
}

}; // Class Matrix


// Compute inverse of complex matrix
void cinv( Matrix RealA, Matrix ImagA, 
			 Matrix& RealAinv, Matrix& ImagAinv ) 
// Inputs
//   RealA  -    Real part of matrix A (N by N)
//   ImagA  -    Imaginary part of matrix A (N by N)
// Outputs
//   RealAinv  - Real part of inverse of matrix A (N by N)
//   ImagAinv  - Imaginary part of A inverse (N by N)
{

  int N = RealA.nRow();
  assert( N == RealA.nCol() && N == ImagA.nRow() 
	                        && N == ImagA.nCol());
    RealAinv = RealA; // Copy matrices to ensure they are same size
  ImagAinv = ImagA;
  
  int i, j, k;
  Matrix scale(N);	 // Scale factor
  int *index;  index = new int [N+1];

  //* Matrix B is initialized to the identity matrix
  Matrix RealB(N,N), ImagB(N,N);
  RealB.set(0.0);  ImagB.set(0.0);
  for( i=1; i<=N; i++ )
    RealB(i,i) = 1.0;

  //* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
  for( i=1; i<=N; i++ ) {
    index[i] = i;			  // Initialize row index list
    double scaleMax = 0.;
    for( j=1; j<=N; j++ ) {
	  double MagA = RealA(i,j)*RealA(i,j) + ImagA(i,j)*ImagA(i,j);
      scaleMax = (scaleMax > MagA) ? scaleMax : MagA;
    }
	scale(i) = scaleMax;
  }

  //* Loop over rows k = 1, ..., (N-1)
  for( k=1; k<=N-1; k++ ) {
	//* Select pivot row from max( |a(j,k)/s(j)| )
    double ratiomax = 0.0;
	int jPivot = k;
    for( i=k; i<=N; i++ ) {
	  double MagA = RealA(index[i],k)*RealA(index[i],k) + 
		            ImagA(index[i],k)*ImagA(index[i],k);
      double ratio = MagA/scale(index[i]);
      if( ratio > ratiomax ) {
        jPivot=i;
        ratiomax = ratio;
      }
    }
	//* Perform pivoting using row index list
	int indexJ = index[k];
	if( jPivot != k ) {	          // Pivot
      indexJ = index[jPivot];
      index[jPivot] = index[k];   // Swap index jPivot and k
      index[k] = indexJ;
	}
	//* Perform forward elimination
    for( i=k+1; i<=N; i++ ) {
	  double denom = RealA(indexJ,k)*RealA(indexJ,k) 
		           + ImagA(indexJ,k)*ImagA(indexJ,k);
      double RealCoeff = (RealA(index[i],k)*RealA(indexJ,k)
		               + ImagA(index[i],k)*ImagA(indexJ,k))/denom;
      double ImagCoeff = (ImagA(index[i],k)*RealA(indexJ,k)
		               - RealA(index[i],k)*ImagA(indexJ,k))/denom;
      for( j=k+1; j<=N; j++ ) {
        RealA(index[i],j) -= RealCoeff*RealA(indexJ,j)
		                   - ImagCoeff*ImagA(indexJ,j);
        ImagA(index[i],j) -= RealCoeff*ImagA(indexJ,j)
		                   + ImagCoeff*RealA(indexJ,j);
      }
	  RealA(index[i],k) = RealCoeff;
	  ImagA(index[i],k) = ImagCoeff;
      for( j=1; j<=N; j++ ) {
        RealB(index[i],j) -= RealA(index[i],k)*RealB(indexJ,j)
			               - ImagA(index[i],k)*ImagB(indexJ,j);
        ImagB(index[i],j) -= RealA(index[i],k)*ImagB(indexJ,j)
			               + ImagA(index[i],k)*RealB(indexJ,j);
	  }
    }
  }

  //* Perform backsubstitution
  for( k=1; k<=N; k++ ) {
	double denom = RealA(index[N],N)*RealA(index[N],N) 
		         + ImagA(index[N],N)*ImagA(index[N],N);
    RealAinv(N,k) = (RealB(index[N],k)*RealA(index[N],N) 
		          + ImagB(index[N],k)*ImagA(index[N],N))/denom;
    ImagAinv(N,k) = (ImagB(index[N],k)*RealA(index[N],N) 
		          - RealB(index[N],k)*ImagA(index[N],N))/denom;
    for( i=N-1; i>=1; i--) {
      double RealSum = RealB(index[i],k);
      double ImagSum = ImagB(index[i],k);
      for( j=i+1; j<=N; j++ ) {
        RealSum -= RealA(index[i],j)*RealAinv(j,k)
			     - ImagA(index[i],j)*ImagAinv(j,k);
        ImagSum -= RealA(index[i],j)*ImagAinv(j,k)
			     + ImagA(index[i],j)*RealAinv(j,k);
	  }
	  double denom = RealA(index[i],i)*RealA(index[i],i) 
		           + ImagA(index[i],i)*ImagA(index[i],i);
      RealAinv(i,k) = (RealSum*RealA(index[i],i) 
		            + ImagSum*ImagA(index[i],i))/denom;
      ImagAinv(i,k) = (ImagSum*RealA(index[i],i) 
		            - RealSum*ImagA(index[i],i))/denom;
    }
  }

  delete [] index;	// Release allocated memory
}

#endif