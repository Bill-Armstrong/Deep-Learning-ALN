// ALN Library
// Copyright (C) 2018 William W. Armstrong.
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// Version 3 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
// 
// For further information contact 
// William W. Armstrong
// 3624 - 108 Street NW
// Edmonton, Alberta, Canada  T6J 1B4

// Some matrix routines are taken from Eigen and are under the Mozilla Public License MPL2

// calccovariance.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif
#include <aln.h>
#include <alnpriv.h>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

#define TOL 1.0e-12

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
// helpers for filling and dumping matrices 
//void fillRMA(MatrixXd &X, double* adblX);
//void ALNAPI dumpRMA(const MatrixXd &X, double* adblX);
//void ALNAPI fillCMA(MatrixXd &X, double* adblX);
//void ALNAPI dumpCMA(const MatrixXd &X, double* adblX);

void fillRMA(MatrixXd &X, double* adblX)
{
	// The array adblX is organized in row-major order
	unsigned int nRows, nr, nCols, nc;
	nRows = X.rows();
	nCols = X.cols();
	for(nr = 0; nr < nRows; ++nr)
	{
		for(nc = 0; nc < nCols; ++nc)
		{
			X(nr,nc) = adblX[nr*nCols + nc];
		}
	}
}

void ALNAPI dumpRMA(const MatrixXd &X, double* adblX)
{
	// The array adblX is organized in row-major order
	int nRows, sr, nCols, sc;
	nRows = X.rows();
	nCols = X.cols();
	for(sr = 0; sr < nRows; ++sr)
	{
		for(sc = 0; sc < nCols; ++sc)
		{
			adblX[sr*nCols + sc] = X(sr,sc) ;
		}
	}
}

void ALNAPI fillCMA(MatrixXd& X, double* adblX)
{
	// The array adblX is organized in column-major order
	int nRows, nr, nCols, nc;
	nRows = X.rows();
	nCols = X.cols();
	for(nc = 0; nc < nCols; ++nc)
	{
		for(nr = 0; nr < nRows; ++nr)
		{
			X(nr, nc) = adblX[nc*nRows + nr];
		}
	}
}

void ALNAPI dumpCMA(const MatrixXd &X, double* adblX)
{
	// The array adblX is organized in column-major order
	int nRows, nr, nCols, nc;
	nRows = X.rows();
	nCols = X.cols();
	for(nc = 0; nc < nCols; ++nc)
	{
		for(nr = 0; nr < nRows; ++nr)
		{
			adblX[nc*nRows + nr] = X(nr, nc);
		}
	}
}


BOOL ALNAPI CalcCovariance(int nCols, // number of input vars = number of ALN inputs - 1
                    int nRows,        // number of rows, i.e. input vectors for a given LFN
                    double* adblX,    // RMA based input vectors (nRows * nCols)
                    double* adblY,    // RMA based result vector (nRows)
                    double* adblC,    // RMA covariance matrix as array (nCols * nCols), is returned
                    double* adblA,    // RMA fitted parameter vector (nCols)
                    double* adblS,    // RMA std dev vector (nRows) (a weighting on points) 
                    double* adblU,    // RMA U matrix (nRows * nRows)
                    double* adblV,    // RMA V matrix (nCols * nCols)
                    double* adblW,    // RMA W array (singvals = min(nRows,nCols))containing singular values
                    double& dblChiSq) // chi square of fit
{
	// implementation based on SVD in Eigen
	// note that matrices below are stored in row-major order (RMA)
  ASSERT(adblX && adblY && adblC && adblA && adblS && adblU && adblV && adblW);
  BOOL bSuccess = TRUE;
	MatrixXd X(nRows,nCols);
	fillRMA(X,adblX);
	MatrixXd U(nRows,nRows);
	MatrixXd V(nCols,nCols);
	int singvals = std::min(nRows, nCols);
	VectorXd  s(singvals);  // the singular values s
	VectorXd  sinv(nCols);  
	// inverses of the singular values s, but inverse of a very small s is 0 
	// We should solve for fitting parameters A then compare it to the ALNresult
	// We get A = V W^(-1)U^T Y.  This works for square X but here
	// we would get UWV^TA = UWV^TVW^(-1)U^T Y =Y except for the zeroed  W entries

	try
	{
		double smax, tmp, thresh;
		// in a try/catch block to free b in case exception occurs
		// now we do the SVD on X
		JacobiSVD<MatrixXd> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV();
		s = svd.singularValues();
		// put result of the SVD into the arrays in the call
		dumpRMA(U, adblU);
		dumpRMA(V, adblV);
		// find the maximal singular value
		smax = 0.0;
		for (int j = 0; j < singvals; j++) {
			if (s(j) > smax) smax = s(j);
		}
		// if any (non-negative) singular value is less than a threshold, set it to zero
		thresh = TOL*smax;
		for (int j = 0; j < singvals; j++){
				if (s(j) < thresh)	s(j) = 0.0;
				adblW[j] = s(j);
		}
		// Create the inverses, extended by 0 where necessary to make nCols
		for (int j = 0; j < nCols; j++) {
			if( j >= singvals || s(j) == 0.0){
				sinv(j)=0.0;  // surprising but the right way to do it
			}
			else
			{
				sinv(j) = 1.0/s(j);
			}
		}

		// compute the fitting parameters for rNows of output values
		// First make a vector out of Y
		VectorXd Y(nRows);
		for (int sr = 0; sr < nRows; sr++)
		{
			Y(sr) = adblY[sr]; // this Y is incorrect for the entries j = 0 to j = 4  Something is overwriting the Y array!!!!
		}
		MatrixXd UT(nRows, nRows);
		UT = U.transpose();  // This is the inverse of U
		VectorXd tempY(nRows);
		tempY = UT * Y; // This should be the output of the diagonal matrix in the middle of the svd
		VectorXd tempV(nCols); 
		// Instead of taking the inverse of the diagonal matrix and multiplying by tempY we do the components
		for (int sr = 0; sr < nCols; sr++){
			if(sr < singvals){
				tempV(sr) = sinv(sr) * tempY(sr);
			}
			else {
				tempV(sr) = 0;
			}
		}
		VectorXd A(nCols); // A is the solution of the least-squares problem, i.e. the ideal weights on the LFN
		A = V * tempV; // V is the inverse of V^T
		// return the values in A to the pointer structure
		for (int i = 0; i < nCols; i++)
		{
			adblA[i] = A(i);
		}
		// finding the measure of error of the fit by chi-squared 
		double dblChiSq=0.0;
		for (int i = 0; i < nRows; i++) 
		{
			double sum = 0.0;
			for (int j = 0; j < nCols; j++)
			{
				sum += X(i,j)* A(j);
			}
			tmp=(Y(i) - sum)/adblS[i]; // adblS[i] is a weight on input vector i, usually 1.0
			dblChiSq += tmp*tmp;
		}

    // calc covariance matrix C
		for(int k = 0; k < nCols; k++)
		{
			for(int i = 0; i <= k; i++)		// use symmetry
			{
				double sum = 0.0;
				for(int j = 0; j < nCols; j++)
				{
					sum += V(k,j) * V(i,j) * sinv(j) * sinv(j);
				}
				adblC[i*nCols+k] = adblC[k*nCols+i] = sum; 
			}
		}
  }

  catch(string SVDdecomposition)
  { 
    bSuccess = FALSE; // something failed!
  }  
  return bSuccess;
}
