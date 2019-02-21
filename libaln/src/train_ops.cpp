// ALN Library
// file train_ops.cpp
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

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// include files
#include <stdafx.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <alnpp.h>
#include <dtree.h>
#include <datafile.h>
#include <malloc.h>
#include ".\cmyaln.h" 
#include "alnextern.h"
#include "alnintern.h"
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
double dist(double*, double*); // calculates the distance between domain points
#define PrintInterval 25
//defines used to set up TRfile and NVfile
#define LINEAR_REGRESSION 0
#define NOISE_VARIANCE 1
#define APPROXIMATION 2
#define BAGGING 3

// We use dblRespTotal in two ways and the following definition helps.
#define DBLNOISEVARIANCE dblRespTotal

// files used in training operations
CDataFile TRfile; // Training data, changed for different purposes.
CDataFile NVfile; // Becomes a file of noise variance samples.

//routines
void ALNAPI doLinearRegression(); // Determines an upper bound on error, and provides a start for other training.
void ALNAPI createNoiseVarianceFile(); // This creates the noise variance file NVfile.
void ALNAPI approximate(); // Actually does training avoiding overtraining using samples in NVfile.
void ALNAPI outputTrainingResults();
void ALNAPI trainAverage(); // Takes several ALNs created in approximate() and creates an ALN of their average
void ALNAPI constructDTREE(int nMaxDepth); // Takes the average ALN and turns it into a DTREE
void ALNAPI cleanup(); // Destroys ALNs etc.
void fillvector(double * adblX, CMyAln* paln); // Sends a vector to training from a file, online or averaging.
void ALNAPI createTR_NVfiles(int nChoose);
void ALNAPI trainNoiseVarianceALN(); // Trains on the samples in NVfile to create NV_ALN
void createSamples(CMyAln* pNV_ALN); // Creates smoothed noise variance samples using NV_ALN
void prepareQuickStart(CMyAln* pALN); // This sets up ALN training with values already determined by Linear Regression
double scalarProduct(double* adblA, double* adblB);
// ALN pointers

static CMyAln* pALN = NULL; // declares a pointer to an ALN used in linear regression
static CMyAln* pNV_ALN = NULL; // an ALN that smooths out the noise variance samples to get the NV function
static CMyAln** apALN = NULL;  // an array of pointers to ALNs used in approximate()
static CMyAln* pAvgALN = NULL;      // an ALN representing the bagged average of several ALNs trained on the TVfile with different random numbers

// Some global variables
double dblMinRMSE = 0; // stops training when the training error is smaller than this
double dblLearnRate = 0.2;  // roughly, 0.2 corrects most of the error for if we make 15 passes through TRfile
int nMaxEpochs = 10; // if the learnrate is 0.2, then one will need 5 or 10 roughly to almost correct the errors
long nRowsTR; // the number of rows in the current training set loaded into TRfile
long nRowsNV; // the number of rows in the noise variance file.  When approximation starts, this should be nRowsTV
double dblLimit = -1  ;// A negative value splits pieces based on an F test, otherwise they split if training MSE < dblLimit.
BOOL bALNgrowable = TRUE; // FALSE for linear regression, TRUE when the ALN has to grow.
BOOL bTrainNV_ALN = FALSE; // Controls setup of training for to fit noise samples in NVfile
BOOL bStopTraining = FALSE; // Set to TRUE and becomes FALSE if any (active) linear piece needs training
int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO; // used with callbacks
double * adblX = NULL;
double noiseSampleSum = 0;
double noiseSampleMax = 0;
double noiseSampleMin = DBL_MAX; // This allows us to avoid negative noise after training.
int Nearby; // The number of nearby samples; set near start of createNoiseVarianceFile().

using namespace std;

void ALNAPI doLinearRegression() // routine
{
  fprintf(fpProtocol,"\n****** Linear regression begins: it gains useful information for further training *****\n");
	fflush(fpProtocol); 
	// If something goes wrong, we can perhaps see it in the protocol file if we flush before the crash.
	// This iterative algorithm is not an accepted method in general, but works well here.
	// The reason for using it for ALN approximation is that the linear regression
	// problem for a linear piece (leaf node = LFN) is constantly changing. A piece gains or loses input vectors (samples).
	// for which it determines the output value via the ALN tree of max and min operators.
	// This constantly changes the samples which are to be least-squares fitted for a given linear piece.
	// Linear regression helps get the centroid and weights of a linear piece in neighborhood of
	// good values to start other training later.
	// Set up a sample buffer.
	adblX = (double *)malloc((nDim) * sizeof(double));
	// Set up the ALN
	pALN = new CMyAln;
	if (!(pALN->Create(nDim, nDim-1))) // nDim is the number of inputs to the ALN plus one for
		//the ALN output. The default output variable is nDim-1.
	{
		fprintf(fpProtocol,"Stopping: linear regression ALN creation failed!\n");
    fflush(fpProtocol);
		exit(0);
	}
	// NB The region concept has not been completed.  It allows the user to impose constraints
	// e.g. on slopes(weights) which differ in different parts of the domain.  All we have is region 0 now.
	bALNgrowable = FALSE; // The ALN consists of one leaf node LFN for linear regression.
	bTrainNV_ALN = FALSE; // TRUE only during training of the noise variance ALN.
	bTrainingAverage = FALSE; // Switch to tell fillvector whether get a training vector or compute an average.
	prepareQuickStart(pALN);
	nMaxEpochs = 20; // nMaxEpochs is the number of passes through the data (epochs) for each call to Train.
	// If the tree is growable, this is epochs before each splitting.
	// For linear regression, it is just a number of epochs.
	// A learning rate of 0.2 means that the error for a training sample is reduced by 20% each adapt.
	// Reducing the error on one training point may make the error of a different sample greater.
	// The theoretical best would be that 0.2 in 20 epochs reduces the error to 1.15% of what it was.
	dblMinRMSE = 0; // Don't stop early because of low training error.
	dblLearnRate = 0.2;  // This rate seems ok for linear regression.
	dblLimit = 0.001 * adblStdevVar[nDim - 1]; // Training MSE below which splitting stops is 0.1% st. dev. of the output.
	(pALN->GetRegion(0))->dblSmoothEpsilon = 0.0; // No smoothing because there is no splitting.
	// Set up the data
	createTR_NVfiles(LINEAR_REGRESSION);
	int nRowsTR = TRfile.RowCount();	// nEpochsize gives the number of training samples. Later nRowsNV=nRowsTR.
	int nColumns = TRfile.ColumnCount(); // This is always nDim for training.
	ASSERT(nColumns == nDim);
	const double* adblData = TRfile.GetDataPtr(); // This is where training gets samples.
	// The third parameter in thefollowing could also set to NULL instead of adblData.
	// Then, instead of using FillInputVector(), the program uses fillvector()
	// for setting up the input vectors to the ALN.  fillvector() allows the
	// system to choose training vectors more flexibly (even online with proper programming).
	// The advantage of giving the pointer adblData instead of NULL is that training permutes
	// the order of the samples and goes through all samples exactly once per epoch.
	pALN->SetDataInfo(nRowsTR, nDim, adblData, NULL);
	int nIterations = 3;
	// TRAIN FOR LINEAR REGRESSION   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	// The reason for iterations is so that we can monitor progress in the ... TrainProtocol.txt file,
	// and set new parameters for subsequent training.

	for (int iteration = 0; iteration < nIterations; iteration++)  // experimental
	{
		if (!pALN->Train(nMaxEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask))
		{
			fprintf(fpProtocol, "Linear regression training failed!\n");
			fflush(fpProtocol);
			exit(0);
		}
		fprintf(fpProtocol, "Iteration %d ",iteration);
		if (iteration == (nIterations * 1) / 5) dblLearnRate = 0.15;
		if (iteration == (nIterations * 2) / 5) dblLearnRate = 0.10;
		if (iteration == (nIterations * 3) / 5) dblLearnRate = 0.05;
		if (iteration == (nIterations * 4) / 5) dblLearnRate = 0.01;
	}
	fprintf(fpProtocol, "Linear regression training succeeded!\n");
	fflush(fpProtocol);
  // We should now have a good linear regression fit and we harvest the important results.
  // Find the weights on the linear piece using an evaluation at the 0 vector (could be any other place!).
	ALNNODE* pActiveLFN;
  for(int m=0; m < nDim; m++)
	{
		adblX[m] = 0;
	}
	double dummy = pALN->QuickEval(adblX, &pActiveLFN); // this finds a pointer to the linear piece
	fprintf(fpProtocol,"Linear regression weights on ALN\n");
	for (int m = 0; m < nDim-1; m++)
	{
    if(nLag[m] == 0) // any data we are looking at can have nonzero lags
    {
			// Note that adblW is stored in a different way, the weight on axis 0 is in adblW[1] etc.
  		fprintf(fpProtocol,"Linear regression weight on %s is %f centroid is %f\n",\
			varname[nInputCol[m]] , ((pActiveLFN)->DATA.LFN.adblW)[m+1], ((pActiveLFN)->DATA.LFN.adblC)[m]); 
		}
    else
    {
  		fprintf(fpProtocol,"Linear regression weight on ALN input %s@lag%d is %f centroid is %f\n",\
			varname[nInputCol[m]], nLag[m] , ((pActiveLFN)->DATA.LFN.adblW)[m+1], ((pActiveLFN)->DATA.LFN.adblC)[m]); 
    }
		adblLRC[m] = ((pActiveLFN)->DATA.LFN.adblC)[m]; // centroid of linear piece for axis m
		adblLRW[m+1] = ((pActiveLFN)->DATA.LFN.adblW)[m+1]; // linear regression weight for axis m
  }
	fprintf(fpProtocol,"Value of the linear piece at the origin is %f; the value at the centroid is %f\n",\
		((pActiveLFN)->DATA.LFN.adblW)[0], ((pActiveLFN)->DATA.LFN.adblC)[nDim - 1] ); 
	fflush(fpProtocol);
	// Put the negative sum of weighted centroids into the 0 position
	// compress the weighted centroid info into W[0], i.e. take the centroids out of all terms like
	// ... + weight_m*(x_m-centroid_m) + ... to have just one number in W[0]
	adblLRW[0] = adblLRC[nDim -1] = ((pActiveLFN)->DATA.LFN.adblC)[nDim -1];
  for (int m = 0; m < nDim - 1; m++)
  {
    adblLRW[0] -= adblLRW[m+1] * adblLRC[m];
  }
	adblLRW[nDim] = -1.0;  // the -1 weight on the ALN output may make a difference somewhere, so we set it here.
	// The idea of weight = -1 is that the hyperplane equation can be written 
	// a*x + b*y -1*z = 0 like an equation instead of z = a*x + b*y like a function.
	// This representation is good if we want to invert the ALN, i.e. get x as a function of y and z.
	// How to do it? Just multiply the first equation by -1/a in all linear pieces. The weight on the output becomes -1.0
	// This is not fully implemented yet.  The output of the ALN is now always the highest index ALN variable nDim -1.

	double desired, predict, se;
	int nCount; 
	nCount = 0;
	for (int j = 0; j < nRowsTR; j++)
	{
		for (int i = 0; i < nDim; i++)
		{
			adblX[i] = TRfile.GetAt(j, i, 0);
		}
		desired = adblX[nDim - 1]; // get the desired result
		adblX[nDim - 1] = 0; // not used in evaluation by QuickEval
		predict = pALN->QuickEval(adblX, &pActiveLFN);
		se = (predict - desired) * (predict - desired);
		nCount++;
		dblLinRegErr +=se;
	} // end loop over  the training file
	dblLinRegErr = sqrt(dblLinRegErr / nCount); 
	fprintf(fpProtocol, "Root Mean Square Linear Regression Training Error is %f \n", dblLinRegErr);
	fprintf(fpProtocol, "The number of data points used in determining linear regression error is %d \n",nCount);
	fflush(fpProtocol);
	pALN->Destroy();
	free(adblX);
	// We are finished with that ALN and have destroyed it
	fprintf(fpProtocol,"Linear regression complete\n");
	fflush(fpProtocol);
}

void ALNAPI createNoiseVarianceFile()
{
	fprintf(fpProtocol, "\n ********* Begin Creation of Noise Variance File ********\n");
	Nearby = 8 * nDim; // This may have to be adjusted. It should be set to the largest
	// number of samples nearby some sample where the ideal function is almost flat.
	fprintf(fpProtocol, "\nUsed to determine noise variance are %d nearby samples.\n", Nearby);
	double* aXcentral = NULL;
	double* aXnearby = NULL;
	double* aXmid = NULL;
	struct disc // each sample (X,y) index i in NVfile has a disc struct at index i
	{
		long*   aXX; // These are indices of samples in NVfile close to (X,y)
		double* aDD; // These are the corresponding distances of those samples to (X,y)
		int maxloc;  // Index above which data in aXX and aDD are meaningless
	};
	disc* adisc = NULL;
	// Example: domain 8 dimensional, samples 9 components, nDim = 9.
	// Each sample i in turn takes on the role of central sample (X,y). 
	// We go through all samples j to find a set of Nearby closest samples to sample i as follows:
	// If j is the index i of the central sample, we go on to the next sample.
	// We compute the distance ds of sample j to (X,y).
	// If there are fewer than Nearby nearby samples, we sort sample j into its proper place in increasing distance
	// If there are Nearby nearby samples, and ds >= the maximal distance of a nearby sample we go on to the next j
	// Otherwise, we get rid of the sample at maximum distance and sort j into its proper place in increasing distance.
	aXcentral = (double*)malloc((nDim - 1) * sizeof(double)); // these hold domain points
	aXnearby = (double*)malloc((nDim - 1) * sizeof(double));
	aXmid = (double*)malloc((nDim - 1) * sizeof(double));
	adisc = (disc*)malloc(nRowsNV * sizeof(disc));
	//We use the NVfile here to set up samples in TRfile to train the NV_ALN.
	long i, j, k;
	double ds;
	for (i = 0; i < nRowsNV; i++) // i is the central sample its struct will be worked on.
	{
		// First we have to create the struct for i, we allow up to Nearby samples closest to sample i.
		adisc[i].aXX = (long*)malloc(Nearby * sizeof(long));
		adisc[i].aDD = (double*)malloc(Nearby * sizeof(double));
		adisc[i].maxloc = -1; // All aXX[k], aDD[k] with k > maxloc are always
		// meaningless, so at initialization all the aXX and aDD data are meaningless
		// So, of course, we don't have to initialize the aXX and aDD arrays,
		// but we do it anyway.
		for (k = 0; k < Nearby; k++)
		{
			adisc[i].aXX[k] = 0;
			adisc[i].aDD[k] = 0;
		}
		// Get the central sample aXcentral with index i.
		for (k = 0; k < nDim - 1; k++)
		{
			aXcentral[k] = NVfile.GetAt(i, k, 0); // just the domain components
		}
		// We have to go through ALL samples j now to get Nearby closest samples to i.
		// This can be made faster than O(n^2), maybe O(n log n).
		for (j = 0; j < nRowsNV; j++)
		{
			if (j != i) // sample i is not included in the nearby ones
			{
				// Get the domain point for the j-th sample into aXnearby
				for (k = 0; k < nDim - 1; k++)
				{
					aXnearby[k] = NVfile.GetAt(j, k, 0);
				}
				ds = dist(aXcentral, aXnearby);
				int jtemp;
				double dstemp;
				int imaxloc = adisc[i].maxloc;
				if (imaxloc < Nearby - 1) 
				{
					// add the new sample at the next available place
					++imaxloc;
					adisc[i].aXX[imaxloc] = j;
					adisc[i].aDD[imaxloc] = ds;
					adisc[i].maxloc = imaxloc;
				}
				else if (ds < adisc[i].aDD[imaxloc])
				{
					ASSERT(imaxloc == Nearby - 1);
					// imaxloc = Nearby -1 and since ds is smaller than for that sample,
					// we replace the latter with the j sample.
					// We assume the nearby samples are sorted; that imaxloc
					// is the index of a sample at maximum distance
					// insert j at imaxloc
					adisc[i].aXX[imaxloc] = j;
					adisc[i].aDD[imaxloc] = ds;
				}
				// now sort the new sample into place in order of increasing ds
				k = imaxloc;
				while ((k > 0) && (adisc[i].aDD[k] < adisc[i].aDD[k - 1]))
				{
					// interchange the k and k - 1 samples
					jtemp = adisc[i].aXX[k];
					dstemp = adisc[i].aDD[k];
					adisc[i].aXX[k] = adisc[i].aXX[k - 1];
					adisc[i].aDD[k] = adisc[i].aDD[k - 1];
					adisc[i].aXX[k - 1] = jtemp;
					adisc[i].aDD[k - 1] = dstemp;
					k--;
				}
			} // drop samples with j == i
		} // end of j loop
	} // end loop over i

	// We now have samples indexed 0 to Nearby - 1 closest to each sample index i
	// in order of increasing distance from the central sample.
	// The following computes the error sum of squares ESS
	long index;
	double NV;
	fprintf(fpProtocol, "nDim = %d, Nearby = %d\n",	nDim, Nearby);
	// the following is adapted from Ken Cogger's write-up.
	nRowsNV = NVfile.RowCount();
	MatrixXd X(Nearby+1, nDim);  // The X matrix as used in Ken's writeup.
	MatrixXd XtX(nDim, nDim); // X transposed times X
	VectorXd y(Nearby+1); // a column of Nearby sample values
	VectorXd Xty(nDim); // X transposed times y
	VectorXd betaHat(nDim); // a helper vector
	VectorXd yy(Nearby+1); // a helper meaning y - yHat = y - X * betaHat
	for (i = 0; i < nRowsNV; i++)
	{
		double val = 0;
		double ESS = 0; // The error sum of squares over the central + Nearby samples
		for (k = 0; k < nDim - 1; k++)
		{
			aXmid[k] = 0; //aXmid will be in the middle of the Nearby domain points.
		}
		//Get the nearby sample domains and values
		for (j = 0; j < Nearby; j++)
		{
			index = adisc[i].aXX[j];
			for (k = 0; k < nDim - 1; k++)
			{
				val = NVfile.GetAt(index, k, 0);
				X(j, k) = val; // put index sample's domain vector into the j-th row of X
				aXmid[k] += val;
			}
			X(j, nDim - 1) = 1.0;
			y(j) = NVfile.GetAt(index, nDim - 1, 0); // put sample value into y(j)
		}
		// Now include the central sample
		for (k = 0; k < nDim - 1; k++)
		{
			val = NVfile.GetAt(i, k, 0);
			X(Nearby, k) = val; // put index sample's domain vector into the j-th row of X
			aXmid[k] += val;
			aXmid[k] /= Nearby + 1;
		}
		X(Nearby, nDim - 1) = 1.0;
		y(Nearby) = NVfile.GetAt(i, nDim - 1, 0);
		XtX = X.transpose() * X;
		Xty = X.transpose() * y;
		// Now we have to, in effect, invert XTX to compute yHat = X * betaHat where
		// betaHat = inverse of (XTX) * Xty , but Eigen doesn't invert matrices!
		// Instead we solve  XtX * betaHat = Xty:
		betaHat = XtX.colPivHouseholderQr().solve(Xty);
		// The error sum of squares, ess, is (y - yHat).transpose() * (y - yHat),
		yy = (y - X * betaHat);
		ESS = yy.transpose() * yy;
		NV = ESS /(Nearby + 1 - nDim); // ESS is the sum of squared errors.
		// An ALN's mean squared error should be close to NV at this place.
		// The noise variance is probably ESS/(Nearby + 1 - nDim) since the
		// fitted linear piece removes nDim degrees of freedom)
		noiseSampleSum += NV;
		if (NV < noiseSampleMin)noiseSampleMin = NV;
		if (NV > noiseSampleMax)noiseSampleMax = NV;
		// We set up TRfile for training on the NV samples
		for (k = 0; k < nDim - 1; k++)
		{
			TRfile.SetAt(i, k, aXmid[k], 0);
		}
		TRfile.SetAt(i, nDim - 1, NV , 0);
	} // end i loop over nRowsNV
	fprintf(fpProtocol, "The average, minimum and maximum Noise Variance are %f %f %f \n",
		noiseSampleSum/nRowsNV, noiseSampleMin, noiseSampleMax);
	TRfile.Write("DiagnoseTRfileBeforeNV_ALN_train.txt");
	for (i =0; i < nRowsNV ; i++)
	{
		// cleanup
		free(adisc[i].aXX);
		free(adisc[i].aDD);
	}
	free(aXcentral);
	free(aXnearby);
	free(adisc);
	trainNoiseVarianceALN();
}
void ALNAPI trainNoiseVarianceALN()
{
	fprintf(fpProtocol, "\n**************Training Noise Variance ALN begins ********\n");
	fflush(fpProtocol);
	pNV_ALN = (CMyAln*)malloc(sizeof(CMyAln*));
	pNV_ALN = new CMyAln; // NULL initialized ALN
	// Set up the sample buffer.
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	// Set up the ALN
	if (!pNV_ALN->Create(nDim, nDim - 1))
	{
		fprintf(fpProtocol, "ALN creation failed!\n");
		fflush(fpProtocol);
		exit(0);
	}
	// Make it growable
	if (!pNV_ALN->SetGrowable(pNV_ALN->GetTree()))
	{
		fprintf(fpProtocol, "Setting NV_ALN growable failed!\n");
		fflush(fpProtocol);
		exit(0);
	}
	bALNgrowable = TRUE;
	bTrainingAverage = FALSE;
	bTrainNV_ALN = TRUE;
	(pNV_ALN->GetRegion(0))->dblSmoothEpsilon = 0.0;// There is room for experimentation.
	nMaxEpochs = 20; // Splitting of linear pieces occurs at the end of these epochs.
	dblMinRMSE = dblLinRegErr * 1e-6; // Stops training when the error is tiny 
	nNumberLFNs = 1;  // initialize at 1
	nRowsNV = NVfile.RowCount();
	ASSERT(nRowsNV == nRowsTV);
	// We now put the smoothed TRfile content into NVfile to use in an F-test
	createTR_NVfiles(NOISE_VARIANCE); // We shouldn't have to always train on TRfile because of splitops
	const double* adblData = TRfile.GetDataPtr();
	pNV_ALN->SetDataInfo(nRowsTR, nDim, adblData, NULL); // Not possible yet to train on NVfile
	dblLimit = -1.0;  // Use a special F-test where NVfile is used differently in split_ops
	fprintf(fpProtocol, "----------  Training NV_ALN  ------------------\n");
	fflush(fpProtocol);
	dblLearnRate = 0.2; // the noise of the noise is large
	fprintf(fpProtocol, "The smoothing for training NV_ALN is %f, initial learning rate %f\n",
		(pNV_ALN->GetRegion(0))->dblSmoothEpsilon, dblLearnRate);
	for (int iteration = 0; iteration < 30; iteration++)
	{
		fprintf(fpProtocol, "Iteration %d ", iteration);
		// TRAIN NOISE VARIANCE ALN   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
		if (!pNV_ALN->Train(nMaxEpochs, dblMinRMSE, dblLearnRate, bJitter, nNotifyMask))
		{
			fprintf(fpProtocol, "Training failed!\n");
			fflush(fpProtocol);
			exit(0);
		}
		if (bStopTraining == TRUE)
		{
			fprintf(fpProtocol, "This training stopped because all leaf nodes have stopped changing!\n");
			fprintf(fpProtocol, "\nTraining of NV_ALN completed at iteration %d \n", iteration);
			bStopTraining = FALSE;
			fflush(fpProtocol);
			break;
		}
		fflush(fpProtocol);
		if(iteration == 3) dblLearnRate = 0.15;
		if (iteration == 6) dblLearnRate = 0.125;
		if (iteration == 10) dblLearnRate = 0.10;
		if (iteration == 15) dblLearnRate = 0.075;
		if (iteration == 20) dblLearnRate = 0.05;
	} // end of loop of training iterations for NV_ALN
	free(adblX);
	bTrainNV_ALN = FALSE;
	createSamples(pNV_ALN); // Changes values in NVfile using NV_ALN
}

void ALNAPI approximate() // routine
{
	fprintf(fpProtocol, "\n**************Approximation with one or more ALNs begins ********\n");
	fflush(fpProtocol);
	int nalns = nALNs;  // The number of ALNs over which we average (for "bagging")
	createTR_NVfiles(APPROXIMATION);  // prepares for using the whole TVfile and the whole NVfile with noise variance samples for training and stopping
	fprintf(fpProtocol,"Training %d approximation ALNs starts with the goal of avoiding overtraining\n", nalns);
	fflush(fpProtocol);
	if(bJitter)
	{
		fprintf(fpProtocol, "Jitter is used during approximation\n");
	}
	else
	{
		fprintf(fpProtocol, "Jitter is not used during approximation\n");
	}
	fflush(fpProtocol);

	dblLimit = -1.0; // Negative to split leaf nodes according to an F test;
	// positive to split if training MSE > dblLimit.
	if (dblLimit <= 0)
	{
		fprintf(fpProtocol, "An F test is used to decide whether to split a piece depending on hit count. \n");
	}
	else
	{
		fprintf(fpProtocol, "A manually set limit, %f, is used to decide whether to split a piece. \n", dblLimit);
	}
	fflush(fpProtocol);
	// ************ SET UP THE ARRAY OF POINTERS TO ALNS FOR TRAINING **********
	if(bTrain)
	{
		apALN = (CMyAln**) malloc(nALNs * sizeof(CMyAln*));
	}
	fflush(fpProtocol);

	for (int n = 0; n < nalns; n++)
	{
		 apALN[n] = new CMyAln; // NULL initialized ALN
	}
	// Set up the sample buffer.
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	// Train all the ALNs
	for (int n = 0; n < nalns; n++)
	{
		// Set up the ALN, index n
		if (!(apALN[n]->Create(nDim, nDim-1)))
		{
	    fprintf(fpProtocol,"ALN creation failed!\n");
      fflush(fpProtocol);
			exit(0);
		}
		// Now make the tree growable
		if (!apALN[n]->SetGrowable(apALN[n]->GetTree()))		
		{
	    fprintf(fpProtocol,"Setting ALN %d growable failed!\n", n);
      fflush(fpProtocol);
      exit(0);
		}
		bALNgrowable = TRUE; // now the nodes can split
		// Set constraints on variables for ALN index n
		prepareQuickStart(apALN[n]);
		bTrainingAverage = FALSE;
		bTrainNV_ALN = FALSE;
		if(bClassify)
    {
     // TO DO
    }
		(apALN[n]->GetRegion(0))->dblSmoothEpsilon = 0.0;
		fprintf(fpProtocol, "The smoothing for training each approximation is %f\n", 0.0); 
		nMaxEpochs = 20;
		dblMinRMSE = dblLinRegErr * 1e-6; // Stops training when the error is tiny 
		dblLearnRate = 0.2;
		bStopTraining = FALSE; // Set TRUE in alntrain.cpp. 
		// Set FALSE by any piece needing more training. 
    nNumberLFNs = 1;  // initialize at 1
		// Set up the data
		// Tell the training algorithm the way to access the data using fillvector
		nRowsTR = TRfile.RowCount();
		ASSERT(nRowsTR == nRowsTV);
		const double* adblData = TRfile.GetDataPtr();
		apALN[n]->SetDataInfo(nRowsTR, nDim, adblData, NULL);
		fprintf(fpProtocol,"----------  Training approximation ALN %d ------------------\n",n);
		fflush(fpProtocol);
		for(int iteration = 0; iteration < 25; iteration++) // 30 iterations seems enough
		{
			fprintf(fpProtocol, "\nIteration %d ALN %d ", iteration, n);
			fflush(fpProtocol);

			// TRAIN ALNS WITHOUT OVERTRAINING   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
			if (!apALN[n]->Train(nMaxEpochs, dblMinRMSE, dblLearnRate, bJitter, nNotifyMask))
			{
			  fprintf(fpProtocol,"Training failed!\n");
				fflush(fpProtocol);
        exit(0);
			}
			if (bStopTraining == TRUE)
			{
				fprintf(fpProtocol, "\nTraining of approximation ALN %d completed at iteration %d \n", n, iteration);
				fprintf(fpProtocol, "This training stopped because all leaf nodes have stopped changing!\n");
				fflush(fpProtocol);
				break;
			}
			fprintf(fpProtocol, "Learning rate is %f\n", dblLearnRate);
			fflush(fpProtocol);
			if (iteration == 3) dblLearnRate = 0.15;
			if (iteration == 6) dblLearnRate = 0.125;
			if (iteration == 10) dblLearnRate = 0.10;
			if (iteration == 15) dblLearnRate = 0.075;
			if (iteration == 20) dblLearnRate = 0.05;
			} // end of loop of training interations over one ALN
  } // end of the loop for n = ALN index
	free(adblX);
	// we don't destroy the ALNs because they are needed for further work in reporting
}

void ALNAPI outputTrainingResults() // routine
{
	fprintf(fpProtocol, "\n**** Analyzing results of approximation begins ***\n");
	// all the ALNs have been trained, now report results
	int i, j, k, n;
	double desired, average, sum;
	int nalns; // lower case indicates a value on the stack
	nalns = nALNs;
	ALNNODE* pActiveLFN = NULL;
	// test the average of the ALNs against data in the TV set
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	double * adblWAcc = (double *)malloc((nDim) * sizeof(double));
	double * adblAbsWAcc = (double *)malloc((nDim) * sizeof(double));
	for (k = 0; k < nDim; k++)
	{
		adblWAcc[k] = 0; // zero the accumulator for the average weight
		adblAbsWAcc[k] = 0; // zero the accumulator for the average absolute weight
	}
	double se = 0; // square error accumulator

	int	nClassError = 0;  // for classification problems
	for (j = 0; j < nRowsTV; j++)
	{
		sum = 0;

		for (n = 0; n < nalns; n++)
		{
			for (i = 0; i < nDim; i++)
			{
				adblX[i] = TVfile.GetAt(j, i, 0);
			}
			double dblValue = apALN[n]->QuickEval(adblX, &pActiveLFN);
			sum += dblValue;
			for (int k = 0; k < nDim; k++)
			{
				adblWAcc[k] += ((pActiveLFN)->DATA.LFN.adblW)[k + 1]; //the adblW vector has the bias in it
																												// so the components are shifted
				adblAbsWAcc[k] += fabs(((pActiveLFN)->DATA.LFN.adblW)[k + 1]);
			}
		}
		average = sum / (double)nalns; // this is the result of averaging [0,1]-limited ALNs???
		desired = TVfile.GetAt(j, nDim - 1, 0); // get the desired result	
		se += (average - desired) * (average - desired);
		if (fabs(desired - average) > 0.5)  nClassError++; // desired must be integer
	}
	double rmse = sqrt(se / ((double)nRowsTV - 1.0)); // frees se for use below.
	// get the average weight on all variables k
	for (k = 0; k < nDim; k++)
	{
		adblWAcc[k] /= (nRowsTV * nalns);
		adblAbsWAcc[k] /= (nRowsTV * nalns);
	}
	fprintf(fpProtocol, "Size of datasets PP TV Test %d  %d  %d \n", nRowsPP, nRowsTV, nRowsTS);
	fprintf(fpProtocol, "Root mean square error of the average over %d ALNS is %f \n", nalns, rmse);
	fprintf(fpProtocol, "Warning: the above result is optimistic, see results on the test set below\n");
	fprintf(fpProtocol, "Importance of each input variable:\n");
	fprintf(fpProtocol, "Abs imp = stdev(input var) * average absolute weight / stdev(output var) \n");
	fprintf(fpProtocol, "Abs imp is numerical and indicates ups and downs in output when the given input varies.\n");
	fprintf(fpProtocol, "For example a sawtooth function with six teeth would have importance 12.\n");
	fprintf(fpProtocol, "First we have to compute the standard deviation of the output variable.\n");
	fflush(fpProtocol);
	//compute the average of the output variable in the TVset
	k = nDim - 1;
	desired = 0;
	for (j = 0; j < nRowsTV; j++)
	{
		desired += TVfile.GetAt(j, k, 0);
	}
	desired /= nRowsTV; // now desired holds the average for variable k

	// compute the standard deviation of the output variable in the TVset
	se = 0;
	double temp;
	for (j = 0; j < nRowsTV; j++)
	{
		temp = TVfile.GetAt(j, k, 0);
		se += (temp - desired) * (temp - desired);
	}
	se /= ((double)nRowsTV - 1.0); // sample variance of the output variable
	double stdevOutput = sqrt(se);
	fprintf(fpProtocol, "\nStandard deviation of the output in the TVfile %f\n", stdevOutput);
	if (fabs(stdevOutput) < 1e-10)
	{
		fprintf(fpProtocol, "\nStopping: The standard deviation of the output on the TV set is near 0.\n");
		fprintf(fpProtocol, "Average output value is %f to great accuracy. \nEnd of this file.\n", desired);
		fclose(fpProtocol);
		exit(0);
	}
	// we compute the variance of each column of TV
	se = 0;
	for (k = 0; k < nDim - 1; k++) // do each variable k
	{
		//compute the average of variable k in TVset
		desired = 0;
		for (j = 0; j < nRowsTV; j++)
		{
			desired += TVfile.GetAt(j, k, 0);
		}
		desired /= nRowsTV; // now desired holds the average for variable k

		// compute the standard deviation of variable k in TVset
		se = 0;
		double temp;
		for (j = 0; j < nRowsTV; j++)
		{
			temp = TVfile.GetAt(j, k, 0);
			se += (temp - desired) * (temp - desired);
		}
		se /= ((double)nRowsTV - 1.0); // sample variance of variable k
		dblImportance[k] = sqrt(se) * adblAbsWAcc[k] / stdevOutput;
		if (nLag[k] == 0)
		{
			fprintf(fpProtocol, "Variable %s: stdev = \t%f; avg.wt = \t%f; abs imp = \t%f\n",
				varname[nInputCol[k]], sqrt(se), adblWAcc[k], dblImportance[k]);
		}
		else
		{
			fprintf(fpProtocol, "Variable %s@lag%d: stdev = \t%f; avg.wt = \t%f; abs imp = \t%f\n",
				varname[nInputCol[k]], nLag[k], sqrt(se), adblWAcc[k], dblImportance[k]);
		}
		// we use the product of the variance of k and the average absolute weight as a measure of importance
	}
	if (bClassify)
	{
		fprintf(fpProtocol, "Number of TV file cases misclassified = %d out of %d\n", nClassError, nRowsTV);
		fprintf(fpProtocol, "Percentage of TV file cases misclassified = %f", 100.0*(double)nClassError / (double)nRowsTV);
	}
	fflush(fpProtocol);
	free(adblX);
	free(adblWAcc);
	free(adblAbsWAcc);
}

void ALNAPI trainAverage() // routine
{
  int nalns;
  nalns = nALNs;
	if (nalns == 1) return; // In this case we skip bagging (averaging several ALNs)
	bTrainingAverage = TRUE;
	fprintf(fpProtocol,"\n**** Training an ALN by resampling the average of approximations ******\n");
	dblLimit = dblLinRegErr * 1e-5;// Making this larger than 0 may reduce bagging time.
	if (dblLimit > 0)
	{
		fprintf(fpProtocol, "The value used for stopping splitting of the average ALN is %f \n", dblLimit);
	}
	else
	{
		fprintf(fpProtocol, "An F-test is used for stopping splitting of the average ALN\n");
	}
	// The noise variance values in NVfile are the previous ones divided by nalns.
	createTR_NVfiles(BAGGING);
	pAvgALN = new CMyAln; // NULL initialized ALN
	if (!(pAvgALN->Create(nDim, nDim-1) &&
				pAvgALN->SetGrowable(pAvgALN->GetTree())))
	{
    fprintf(fpProtocol,"Stopping: Growable average ALN creation failed!\n");
    exit(0);
	}
	bALNgrowable = TRUE;
	// Get epsilons, centroids and weights in neighborhood of good values to start training
	prepareQuickStart(pAvgALN);
	// We can get rid of the quick start data
	delete[] adblLRC; // delete these from free store
	delete[] adblLRW;
	adblLRC = NULL; // set the pointers to NULL
	adblLRW = NULL;
	fprintf(fpProtocol,"Smoothing epsilon same as for approximation\n\n");
	fflush(fpProtocol);
	nRowsTR = TRfile.RowCount();
	const double* adblData = TRfile.GetDataPtr();
	pAvgALN->SetDataInfo(nRowsTR, nDim, adblData, NULL);
	nMaxEpochs = 20;
	nNumberLFNs = 1;  // initialize at 1
	dblMinRMSE = dblLinRegErr * 1e-6; // Stops training when the error is tiny 
	// Stopping splitting uses the F-test with a noise variance in NVfile
	// that has been divided by nALNs.
	dblLearnRate = 0.3;
	int iterations = 30; // This may have to be adjusted.
	if(bEstimateNoiseVariance)
	{
		for (int iteration = 0; iteration < iterations; iteration++)
		{
			fprintf(fpProtocol, " %d ", iteration);
			// TRAIN AVERAGE ALN vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
			if (!pAvgALN->Train(nMaxEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // fast change to start, dblLearnRate 1.0 for 2 epochs
			{
				fprintf(fpProtocol, "Average ALN training failed!\n");
			}
			fflush(fpProtocol);
			if (bStopTraining == TRUE)
			{
				fprintf(fpProtocol, "This training stopped because all leaf nodes have stopped changing!\n");
				bStopTraining = FALSE;
				fprintf(fpProtocol, "\nTraining of the average ALN completed at iteration %d \n", iteration);
				fflush(fpProtocol);
				break;
			}
			if(iteration == 5)dblLearnRate = 0.2;
			if(iteration == 15)dblLearnRate = 0.15;
			if(iteration == 25)dblLearnRate = 0.1;
			fprintf(fpProtocol, "Learning rate = %f\n", dblLearnRate);
		}
	}
	fflush(fpProtocol);
}

void ALNAPI constructDTREE(int nMaxDepth) // routine
{
	// ******************  CONSTRUCT A DTREE FOR THE AVERAGE ALN *******************************
  fprintf(fpProtocol,"\n***** Constructing an ALN decision tree from the average ALN *****\n");
	DTREE* pAvgDTR;

	// Create a single-layer DTREE directly with ConvertDtree.
	// Setting nMaxDepth to a higher value than 1
	// requires a lot of time to generate the DTREE, but
	// allows much faster evaluation. This could turn out to be
	// useful for extremely demanding real-time tasks like
	// controlling nuclear fusion in ITER.
	if (nALNs > 1)
	{
		pAvgDTR = pAvgALN->ConvertDtree(nMaxDepth);
	}
	else
	{
		pAvgDTR = apALN[0]->ConvertDtree(nMaxDepth);
	}
  if(pAvgDTR == NULL)
	{
		fprintf(fpProtocol,"No DTREE was generated from the average ALN. Stopping. \n");
    exit(0);
	}
	else
	{
		WriteDtree(szDTREEFileName,pAvgDTR);
		fprintf(fpProtocol,"The DTREE of the average of ALNs %s  was written.\n",szDTREEFileName );
	  fflush(fpProtocol);
  }
}

void ALNAPI cleanup() // routine
{
  if(bTrain)
  {
		// cleanup what was allocated for training in approximation
		for (int n = 0; n < nALNs; n++)
		{
			apALN[n]->Destroy();
		}
		free(apALN);
		if(nALNs > 1 ) pAvgALN->Destroy();
    // the TV file is not created for evaluation
		TVfile.Destroy();
		TSfile.Destroy();
    TRfile.Destroy();
		NVfile.Destroy();
    free(adblEpsilon);
	}
  
  // clean up the rest
	PreprocessedDataFile.Destroy();                                       
	free(adblMinVar);
	free(adblMaxVar);
	free(adblStdevVar);
  fclose(fpOutput);
  fclose(fpProtocol);
}

// file task.cpp

void fillvector(double * adblX, CMyAln* paln) // routine
// This is called by the callback in CMyAln to fill in a data vector for training.
// If adblData is set by TRfile.GetDataPtr() and used in paln.SetDataInfo(nPoints,nCols,adblData)
// then the callback does not need to call this routine to get a training vector.
// This routine used to be used to do bagging, but now it uses FillInputVector(...)
// instead, a routine that now trains on an average of ALNs in TRfile.
{
	/* IMPORTANT: This routine can be adapted to create training vectors
	// for real-time applications.
	long nRow;
	nRow = (long)floor(ALNRandFloat() * (double) nRowsTR); // This is where the TRfile is indicated for training.
	for(int i = 0; i < nDim; i++)
	{
		adblX[i] = TRfile.GetAt(nRow,i,0); // Notice that TRfile is fixed.
		// Global nRowsTR is not, so that can be varied.
		// To get data in real time, you have to write new code.
	}
	if(bTrainingAverage) // In this case, we use only the domain components from TRfile.
	{
		// To average several ALNs, we create new samples around those in TRfile and average the ALN outputs
		// The average ALN is created with negligible smoothing
		const ALNCONSTRAINT* pConstr;
		int nalns = nALNs;
		double dblValue;
		ALNNODE* pActiveLFN;
		for(int i = 0; i < nDim -1; i++)
		{
			// The idea here is the same as jitter, and gives better sampling.
			// The chosen point is triangularly distributed around the initial point
			pConstr = paln->GetConstraint(i, 0);
			ASSERT(pConstr != NULL);
			adblX[i] += (ALNRandFloat() - ALNRandFloat()) * pConstr->dblEpsilon;
			// this dblEpsilon is the size of a box "belonging to" a point in the i-th axis
		}
		double sum = 0;
		for (int n = 0; n < nalns; n++)
		{
			dblValue = apALN[n]->QuickEval(adblX, &pActiveLFN);
			sum += dblValue;
		}
		// put the bagging result into the training vector
		adblX[nDim-1] = sum / (double) nalns; // this is the result of averaging ALNs for vector j
	}
	*/
}


void ALNAPI createTR_NVfiles(int nChoose) // routine
{
	// This routine uses the TVfile to set up TRfile and NVfile in various ways.
	// The TVfile is all of the PreprocessedDataFile which is not used for testing in TSfile.
	// 
	// nChoose = 0: LINEAR_REGRESSION. The TVfile  is copied, in a different order
	// (destroying a  possibly unsuitable order) into TRfile and NVfile.
	// TRfile is used for linear regression and creating noise variance samples.
	// Ready for the approximation step, the NVfile will contain smoothed noise variance samples.
	// 
	// nChoose = 1 NOISE_VARIANCE. The NV file is initially as was copied from TVfile,
	// and serves to generate noise variance samples in the TV file. Those samples
	// are not at the same places in the domain. In this routine, NV values are squared and
	// multiplied by a constant to get the "noise variance of the noise variance"
	// and training of a noise variance ALN is carried out.
	// 
	// nChoose = 2 APPROXIMATION. The TRfile is reloaded from TVfile
	// nChoose = 3: BAGGING. For averaging several ALNs, the values of the noise
	// variance samples are divided by the number of ALNs averaged.
	// Again, this avoids overtraining.
	fprintf(fpProtocol, "TVfile gives rise to TRfile for training and NVfile for noise variance.\n");
	double dblValue;
	long i;
	int j;
	if (nChoose == 0) // LINEAR_REGRESSION
	{
		// Create the files
		fprintf(fpProtocol, "Setting up TRfile and NVfile for linear regression\n");
		nRowsTR = nRowsNV = TVfile.RowCount();
		TRfile.Create(nRowsTR, nDim);
		NVfile.Create(nRowsNV, nDim);
		// First we fill TRfile and NVfile from TVfile
		for (i = 0; i < nRowsTV; i++)
		{
			for (j = 0; j < nDim; j++) // ... to the front or
			{
				dblValue = TVfile.GetAt(i, j, 0);
				TRfile.SetAt(i, j, dblValue, 0);
				NVfile.SetAt(i, j, dblValue, 0);
			}
		} // end of i loop
		//if (bPrint && bDiagnostics) NVfile.Write("DiagnoseNVfileLR.txt");
		//fflush(fpProtocol);
	}	// end if(nChoose == 0) LINEAR_REGRESSION

	if (nChoose == 1) // NOISE_VARIANCE
	{
		// At the start of the noise variance routine, both TRfile and NVfile are
		// the same as they were during linear regression, i. e. the same as TVfile.
		// Then NVfile is used to create noise variance samples in TRfile. Those
		// samples also have their domain location slightly changed. They
		// are to be used for training a noise variance ALN NV_ALN. But before that
		// training can occur, the "noise of the noise" has to be put into NVfile
		// by the present code.
		double dblFactor = 1.0 / (Nearby + 1.0 - nDim);
		fprintf(fpProtocol,
			"Setting up TRfile and NVfile for training the noise variance ALN\n");
		for (i = 0; i < nRowsTV; i++)
		{
			// first we change the domain point to the midpoint used for i
			for (j = 0; j < nDim - 1; j++)
			{
				dblValue = TRfile.GetAt(i, j, 0);
				NVfile.SetAt(i, j, dblValue, 0);
			}
			// then we put in a correct value for the noise of the noise
			dblValue = TRfile.GetAt(i, nDim - 1, 0);
			NVfile.SetAt(i, nDim - 1, dblFactor * dblValue * dblValue, 0);
		} // end of i loop
	}// end if(nChoose == 1) NOISE_VARIANCE

	if (nChoose == 2) // APPROXIMATION 
	{
		fprintf(fpProtocol, "Setting up TRfile and NVfile for approximation\n");
		// We must use TRfile for training because of the way we did splitops.
		// This should be corrected so we can use any file for training.
		for (i = 0; i < nRowsTV; i++)
		{
			for (j = 0; j < nDim; j++)
			{
				dblValue = TVfile.GetAt(i, j, 0);
				TRfile.SetAt(i, j, dblValue, 0);
				// NVfile was set up using createSamples
			}
		}
	}

	if (nChoose == 3) // BAGGING
	{
		fprintf(fpProtocol, "Setting up TRfile and NVfile for bagging\n");
		int nalns = nALNs;
		double reciprocalalns = 1.0 / (double)nalns;
		double reciprocalalnsm1 = 1.0 / ((double)nalns - 1.0);
		double dblMaxVar = 0;
		// Now we generate a new training file based on averaging the
		// values of the ALNs generated during approximation.
		// We first try this just averaging the values at existing
		// samples in the TRfile. (The old way used fillvector and jittered inputs.)
		// While we are doing this, we determine the variance among the averaged ALN
		// values and put a constant times that in NVfile at the same location.
		ALNNODE* pActiveLFN = NULL;
		double * adblX = (double *)malloc((nDim) * sizeof(double));
		for (i = 0; i < nRowsTR; i++)
		{
			for (int k = 0; k < nDim; k++)
			{
				adblX[k] = TRfile.GetAt(i, k, 0); // the nDim - 1 component is not used
			}
			double sum = 0;
			double sumsq = 0;
			for (int n = 0; n < nalns; n++)
			{
				dblValue = apALN[n]->QuickEval(adblX, &pActiveLFN);
				sum += dblValue;
				sumsq += dblValue * dblValue;
			}
			double var;
			var = (sumsq + reciprocalalns * sum * sum) * reciprocalalnsm1;
			if (var > dblMaxVar) dblMaxVar = var;
			// put the average into the training file
			TRfile.SetAt(i, nDim - 1, sum * reciprocalalns, 0);
			// At this point we just interpolate
			NVfile.SetAt(i, nDim - 1, 0.0, 0);
		}
		free(adblX);
		fprintf(fpProtocol, "Maximum variance at sample points when averaging is %f\n",
			dblMaxVar);
		fprintf(fpProtocol, "The values in NVfile are set to 0 for bagging, i.e. we are interpolating.");
	} // This ends if (nChoose == 3) BAGGING
} // This ends createTR_NVfiles

void prepareQuickStart(CMyAln* pALN)
{
	// We use information from Linear Regression (LR) to save a bit of training time.
	// Set constraints on variables for the ALN
	// NB The following loop excludes the output variable of the ALN, index nDim -1.
	for (int m = 0; m < nDim - 1; m++) 
	{
		pALN->SetEpsilon(adblEpsilon[m], m);
		if (adblEpsilon[m] == 0)
		{
			fprintf(fpProtocol, "Stopping: Variable %d appears to be constant. Try removing it.\n", m);
			fflush(fpProtocol);
			exit(0);
		}
		// The minimum value of the domain is a bit smaller than the min of the data points
		// in TVfile and the maximum is a bit larger.
		pALN->SetMin(adblMinVar[m] - 0.1 * adblStdevVar[m], m);
		pALN->SetMax(adblMaxVar[m] + 0.1 * adblStdevVar[m], m);
	
		// A rough value for the range of output (for a uniform dist.) divided by the
		// likely distance between samples in axis m.
		pALN->SetWeightMin(-pow(3.0, 0.5) * adblStdevVar[nDim - 1] / adblEpsilon[m], m);
		pALN->SetWeightMax(pow(3.0, 0.5) * adblStdevVar[nDim - 1] / adblEpsilon[m], m);
	
		// Impose the a priori bounds on weights which have been set by the user
		if (dblMinWeight[m] > pALN->GetWeightMin(m))
		{
			pALN->SetWeightMin(dblMinWeight[m], m);
		}
		if (dblMaxWeight[m] < pALN->GetWeightMax(m))
		{
			pALN->SetWeightMax(dblMaxWeight[m], m);
		}
	}
	if(bALNgrowable && !bTrainNV_ALN) // This setup is not used for LR or training NV_ALN
	{
		(pALN->GetRegion(0))->dblSmoothEpsilon = 0;
		// use the weights and centroids from linear regression
		ALNNODE* pActiveLFN;
		pActiveLFN = pALN->GetTree();
		for (int m = 0; m < nDim - 1; m++)
		{
			((pActiveLFN)->DATA.LFN.adblC)[m] = adblLRC[m];
			((pActiveLFN)->DATA.LFN.adblW)[m + 1] = adblLRW[m + 1];
		}
		((pActiveLFN)->DATA.LFN.adblC)[nDim - 1] = adblLRC[nDim - 1];
		((pActiveLFN)->DATA.LFN.adblW)[0] = adblLRW[0];
		((pActiveLFN)->DATA.LFN.adblW)[nDim] = -1.0;
	}
}

double dist(double* adblA, double* adblB)
{
	double sum = 0.0;
	for (int j = 0; j < nDim - 1; j++) // the sample value at nDim - 1 doesn't matter
	{
		sum += pow(adblA[j] - adblB[j], 2);
	}
	return pow(sum, 0.5);
}

double scalarProduct(double* adblA, double* adblB)
{
	double sum = 0.0;
	for (int j = 0; j < nDim - 1; j++)
	{
		sum += adblA[j] * adblB[j];
	}
	return sum;
}

void createSamples(CMyAln* pNV_ALN)  // routine
{
	// The goal is to smooth the samples in NVfile for approximation
	// using the result of training of log10 NV_ALN in TRfile.
	ASSERT(pNV_ALN);
	ALNNODE* pActiveLFN;
	double dblALNValue;
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	long nRowsNV = NVfile.RowCount();
	double nRowsNVrecip = 1.0 / (double) nRowsNV;
	double ReductionFactor = 1.0;// A value slightly less than 1.0
	// will cause splitting to continue longer, somewhat compensating for
	// curvature and slope of the function to be learned, which make
	// NV larger than it should be. This might result in a better fit.
	fprintf(fpProtocol,"The reduction factor of noise variance for training is %f\n",
		ReductionFactor);
	for (long i = 0; i < nRowsNV; i++)
	{
		for (int j = 0; j < nDim; j++)
		{
			adblX[j] = TRfile.GetAt(i, j, 0);
			NVfile.SetAt(i, j, adblX[j], 0);
		}
		dblALNValue = pNV_ALN->QuickEval(adblX, &pActiveLFN);
		NVfile.SetAt(i, nDim - 1, max(ReductionFactor * dblALNValue, noiseSampleMin), 0); // MYTEST put factor 0.5 to compare to below
		// You can use ReductionFactor = 3/4 on example Noisysin20000.txt as follows.
		// NVfile.SetAt(i, nDim - 1, pow(1.00461579,2.0 * adblX[0]), 0);
		// The noise added is 2 * (1.00461579^A1) * 2 * (RAND() - 0.5)
		// with variance 4* 1.00461579^(2*A1) * 4 * 0.5^2 /3.
	}
	free(adblX);
	pNV_ALN->Destroy();
	// Now check to see the global noise variance 
	// Check a case of known constant noise variance!
	double arithAvg = 0;
	for(long i = 0; i < nRowsNV; i++)
	{
		arithAvg += NVfile.GetAt(i, nDim - 1, 0);
	}
	dblSetTolerance = arithAvg / nRowsNV;
	fprintf(fpProtocol, "Arithmetic average of noise variance samples = %f\n", dblSetTolerance);
	fflush(fpProtocol);
	if(bPrint && bDiagnostics) NVfile.Write("DiagnoseNVfileAfterNV_ALN_train.txt");
	if(bPrint && bDiagnostics) fprintf(fpProtocol, "Diagnose NVfileAfterNV_ALN_train.txt written\n");
}
