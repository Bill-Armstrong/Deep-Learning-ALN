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
void ALNAPI createNoiseVarianceFile();
double dist(double*, double*); // calculates the distance between domain points
#define PrintInterval 25
//defines used to set up TRfile and VARfile
#define LINEAR_REGRESSION 0
#define NOISE_VARIANCE 1
#define APPROXIMATION 2
#define BAGGING 3

// We use dblRespTotal in two ways and the following definition helps.
#define DBLNOISEVARIANCE dblRespTotal

// files used in training operations
CDataFile TRfile; // Training data, changed for different purposes.
CDataFile VARfile; // Becomes a file of noise variance samples.

//routines
void ALNAPI doLinearRegression(); // Determines an upper bound on error, and provides a start for other training.
void ALNAPI createNoiseVarianceFile(); // This creates the noise variance file VARfile.
void ALNAPI approximate(); // Actually does training avoiding overtraining using samples in VARfile.
void ALNAPI outputTrainingResults();
void ALNAPI trainAverage(); // Takes several ALNs created in approximate() and creates an ALN of their average
void ALNAPI constructDTREE(int nMaxDepth); // Takes the average ALN and turns it into a DTREE
void ALNAPI cleanup(); // Destroys ALNs etc.
void fillvector(double * adblX, CMyAln* paln); // Sends a vector to training from a file, online or averaging.
void ALNAPI createTR_VARfiles(int nChoose);
//void ALNAPI trainNoiseVarianceALN(); // Trains on the samples in VARfile to create NV_ALN
//void createSamples(CMyAln* pNV_ALN); // Creates smoothed noise variance samples using NV_ALN
void prepareQuickStart(CMyAln* pALN); // This sets up ALN training with values already determined by Linear Regression

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
long nRowsVAR; // the number of rows in the noise variance file.  When approximation starts, this should be nRowsTV
double dblLimit = -1  ;// A negative value splits pieces based on an F test, otherwise they split if training MSE < dblLimit.
int nEpochSize; // the number of samples in the current training set; one epoch = one pass through the set
BOOL bALNgrowable = TRUE; // FALSE for linear regression, TRUE when the ALN has to grow.
//BOOL bTrainNV_ALN = FALSE; // Controls setup of training for to fit noise samples in VARfile
BOOL bStopTraining = FALSE; // Set to TRUE and becomes FALSE if any (active) linear piece needs training
int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO; // used with callbacks
double * adblX = NULL;

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
	//bTrainNV_ALN = FALSE; // TRUE only during training of the noise variance ALN.
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
	createTR_VARfiles(LINEAR_REGRESSION);
	int nRowsTR = TRfile.RowCount();	// nEpochsize gives the number of training samples. Later nRowsVAR=nRowsTR.
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
	int nIterations = 30;
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
	double value = 0;
	double* aX = NULL;
	double* aY = NULL;
	double noiseSampleSum = 0;
	long inside = 0; // counts samples inside their symplex of closest nDim points.
	double yAtBarycentre;  // using the value at the Barycentre to subtract from the sample value
	double Correction = ((float)nDim) / ((float)nDim + 1); // a correction factor to adjust variance (Barycentre!)
	struct disc // each sample (X,y) has a disc struct at the same index
	{
		long*   aXX; // These are indices of samples nearby
		double* aDD; // These are the distances of sample j to (X,y)
		int maxloc;      // Location in aXX of a sample at the maximum distance from the central one
	};
	disc* adisc = NULL;
	// Example: domain 8 dimensional, samples 9 components, nDim = 9, we seek any 9 closest samples in the domain.
	// ...TO DO Preprocessing: N samples with equal X are grouped and have the average y and N attached.
	// Processing:
	// Each sample in turn has the role of central sample X. It has a struct with an array of nDim sample indices
	// of other samples. The other samples have distances from the central sample recorded and their
	// maximum distance from the central sample is tracked. Initially, the distances are all extremely large
	// as is the maximum distance.  The initial indices don't matter because they will all be replaced.
	// We go through all samples Y to find any set of nDim closest. If we come across the index of the central
	// sample, we go on to the next sample. Otherwise we compute the distance d1 of Y to X.
	// With preprocessing d1 can't be zero; but here we just ignore the sample. If d1 >= maxdist we also
	// ignore it. If d1 < maxdist, the dist array is searched and Y replaces any sample at maxdist.
	// The dist array and maxdist are then adjusted.
	aX = (double*)malloc(nDim * sizeof(double));
	aY = (double*)malloc(nDim * sizeof(double));
	adisc = (disc*)malloc(nRowsTR * sizeof(disc));
	long i, j, k;
	double d1;
	for (i = 0; i < nRowsTR; i++) // i is the central sample its struct will be worked on.
	{
		// First we have to create the struct for i
		adisc[i].aXX = (long*)malloc(nDim * sizeof(long));
		adisc[i].aDD = (double*)malloc(nDim * sizeof(double));
		// Get the struct data associated with the central sample aX with index i.
		for (j = 0; j < nDim; j++)
		{
			aX[j] = TRfile.GetAt(i, j, 0); // Get the central sample's domain components.
			adisc[i].aXX[j] = 0;          // This 0 will be replaced if there are at least nDim + 1 samples.
			adisc[i].aDD[j] = DBL_MAX;    // Set all distances of the i-th struct to the max possible.
		}
		adisc[i].maxloc = nDim - 1;
		// We fill the struct with samples now, we have to go through ALL samples to get
		// the nDim closest samples. This can be made faster than O(n^2), maybe O(n log n).
		for (j = 0; j < nRowsTR; j++)
		{
			if (j != i) // insert only samples j into the list of nearby samples which are not the central sample
			{
				// Get the domain point for the j-th sample into aY
				for (k = 0; k < nDim - 1; k++)
				{
					aY[k] = TRfile.GetAt(j, k, 0);
				}
				d1 = dist(aX, aY);
				int imaxloc = adisc[i].maxloc;
				if (d1 < adisc[i].aDD[imaxloc])
				{
					// insert the point closer to X at the place indicated by maxloc
					// we replace the index of a sample at the maximum distance with the index of aY
					adisc[i].aXX[imaxloc] = j;
					adisc[i].aDD[imaxloc] = d1; // the distance at m becomes d1
					// find the new maximum distance and a sample index and where it occurs
					double current = 0;
					for (int check = 0; check < nDim; check++)
					{
						if (adisc[i].aDD[check] > current)
						{
							current = adisc[i].aDD[check];
							adisc[i].maxloc = check;
						}
					}
				} // drop the samples with d1 > maximum distance
			} // drop samples with j == i
		} // end of j loop
	} // loop over i
	// We have nDim samples closest to each sample index i
	// Now we analyze geometrically the collections of samples in the structs
	// In the following, M and z have additional 1 entries at the bottom row (for barycentric coordinates)
	MatrixXd M(nDim, nDim); // the nDim nearby domain points to the i-th, arranged in columns
	VectorXd z(nDim); // the central sample domain point
	VectorXd y(nDim); // the nDim nearby sample values
	VectorXd alpha(nDim); // barycentric coordinates of central sample w.r.t. nDim points nearby
	long index;
	double q;
	for (i = 0; i < nRowsTR; i++)
	{
		// get the central sample's domain vector as a column with a 1 at the bottom
		for (k = 0; k < nDim - 1; k++)
		{
			z(k) = TRfile.GetAt(i, k, 0); // get the central (i-th) sample's domain point
		}
		z(nDim - 1) = 1.0; // add a 1 at the bottom
		q = TRfile.GetAt(i, nDim - 1, 0); // the value of the central sample
		// get the nearby vectors as columns, synthesize sample value at barycentre
		yAtBarycentre = 0;
		for (j = 0; j < nDim; j++) // j column index over nDim samples
		{
			index = adisc[i].aXX[j];
			for (k = 0; k < nDim - 1; k++) // k row index over domain components
			{
				M(k, j) = TRfile.GetAt(index, k, 0); // the j-th vector as a column in M
			}
			M(nDim - 1, j) = 1.0; // 1 at the bottom
			y(j) = TRfile.GetAt(index, nDim - 1, 0); // the sample values as an nDim x 1 column
			yAtBarycentre += y(j);
		}
		yAtBarycentre /= nDim; // subtracted from the value of sample i to get a noise estimate
		// Find the barycentric coordinates of the central sample -- useful to see if
		// the simplex defined by the closest nDim points to X actually contains X. 
		// alpha = M.colPivHouseholderQr().solve(z);
		// Express the Noise Variance related to dB using the ratio of aplitudes
		double NV = pow(q - yAtBarycentre, 2)* Correction * 0.75; // the 0.75 is a tweak to compensate for slope
		noiseSampleSum += NV;
		VARfile.SetAt(i, nDim - 1, NV, 0); // the domain components should be unchanged
	} // end i loop over nRowsTR
	// Report on VARfile
	fprintf(fpProtocol, "Average Noise Variance before smoothing %f \n ",
		noiseSampleSum / nRowsVAR);
	VARfile.Write("DiagnoseVARfileBeforeSmoothing.txt");
	// Now we smooth the VARfile by averaging over the discs.
	double smooth;
	double Dimp1 = (double)nDim + 1;
	int nIter = 20;
	noiseSampleSum = 0;
	for (int iteration = 0; iteration < nIter; iteration++)
	{
		for (i = 0; i < nRowsVAR; i++)
		{
			smooth = 0;
			for (j = 0; j < nDim; j++) // j column index over nDim samples
			{
				index = adisc[i].aXX[j];
				smooth += VARfile.GetAt(index, nDim - 1, 0);
			}
			smooth += VARfile.GetAt(i, nDim - 1, 0);
			smooth /= Dimp1;
			VARfile.SetAt(i, nDim - 1, smooth, 0);
			if (iteration == nIter - 1) noiseSampleSum += smooth;
		}
	}
	fprintf(fpProtocol, "Average Noise Variance after smoothing %f ",
		noiseSampleSum / nRowsVAR);
	VARfile.Write("DiagnoseVARfileAfterSmoothing.txt");
	/*
	// TEST  put in the real value of noise variance!
	// For this to be correct, there must be no test set.
	for (i = 0; i < nRowsVAR; i++)
	{
		ASSERT(fabs(VARfile.GetAt(i, 0, 0) - i) < 0.0001);
		ASSERT(fabs(TRfile.GetAt(i, 0, 0) - i) < 0.0001);
		VARfile.SetAt(i, 1, 0.01333333 * (i + 1) * (i + 1), 0);
	}
	*/
	free(aX);
	free(aY);
	free(adisc);
}

/*
void ALNAPI trainNoiseVarianceALN()
{
	fprintf(fpProtocol, "\n**************Training Noise Variance ALN begins ********\n");
	fflush(fpProtocol);
	if (bJitter)
	{
		fprintf(fpProtocol, "Jitter is used during approximation\n");
	}
	else
	{
		fprintf(fpProtocol, "Jitter is not used during approximation\n");
	}
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
	// Now make the tree growable
	if (!pNV_ALN->SetGrowable(pNV_ALN->GetTree()))
	{
		fprintf(fpProtocol, "Setting NV_ALN growable failed!\n");
		fflush(fpProtocol);
		exit(0);
	}
	bALNgrowable = TRUE;
	bTrainingAverage = FALSE;
	bTrainNV_ALN = TRUE; // This causes the bounds on weights to be tighter than for other
		//training (see prepareQuickStart routine)
	(pNV_ALN->GetRegion(0))->dblSmoothEpsilon = 0.0;

	nMaxEpochs = 20;
	dblMinRMSE = 0.0;
	dblLearnRate = 0.2;
	fprintf(fpProtocol, "The smoothing for training NV_ALN is %f, and the learning rate is %f\n",
		(pNV_ALN->GetRegion(0))->dblSmoothEpsilon, dblLearnRate);
	bStopTraining = FALSE; // Set TRUE in alntrain.cpp. Set FALSE by any piece needing more training. 
	nNumberLFNs = 1;  // initialize at 1
	nRowsVAR = VARfile.RowCount();
	ASSERT(nRowsVAR == nRowsTV);
	createTR_VARfiles(NOISE_VARIANCE); // We have to fix things up so we don't have to always train on TRfile because of splitops
	const double* adblData = TRfile.GetDataPtr();
	pNV_ALN->SetDataInfo(nRowsTR, nDim, adblData, NULL); // Not possible yet to train on VARfile
	dblLimit = 1.0;  // TEST !!!  How should we do this?????
	fprintf(fpProtocol, "----------  Training NV_ALN  ------------------\n");
	fflush(fpProtocol);
	for (int iteration = 0; iteration < 20; iteration++) // is 10 iterations enough?
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
	} // end of loop of training interations for NV_ALN
	free(adblX);
	bTrainNV_ALN = FALSE;
	createSamples(pNV_ALN); // This changes the values in VARfile according to NV_ALN
	pNV_ALN->Destroy();
}
*/

void ALNAPI approximate() // routine
{
	fprintf(fpProtocol, "\n**************Approximation with one or more ALNs begins ********\n");
	fflush(fpProtocol);
	int nalns = nALNs;  // The number of ALNs over which we average (for "bagging")
	createTR_VARfiles(APPROXIMATION);  // prepares for using the whole TVfile and the whole VARfile with noise variance samples for training and stopping
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

	dblLimit =-1.0; // Negative to split according to an F test, otherwise split if training MSE > dblLimit.
	if (dblLimit <= 0)
	{
		fprintf(fpProtocol, "An F test is used to decide whether to split a piece depending on hit count. \n");
	}
	else
	{
		fprintf(fpProtocol, "A manually set limit, %f, is used to decide whether to split a piece. \n", dblLimit);
	}
	fflush(fpProtocol);
	// ***************** SET UP THE ARRAY OF POINTERS TO ALNS FOR TRAINING ONLY *************************
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
		//bTrainNV_ALN = FALSE;
		if(bClassify)
    {
     // TO DO
    }
		(apALN[n]->GetRegion(0))->dblSmoothEpsilon = 0.0;
		fprintf(fpProtocol, "The smoothing for training each approximation is %f\n", 0.0); 
		nMaxEpochs = 20;
		dblMinRMSE = 0.0;
		dblLearnRate = 0.2;
		bStopTraining = FALSE; // Set TRUE in alntrain.cpp. Set FALSE by any piece needing more training. 
    nNumberLFNs = 1;  // initialize at 1
		// Set up the data
		// Tell the training algorithm the way to access the data using fillvector
		nRowsTR = TRfile.RowCount();
		ASSERT(nRowsTR == nRowsTV);
		const double* adblData = TRfile.GetDataPtr();
		apALN[n]->SetDataInfo(nRowsTR, nDim, adblData, NULL);
		fprintf(fpProtocol,"----------  Training approximation ALN %d ------------------\n",n);
		fflush(fpProtocol);
		for(int iteration = 0; iteration < 15; iteration++) // is 15 iterations enough?
		{
			fprintf(fpProtocol, "\nIteration %d ALN %d ", iteration, n);
			fflush(fpProtocol);

			// TRAIN ALNS WITHOUT OVERTRAINING   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
			if (!apALN[n]->Train(nMaxEpochs, dblMinRMSE, dblLearnRate, bJitter, nNotifyMask))
			{
			  fprintf(fpProtocol,"Training failed!\n");
        exit(0);
			}
			if(bEstimateNoiseVariance)
      {
				fprintf(fpProtocol, " %d ", iteration);
			}
			if (bStopTraining == TRUE)
			{
				fprintf(fpProtocol, "\nTraining of approximation ALN %d completed at iteration %d \n", n, iteration);
				fprintf(fpProtocol, "This training stopped because all leaf nodes have stopped changing!\n");
				bStopTraining = FALSE;
				fflush(fpProtocol);
				break;
			}
			fflush(fpProtocol);
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
		average = sum / (double)nalns; // this is the result of averaging [0,1]-limited ALNs
		desired = TVfile.GetAt(j, nDim - 1, 0); // get the desired result	
		se += (average - desired) * (average - desired);
		if (((desired > 0.5) && (average < 0.5)) || ((desired < 0.5) && (average > 0.5))) nClassError++;
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
	free(adblX);
	free(adblWAcc);
	free(adblAbsWAcc);
	fflush(fpProtocol);
}

void ALNAPI trainAverage() // routine
{
  int nalns;
  nalns = nALNs;
	bTrainingAverage = TRUE;
	// For training the average, we use the TV set to
	// define the region where the data points are located.
	// The values of those points are automatically jittered to cover that region.
	fprintf(fpProtocol,"\n**** Training an ALN by resampling the average of approximations ******\n");
	fprintf(fpProtocol, "The F-limit used for stopping splitting of the average ALN is %f \n", dblLimit);
	// The noise variance values in VARfile are the previous ones divided by nalns.
	createTR_VARfiles(BAGGING);
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
	fprintf(fpProtocol,"Smoothing epsilon same as for approximation\n\n");
	fflush(fpProtocol);
	// Tell the training algorithm about the data, in particular that fillvector must be used (second last NULL)
	nRowsTR = TRfile.RowCount();
	pAvgALN->SetDataInfo(nRowsTR, nDim, NULL, NULL);
	nMaxEpochs = 20;
	dblMinRMSE = 0; // Stopping splitting uses the F-test with a noise variance divided by nALNs.
	dblLearnRate = 0.2;
	int numberIterations = 3;
	if(bEstimateNoiseVariance)
	{
		for (int iteration = 0; iteration < numberIterations; iteration++)
		{
			// TRAIN AVERAGE ALN vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
			if (!pAvgALN->Train(nMaxEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // fast change to start, dblLearnRate 1.0 for 2 epochs
			{
				fprintf(fpProtocol, "Training failed!\n");
			}
			fprintf(fpProtocol, " %d ", iteration);
		}
	}
	else // we have opened a .fit file
	{
		// use the weights and centroids from linear regression to start
		ALNNODE* pActiveLFN;
		pActiveLFN = pAvgALN->GetTree();
		for(int m = 0; m < nDim - 1; m++)
		{
			((pActiveLFN)->DATA.LFN.adblC)[m] = adblLRC[m];
			((pActiveLFN)->DATA.LFN.adblW)[m+1] = adblLRW[m+1];
		}
		((pActiveLFN)->DATA.LFN.adblC)[nDim - 1] = adblLRC[nDim - 1];
		((pActiveLFN)->DATA.LFN.adblW)[0] = adblLRW[0];
		((pActiveLFN)->DATA.LFN.adblW)[nDim] = -1.0;
	}
	//**********	// at this point we have no further use for the starting values from linear regression
	if(bEstimateNoiseVariance)
	{
		delete [] adblLRC; // delete these from free store
		delete [] adblLRW;
		adblLRC = NULL; // set the pointers to NULL
		adblLRW = NULL;
	}
	dblLearnRate = 0.2;
	nEpochSize = nRowsTV; //training average ALN
  nNumberLFNs = 1;  // initialize at 1
  for(int iteration = 0; iteration < 20; iteration++) // is 20 iterations enough?
	{
	  // Call the training function
		bStopTraining = FALSE;
	  if (!pAvgALN->Train(nMaxEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // average
	  {
		  fprintf(fpProtocol,"Average ALN training failed!\n");
      fflush(fpProtocol);
      exit(0);
	  }
		if (bStopTraining == TRUE)
		{
			fprintf(fpProtocol, "This training stopped because all leaf nodes have stopped changing!\n");
			bStopTraining = FALSE;
			fprintf(fpProtocol, "\nTraining of the average ALN completed at iteration %d \n", iteration);
			fflush(fpProtocol);
			break;
		}
		fprintf(fpProtocol,"\nIteration %d of training average ALN, RMSE = %f\n", iteration, dblTrainErr);
		fflush(fpProtocol);
	}
	fflush(fpProtocol);
}

void ALNAPI constructDTREE(int nMaxDepth) // routine
{
	// ******************  CONSTRUCT A DTREE FOR THE AVERAGE ALN *******************************
  fprintf(fpProtocol,"\n***** Constructing an ALN decision tree from the average ALN *****\n");
	DTREE* pAvgDTR;

 // create a single-layer average DTREE directly with ConvertDtree and without splitting
	pAvgDTR = pAvgALN->ConvertDtree(nMaxDepth);
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
		// cleanup just what was allocated for training
		for (int n = 0; n < nALNs; n++)
		{
			apALN[n]->Destroy();
		}
		free(apALN);
		pAvgALN->Destroy();
    // the TV file is not created for evaluation
		TVfile.Destroy();
		TSfile.Destroy();
    TRfile.Destroy();
		VARfile.Destroy();
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
// It uses FillInputVector(...) instead, a routine that can't do the average.
{
	long nRow;
	nRow = (long)floor(ALNRandFloat() * (double) nRowsTR); // This is where the TRfile is indicated for training.
	for(int i = 0; i < nDim; i++)
	{
		adblX[i] = TRfile.GetAt(nRow,i,0); // Notice that TRfile is fixed. Global nRowsTR is not.
		// To get data in real time, you likely have to write new code.
	}
	if(bTrainingAverage) // In this case, we don't use the output component from above. 
	{
		// To average several ALNs, we create new samples around those in TRfile and average the ALN outputs
		// The average ALN is created with negligible smoothing
		const ALNCONSTRAINT* pConstr;
		int nalns = nALNs;
		double dblValue;
		ALNNODE* pActiveLFN;
		for(int i = 0; i < nDim -1; i++) 
		{
			// The idea here is the same as jitter, and gives  much better sampling.
			// The chosen point is triangularly distributed around the initial point
			pConstr = paln->GetConstraint(i, 0);
			ASSERT(pConstr != NULL);
			adblX[i] += (ALNRandFloat() - ALNRandFloat()) * pConstr->dblEpsilon;  // this dblEpsilon is the size of a box "belonging to" a point in the i-th axis
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
}

void ALNAPI createTR_VARfiles(int nChoose) // routine
{
	// This routine uses the TVfile to set up TRfile and VARfile in various ways.
	// The TVfile is all of the PreprocessedDataFile which is not used for testing in TSfile.
	// 
	// nChoose = 0: LINEAR_REGRESSION. The TVfile  is copied, in a different order
	// (destroying a  possibly unsuitable order) into TRfile and VARfile.
	// TRfile is used for linear regression and creating noise variance samples.
	// Ready for the approximation step, the VARfile will contain smoothed noise variance samples.
	// 
	// nChoose = 3: BAGGING. For averaging several ALNs, the values of the noise
	// variance samples are divided by the number of ALNs averaged.
	// Again, this avoids overtraining.
	
	double dblValue;
	long i;
	int j;

	if (nChoose == 0) // LINEAR_REGRESSION
	{
		// Create the files
		nRowsTR = nRowsVAR = TVfile.RowCount();
		TRfile.Create(nRowsTR, nDim);
		VARfile.Create(nRowsVAR, nDim);
		// First we fill TRfile and VARfile from TVfile
		for (i = 0; i < nRowsTV; i++)
		{
			for (j = 0; j < nDim; j++) // ... to the front or
			{
				dblValue = TVfile.GetAt(i, j, 0);
				TRfile.SetAt(i, j, dblValue, 0);
				VARfile.SetAt(i, j, dblValue, 0);
			}
		} // end of i loop
		if (bPrint && bDiagnostics) VARfile.Write("DiagnoseVARfileLR.txt");
		fflush(fpProtocol);
	}	// end if(nChoose == 0) LINEAR_REGRESSION

	if(nChoose == 1) // NOISE_VARIANCE
	{
		// We have to use only TRfile for training because of the way we did splitops.  Should be corrected.
		// For now, the VARfile is copied into the TRfile
		for (i = 0; i < nRowsTV; i++)
		{
			for (j = 0; j < nDim; j++) 
			{
				dblValue = VARfile.GetAt(i, j, 0);
				TRfile.SetAt(i, j, dblValue, 0);
			}
		}
	}

	if (nChoose == 2) // APPROXIMATION 
	{
		// We have to use only TRfile for training because of the way we did splitops.  Should be corrected.
		for (i = 0; i < nRowsTV; i++)
		{
			for (j = 0; j < nDim; j++)
			{
				dblValue = TVfile.GetAt(i, j, 0);
				TRfile.SetAt(i, j, dblValue, 0);
			}
		}
	}

	if (nChoose == 3) // BAGGING
	{
		// Here we again leave the TRfile unchanged but we divide the
		// noise variance values in VARfile by nALNs because of averaging.
		// For averaging, we must use the fillvector routine, which requires setting
		// two parameters to NULL: SetDataInfo(...,..., NULL, NULL)..
		for (i = 0; i < nRowsVAR; i++)
		{
			dblValue = VARfile.GetAt(i, nDim - 1, 0);
			VARfile.SetAt(i, nDim - 1, dblValue/nALNs, 0);
		} // end of i loop
		if (bPrint && bDiagnostics) VARfile.Write("DiagnoseVARfileBAG.txt");
	} // This ends if (nChoose == 3) BAGGING
} // This ends createTR_VARfiles

/*
void createSamples(CMyAln* pNV_ALN)  // routine
{
	// The goal is to smooth the samples in VARfile
	// using the result of training NV_ALN in TRfile
	ASSERT(pNV_ALN);
	ALNNODE* pActiveLFN;
	double dblValue, dblALNValue;
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	long nRowsVAR = VARfile.RowCount();
		for (long i = 0; i < nRowsVAR; i++)
		{
			for (int j = 0; j < nDim; j++)
			{
				adblX[j] = VARfile.GetAt(i, nDim - 1, 0);
			}
			dblValue = adblX[nDim - 1];
			ASSERT(dblValue >= 0);
			dblALNValue = pNV_ALN->QuickEval(adblX, &pActiveLFN);
			if (dblALNValue < 0) dblALNValue = 0;
			VARfile.SetAt(i, nDim - 1, dblALNValue, 0); // keep samples >= 0
		}
	free(adblX);
	// Now check to see the global noise variance (You can comment out what follows if it's proven OK)
	// Check a case of known constant noise variance!
	dblValue = 0;
	for(long i = 0; i < nRowsVAR; i++)
	{
		dblValue += VARfile.GetAt(i, nDim - 1, 0);
	}
	fprintf(fpProtocol, "Average of noise variance samples = %f\n", dblValue / nRowsVAR);
	fflush(fpProtocol);
	if (bPrint && bDiagnostics) VARfile.Write("DiagnoseVARfileAfterNV_ALN.txt");
	if (bPrint && bDiagnostics) fprintf(fpProtocol, "Diagnose VARfileAfterNV_ALN.txt written\n");
}
*/

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
	if(bALNgrowable) // && !bTrainNV_ALN) // This setup is not used for LR or training NV_ALN
	{
		(pALN->GetRegion(0))->dblSmoothEpsilon = adblStdevVar[nDim - 1] / 100.0; // a shot in the dark TEST !!!!
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
