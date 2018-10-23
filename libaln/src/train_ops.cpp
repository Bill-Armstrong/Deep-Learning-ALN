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
//routines
void ALNAPI dolinearregression(); // this determines an upper bound on error
void ALNAPI overtrain(CMyAln * pOT); // this overtrains OTTS and OTVS to obtain the noise variance, which can be used during approximation training later.
void ALNAPI approximate(); // this actually does the training while avoiding overtraining using noise variance samples.
void ALNAPI trainaverage(); // this takes several ALNs created in approximate() and creates an ALN of their average
void ALNAPI constructDTREE(int nMaxDepth); // this takes the average ALN and turns it into a DTREE
void ALNAPI cleanup(); // this destroys ALNs etc.
void fillvector(double * adblX, CMyAln* paln); // this send a vector to training either from a file or, in some applications from online.
// ALN pointers
CMyAln * pALN = NULL; // declares a pointer to an ALN used in linear regression
CMyAln * pOTTS;  // This ALN is overtrained on the original training set (OTTS) to help determine the noise level.
CMyAln * pOTVS; // This ALN is overtrained on the complement (in the TVfile) of the original training set
static CMyAln** apALN = NULL;       // an array of pointers to ALNs used in approximate()
static CMyAln* pAvgALN;      // an ALN representing the bagged average of several ALNs trained on the TVfile with different random numbers

double	dblMinRMSE = 0; // stops training when the training error is smaller than this
double  dblLearnRate = 0.2; // the reciprocal of the learnrate tells how many passes through the data will be required to get close to a fit.
int nNumberEpochs = 10; // if the learnrate is 0.2, then one will need 5 or 10 roughly to almost correct the errors
extern long nRowsTR; // the number of rows in the current training set loaded into TRfile
extern long nRowsVAR; // the number of rows in the noise variance file.  When approximation starts, this should be nRowsTV
extern BOOL bStopTraining; // this boolean is set to FALSE and becomes TRUE if all (active) linear pieces fit well.
double computeGlobalNoiseVariance(); // used at entry to the approximation step, this computes the average of all noise variance samples in VARfile.
int SplitDtree(DTREE** ppDest, DTREE* pSrc, int nMaxDepth);
double dblFlimit =2.59;// this says that splitting of a linear piece is prevented when the mean square
// training error of a piece becomes less than 2.59 times the average of noise variance samples on it.
// This value comes from tables of the F-test for d.o.f. > 7 and probability 90%.
// For 90% with 3 d.o.f the value is 5.39, i.e. with fewer d.o.f. training stops
// when the training error is larger than with a lower F-value.
// We have to IMPROVE the program to use d.o.f. of each piece for both training error and noise variance 
int nEpochSize; // the number of input vectors it takes to correct the errors of position of pieces before splitting

void ALNAPI dolinearregression() // routine
{
  fprintf(fpProtocol,"\n****** Linear regression begins: we want an upper bound on error plus starting weights *****\n");
	// begin linear regression on a single ALN -- this iterative algorithm is not an accepted method in general,
	// but works here. The reason for using it for ALN approximation is that the linear regression
	// problem for a piece is constantly changing. A piece gains or loses input vectors for which it is responsible and
	// for which it computes the output value depending on its value for that input relative to the values of other pieces.
	// This  constantly changes the samples which are to be least-squares fitted.

	// Set up the ALN
	pALN = new CMyAln;
	if (!(pALN->Create(nDim, nDim-1))) // nDim is the number of inputs to the ALN plus one for the ALN output;
																		 // the default output variable is nDim-1
	{
		fprintf(fpProtocol,"Stopping: linear regression ALN creation failed!\n");
    fflush(fpProtocol);
		exit(0);
	}
  fflush(fpProtocol);
	// Set constraints on variables for the ALN	
	for( int m = 0; m < nDim-1; m++) // note that loops that do m < nDim-1 omit the output variable 
		                               // (always the last ALN input, but maybe not connected to the rightmost data file column)
	{
		pALN->SetEpsilon(adblEpsilon[m],m); // adblEpsilon suggests the size along axis m of a box
		// such that the boxes just cover the domain. It can be used to compute jitter, or to
		// estimate the maximum absolute slope of linear pieces that makes sense, given the data.
		// We input all but the output adblEpsilon (index nDim-1) into this ALN. These Epsilons
		// are set in analyzeTV()
		// The values of adblEpsilon show up in the constraints.
    if(adblEpsilon[m] == 0)
    {
      fprintf(fpProtocol,"Stopping: Variable %d appears to be constant. Try removing it.\n",m);
      fflush(fpProtocol);
			exit(0);
    }
		pALN->SetMin(adblMinVar[m] - 0.1 * adblStdevVar[m],m); // we extend the domain a bit beyond
		pALN->SetMax(adblMaxVar[m] + 0.1 * adblStdevVar[m],m); // the range of the data - an evaluation file could have a bigger domain.
		// next we bound weights (partial derivatives) by noting it makes no sense to have a slope that jumps a huge amount between two input data points
		pALN->SetWeightMin(- pow(2.0,0.5) * adblStdevVar[nDim - 1]/adblEpsilon[m],m); // the range of output (for uniform dist.) divided by... 
		pALN->SetWeightMax(  pow(2.0,0.5) * adblStdevVar[nDim - 1]/adblEpsilon[m],m); // ...the likely distance between samples in axis m 
    // impose the a priori bounds on weights
    if(dblMinWeight[m] > pALN->GetWeightMin(m))
    {
		  pALN->SetWeightMin(dblMinWeight[m],m);
    }
    if(dblMaxWeight[m] < pALN->GetWeightMax(m))
    {
		  pALN->SetWeightMax(dblMaxWeight[m],m);
    }
	}
  // createTR_VARfiles(0); done in View; training set TRfile is  about 50% of the TVfile. The rest is for noise variance in VARfile, not used here.
	fprintf(fpProtocol, "TRfile and VARfile created, we start training\n");
	fflush(fpProtocol);
  bTrainingAverage = FALSE; // Switch to tell fillvector whether get a training vector or compute an average
 	(pALN->GetRegion(0))->dblSmoothEpsilon = 0.0; // linear regression: nothing to smooth
	// NB The region concept has not been completed.  It allows the user to impose constraints
	// e.g. on slopes(weights) which differ in different parts of the domain.  All we have is region 0 now.
	int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO;
	//********* dolinearregression
	// Quickstart: get centroids and weights in neighborhood of good values
	dblLearnRate = 0.2;  // roughly, 0.2 corrects about 20% of the error for each pass through the data in nRowsTR.
	int nEpochSize = nRowsTR; // nEpochsize gives the number of training samples. Later nRowsTR=nRowsTV.
	pALN->SetDataInfo(nEpochSize, nDim, NULL, NULL);
	dblMinRMSE = 0; //don't stop early because of low training error
	nNumberEpochs = 100; // we are may have to do tests to see what is sufficient
	if (!pALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask))
	{
		fprintf(fpProtocol,"Initial linear regression training failed!\n");
		fflush(fpProtocol);
		exit(0);
	}
	fprintf(fpProtocol, "Initial linear regression training succeeded!\n");
	fflush(fpProtocol);
	// now we go for a lower learning rate, so the changes become smaller and the piece doesn't jump around so much
	dblLearnRate = 0.05;  // a lower learning rate is called "polishing the ALN"
	// pALN->SetDataInfo(nEpochSize, nDim, NULL, NULL); A lower learning rate usually calls for longer training but we pass here
	dblMinRMSE = 0; // Can be a primitive stopping criterion. Here we don't stop early because of low training error
	nNumberEpochs = 1000; //
	dblFlimit = 0.99;
	if (!pALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask))
	{
		fprintf(fpProtocol, "Polishing linear regression training failed!\n");
		fflush(fpProtocol);
		exit(0);
	}
	fprintf(fpProtocol, "Polishing linear regression training succeeded!\n");
	fflush(fpProtocol);
  // we should now have a good linear regression fit now
  // find the weights on the linear piece using an evaluation at the 0 vector (could be any other place!).
	double * adblX = (double *) malloc((nDim) * sizeof(double));
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
	// now put the negative sum of weighted centroids into the 0 position
	// compress the weighted centroid info into W[0], i.e. take the centroid out of  weight_m*(x_m-centroid_m) + etc.
	// to have just one number in W[0]
	adblLRW[0] = adblLRC[nDim -1] = ((pActiveLFN)->DATA.LFN.adblC)[nDim -1];
  for (int m = 0; m < nDim - 1; m++)
  {
    adblLRW[0] -= adblLRW[m+1] * adblLRC[m];
  }
	adblLRW[nDim] = -1.0;  // the -1 weight on the ALN output may make a difference somewhere, so we set it here
	// The idea of weight = -1 is that the hyperplane equation can be written 
	// a*x + b*y -1*z = 0 like an equation instead of z = a*x + b*y as a function
	// This representation is good if we want to invert the ALN, i.e. get x as a function of y and z.
	// How to do it? Just multiply the first equation by -1/a. The weight on the output is -1.0
	// This is not fully implemented yet.  The output of the ALN is now always the highest index variable.

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

void ALNAPI overtrain(CMyAln * pOT) // this routine overfits the training data in pOT and the complement in pOTVS
{
	// begin fitting a single ALN using new training & variance sets
  // pieces are allowed to break as long as they have over nDim training points on them instead of
	// being constrained to have training error above the noise level as in approximation below
	// We can overtrain using several different splits of the TVfile into two equal parts.
	// This program now uses overtraining on the original training set and its complement.
	fprintf(fpProtocol,"\n*** Overtraining an ALN for use with noise variance estimation begins***\n");
	//  CMyAln * pOT is declared in the header so it can be used elsewhere including in split_ops.cpp
	// Set up the ALN
	pOT = new CMyAln;
	if (!(pOT->Create(nDim, nDim-1)))
	{
		fprintf(fpProtocol,"Stopping: noise estimation overtrained ALN creation failed!\n");
    fflush(fpProtocol);
    exit(0);
	}
  if (!pOT->SetGrowable(pOT->GetTree()))		
	{
	  fprintf(fpProtocol,"Setting overtraining ALN growable failed!\n"); 
    fflush(fpProtocol);
    exit(0);
	}
	// Set constraints on variables for the ALN	
	for( int m = 0; m < nDim-1; m++) // NB Excludes the output variable of the ALN
	{
    if(adblEpsilon[m] == 0) // these are set in analyzeTV(), they are the sides of a box in the domanin being the volume per point
    {
	    fprintf(fpProtocol,"Stopping: Variable %d appears to have 0 variance, i.e. it's constant. Try removing it.\n",m); 
      fflush(fpProtocol);
      exit(0);
    }
		pOT->SetEpsilon(adblEpsilon[m],m); // the Epsilons are set in analyzeTV()
		pOT->SetMin(adblMinVar[m] - 0.1 * adblStdevVar[m],m); // the minimum value of the domain is a bit smaller than the min of the data points in TVfile
		pOT->SetMax(adblMaxVar[m] + 0.1 * adblStdevVar[m],m); // the max is larger by a bit
		pOT->SetWeightMin(- pow(2.0,0.5) * adblStdevVar[nDim - 1]/adblEpsilon[m],m); // the range of output (for uniform dist.) divided by... 
		pOT->SetWeightMax(  pow(2.0,0.5) * adblStdevVar[nDim - 1]/adblEpsilon[m],m); // ...the likely distance between samples in axis m 
    // impose the a priori bounds on weights
    if(dblMinWeight[m] > pOT->GetWeightMin(m))
    {
		  pOT->SetWeightMin(dblMinWeight[m],m);
    }
    if(dblMaxWeight[m] < pOT->GetWeightMax(m))
    {
		  pOT->SetWeightMax(dblMaxWeight[m],m);
    }
	}
	bStopTraining = FALSE; // This becomes TRUE and stops training when all pieces fit well according to the F-test. 
	// We just set the training error here in place of the dblFlimit. The noise variance is set to 1.0 if dblFlimit < 1.0.
	(pOT->GetRegion(0))->dblSmoothEpsilon = 0.0; // for overfitting a single ALN, the smoothing should be zero so as not to interfere.
	int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO;
	fprintf(fpProtocol,"Initial root mean square training error at start of overtraining = %fl  Smoothing = %f \n"
					, dblLinRegErr, (pOT->GetRegion(0))->dblSmoothEpsilon);
	//********* initial weight and centroid components for overtraining
	dblMinRMSE = dblLinRegErr / 1000.0; //we stop overtraining when reaching a low training error which is 1/1000 of the linear regression eror
	// use the weights and centroids from linear regression
	ALNNODE* pActiveLFN;
	pActiveLFN = pOT->GetTree();
	for(int m = 0; m < nDim -1; m++)
	{
		((pActiveLFN)->DATA.LFN.adblC)[m] = adblLRC[m]; // LR just means for linear regression
		((pActiveLFN)->DATA.LFN.adblW)[m+1] = adblLRW[m+1];
		fprintf(fpProtocol,"Weight on %s is %f centroid is %f\n", varname[nInputCol[m]] , adblLRW[m+1] ,adblLRC[m] ); 
	}
	((pActiveLFN)->DATA.LFN.adblC)[nDim - 1] = adblLRC[nDim - 1];
	((pActiveLFN)->DATA.LFN.adblW)[0] = adblLRW[0];
	((pActiveLFN)->DATA.LFN.adblW)[nDim] = -1.0;
	fprintf(fpProtocol,"Weight 0 is %f centroid %d is %f\n", adblLRW[0], nDim -1, adblLRC[nDim - 1] ); 
	dblFlimit = 0.1; // tiny value for the desired training error of overtraining
 	nNumberLFNs = 1;  // starts at 1 and should grow quickly
	dblLearnRate = 0.2; // try to correct all the error on first pass, one more pass to get better values
	// NB the View sets up the TRfile and the VARfile so nRowsTR is known here
	nEpochSize =  nRowsTR; // splitting occurs after an epoch of this size, but not for linear regression
	pOT->SetDataInfo(nEpochSize, nDim, NULL, NULL);
	(pOT->GetRegion(0))->dblSmoothEpsilon = 0.0; // Smoothing should never be used when estimating noise!
	for(int iteration = 0; iteration <20; iteration++)  // assume 20 to speed up
	{
    if(iteration == 3)
    {
			dblLearnRate = 0.1;
			fprintf(fpProtocol,"Learning rate changed to %f, Smoothing set to %f\n", dblLearnRate, (pOT->GetRegion(0))->dblSmoothEpsilon);
    }
    if(iteration == 10)// getting near the end of noise estimation, slowing down learning //was 46
    {
      dblLearnRate = 0.05;
			fprintf(fpProtocol,"Learning rate changed to %f\n", dblLearnRate);
    }
		bStopTraining = FALSE; // Should stop training when all pieces fit well. Doesn't work yet!
		// Jitter should never be used for overfitting, so we set this to false, but leave a TRUE value for approximation
		// call the training function
		nNumberEpochs = 10; //epochs per iteration, 20 iterations, total 200 epochs
		if(!pOT->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // overtraining
		{
		  fprintf(fpProtocol,"Overtraining failed!\n");
		}
		//fprintf(fpProtocol, "Iteration %d of overtraining continues. RMS Error %f\n", iteration, dblTrainErr);
		if (bStopTraining == TRUE)
		{
			fprintf(fpProtocol, "This overtraining stopped because all leaf nodes have stopped changing!\n");
			fprintf(fpProtocol, "\nOvertraining of an ALN completed at iteration %d \n", iteration);
			break;
		}
		fflush(fpProtocol);
	}
	fprintf(fpProtocol,"\nOvertraining of an ALN completed. Training RMSE = %f \n",dblTrainErr);
	// We are not finished with OTTS , we add OTVS and keep them for the noise level determination.
	// In future versions of the program we will create a weight-bounded ALN to learn the noise and store it for future evaluations
	// At present, there is no way to pass on the noise variance information
}

void ALNAPI approximate() // routine
{
	fprintf(fpProtocol, "\n**************Approximation with one or more ALNs begins ********\n");
	fflush(fpProtocol);
	int nalns = nALNs;  // The number of ALNs over which we average (for "bagging")
	double dblGlobalNoiseVariance;
	createTR_VARfiles(2);  // 2 prepares for using the whole TVfile and the whole VARfile with noise variance samples for training and stopping
	dblGlobalNoiseVariance = computeGlobalNoiseVariance();
	fprintf(fpProtocol,"Training %d approximation ALNs starts, using local noise variance to limit splitting\n", nalns);
	fprintf(fpProtocol, "The initial sample noise variance computed from all samples in the VARfile is %f\n", dblGlobalNoiseVariance);
	if(bJitter)
	{
		fprintf(fpProtocol, "Jitter is used during approximation\n");
	}
	else
	{
		fprintf(fpProtocol, "Jitter is not used during approximation\n");
	}
	const double adblFconstant[13]{ 9.00, 5.39, 4.11, 3.45, 3.05, 2.78, 2.59, 2.44, 2.32, 1.79, 1.61, 1.51, 1.40 };
	int dofIndex;
	dofIndex = nDim - 2; // the lowest possible nDim is 2 for one ALN input and one ALN output
	if(nDim > 10) dofIndex = 8;
	if(nDim > 20) dofIndex = 9;
	if(nDim > 30) dofIndex = 10;
	if(nDim > 40) dofIndex = 11;
	if(nDim > 60) dofIndex = 12;

	dblFlimit = adblFconstant[dofIndex]; // This can determine splitting for linear pieces
	//!!!!!!!!!!!!!!!!!
	dblFlimit = 2.59;  // Finally used for stopping splitting, but temporary until we can compute the dof for the pieces
	// have enough training samples and noise variance samples.
	// other values for dblFlimit with other numbers of samples, i.e. degrees of freedom, are:
	// n dblFlimit
	// 2 9.00
	// 3 5.39
	// 4 4.11
	// 5 3.45
	// 6 3.05
	// 7 2.78
	// 8 2.59
	// 9 2.44
	// 10 2.32
	// 20 1.79
	// 30 1.61
	// 40 1.51
	// 60 1.40
	// 120 1.26
	// Note: pieces don't split if they have enough training hits to determine them, so choose n >= nDim = number of ALN inputs plus 1 for the output.
  // REQUIRED IMPROVEMENT  We have to take into account the actual numbers of samples of TSfile and VARfile per linear piece of approximant.
	// As training of the approximant progresses, the dof of pieces decreases and dblFlimit should increase, depending on the minimum of samples per piece
	fprintf(fpProtocol, "nDim is %d and the F-limit used for stopping splitting is %f \n", nDim, dblFlimit);

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
		splitcontrol(apALN[n], dblFlimit); // here F-limit for splitting is set
		// NB it would be better to use the Flimit that is given by the hit counts on the piece that is to be split or not  IMPROVEMENT LATER
		int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO;
		// Set constraints on variables for ALN index n
		for(int m = 0; m < nDim-1; m++)
		{
			apALN[n]->SetEpsilon(adblEpsilon[m],m);
			// Expand the bounds on the variables by adding 10% of their stdevs.
			// These bounds will be needed for creating a multilevel DTREE and for inversion
			apALN[n]->SetMin(adblMinVar[m] - 0.1 * adblStdevVar[m],m);
			apALN[n]->SetMax(adblMaxVar[m] + 0.1 * adblStdevVar[m],m);
			// Bound weights
			// IS USING nDim -1 instead of nOuput A MISTAKE??
		  apALN[n]->SetWeightMin(- pow(2.0,0.5)* adblStdevVar[nDim - 1]/adblEpsilon[m],m); // changed 2009.09.09
		  apALN[n]->SetWeightMax( pow(2.0,0.5) * adblStdevVar[nDim - 1]/adblEpsilon[m],m); // reflects uniform dist. range
      // impose the a priori bounds on weights
      if(dblMinWeight[m] > apALN[n]->GetWeightMin(m))
      {
		    apALN[n]->SetWeightMin(dblMinWeight[m],m);
      }
      if(dblMaxWeight[m] < apALN[n]->GetWeightMax(m))
      {
		    apALN[n]->SetWeightMax(dblMaxWeight[m],m);
      }
		} // end of loop over m

		bTrainingAverage = FALSE;

    if(bClassify)
    {
      // in classification tasks we skip the linear regression and one ALN steps
      // so we have to put in an initial variance error estimate
      dblVarianceErr = 0.288;
    }
    apALN[n]->SetEpsilon(0,nDim -1); //  is this used? 
		// Tell the training algorithm the way to access the data using fillvector
    nEpochSize = nDim * 100; // changed 2009.11.06
		apALN[n]->SetDataInfo(nEpochSize,nDim,NULL,NULL);
		// >>>>>>>>>>>>> HERE IS WHERE SMOOTHING IS SET TO 0.5 OF APPROXIMATION TOLERANCE OR TOZERO BY CHANGING WHICH IS COMMENTED OUT   <<<<<<<<<<<<<<<<<<<<<<<
		// here is where Tolerance is used.  But maybe it should be determined by the noise!!!  The variance error can't be computed if we used the variance
		// set for noise samples!  Lots of work here.
		// (apALN[n]->GetRegion(0))->dblSmoothEpsilon = dblSmoothingFraction * dblTolerance;// approximation: we guess that smoothing, if used, should be about half the level of noise
		// for DEEP LEARNING, SMOOTHING SHOULD BE SET TO ZERO
		(apALN[n]->GetRegion(0))->dblSmoothEpsilon = 0.0;
		fprintf(fpProtocol, "The smoothing used for training each approximation ALN and the bagging average is %f\n", 0.0); 
		if(bEstimateRMSError)
		{
			// use the weights and centroids from linear regression
			ALNNODE* pActiveLFN;
			pActiveLFN = apALN[n]->GetTree();
			for(int m = 0; m < nDim - 1; m++)
			{
				((pActiveLFN)->DATA.LFN.adblC)[m] = adblLRC[m];
				((pActiveLFN)->DATA.LFN.adblW)[m+1] = adblLRW[m+1];
			}
			((pActiveLFN)->DATA.LFN.adblC)[nDim - 1] = adblLRC[nDim - 1];
			((pActiveLFN)->DATA.LFN.adblW)[0] = adblLRW[0];
			((pActiveLFN)->DATA.LFN.adblW)[nDim]= -1.0;
		}
		else
		{
			// Quickstart: get centroids and weights in neighborhood of good values
			dblMinRMSE = 0;
			dblLearnRate = 0.5; // we just want to get a linear piece near the data points
			nNumberEpochs = (int)floor(6.0 / dblLearnRate);
			bStopTraining = FALSE;
			if (!apALN[n]->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // fast change to start, dblLearnRate 1.0 for 2 epochs
			{
				fprintf(fpProtocol, "Training failed!\n");
			}
			if (bStopTraining == TRUE)
			{
				fprintf(fpProtocol, "This approximation training of ALN %d stopped because all leaf nodes have stopped changing!\n",n);
				fprintf(fpProtocol, "\nFirst step of training of approximation ALN %d completed \n", n);
				break;
			}
		}
		dblLearnRate = 0.2;
		nEpochSize = 3 * nRowsTR;
		// We stop when all pieces cease to allow splitting -- IF WE GET IT IMPLEMENTED!
		dblMinRMSE = 0.0; //temporary
    nNumberLFNs = 1;  // initialize at 1
		//int nOldNumberLFNs = 0; // for stopping criterion
    fprintf(fpProtocol,"----------  Training approximation ALN %d ------------------\n",n);
		for(int iteration = 0; iteration < 40; iteration++) // is 20 iterations enough?
		{
			fprintf(fpProtocol,"\nIteration %d of approximation with ALN %d\n", iteration,n);
      fflush(fpProtocol);
			// Call the training function
			bStopTraining = FALSE;
			if (!apALN[n]->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, bJitter, nNotifyMask))
			{
			  fprintf(fpProtocol,"Training failed!\n");
        exit(0);
			}
			if (bStopTraining == TRUE)
			{
				fprintf(fpProtocol, "This training stopped because all leaf nodes have stopped changing!\n");
				fprintf(fpProtocol, "\nOvertraining of approximation ALN number %d completed at iteration %d \n",n, iteration);
				break;
			}
			if(bEstimateRMSError == FALSE) // THIS IS IMPORTANT -- WE NEED TO HAVE A NEW STOPPING RULE
      {
				fprintf(fpProtocol,"Training RMSE = %f\n", dblTrainErr);
				if(dblTrainErr < adblEpsilon[nDim-1])
				{
					fprintf(fpProtocol,"Stopping approximation: Training RMSE = %f less than tolerance\n", dblTrainErr);
					fflush(fpProtocol);
					break;
				}
      }
      else
      {
				fprintf(fpProtocol,"Training RMSE = %f\n", dblTrainErr);
				if(dblTrainErr < adblEpsilon[nDim - 1]) // was dblTolerance
				{
					fprintf(fpProtocol,"Stopping approximation: Training RMSE = %f less than tolerance\n", dblTrainErr);
					fflush(fpProtocol);
					break;
				}
			}
			if(iteration == 10)
			{
				dblLearnRate = 0.1;
			}
			if(iteration < 25)
			{
				// Tree allowed to grow
				fprintf(fpProtocol,"Number of active LFNs = %d. Tree growing\n", nNumberLFNs);
				fflush(fpProtocol);
			}
			else
			{
				if(iteration == 24) dblLearnRate = 0.05;
				if(iteration == 29) dblLearnRate = 0.04;
				if(iteration == 34) dblLearnRate = 0.03;
				if(iteration == 39) dblLearnRate = 0.02;
				fprintf(fpProtocol,"Number of active LFNs = %d. Tree growth stopped. Next learn rate = %f\n", nNumberLFNs, dblLearnRate);				
			}
			//nOldNumberLFNs = nNumberLFNs;
			fflush(fpProtocol);
		} // end of loop of interations over one ALN
  } // end of the loop for n = ALN index
	// we don't destroy the ALNs because they are needed for further work in reporting
}

void ALNAPI outputtrainingresults() // routine
{
	fprintf(fpProtocol, "\n**** Analyzing results on the training/variance set begins ***\n");
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

void ALNAPI trainaverage() // routine
{
	// ************* TRAIN THE AVERAGE ALN **************************************
  int nalns;
  nalns = nALNs;
	bTrainingAverage = TRUE;
	// for training the average, we should use the TV set for points
	// since they best define the region where the data is
	// The values of those points are not used
	fprintf(fpProtocol,"\n********** Training an ALN by resampling the average ******\n");
	fprintf(fpProtocol, "The F-limit used for stopping splitting of the average ALN is %f \n", dblFlimit/nALNs);
	pAvgALN = new CMyAln; // NULL initialized ALN
	if (!(pAvgALN->Create(nDim, nDim-1) &&
				pAvgALN->SetGrowable(pAvgALN->GetTree())))
	{
    fprintf(fpProtocol,"Stopping: Growable average ALN creation failed!\n");
    exit(0);
	}
	splitcontrol(pAvgALN, dblFlimit / nALNs); // here is one place where the F-limit is set from F-tables, averaging reduces noise variance
	for(int k = 0; k < nDim - 1; k++) // do each variable k except the output
	{
		pAvgALN->SetEpsilon(adblEpsilon[k],k);
    // set the boundaries of the DTREE a bit beyond where the data samples are
		pAvgALN->SetMin(adblMinVar[k] - 0.1 * adblStdevVar[k],k);
		pAvgALN->SetMax(adblMaxVar[k] + 0.1 * adblStdevVar[k],k);
    // first put in the maximum and minimum weights that make sense
    // we expect about 10 samples across the box in dimension k, so
    // we use a factor of 0.1, which could be questioned
		pAvgALN->SetWeightMin(- pow(2.0,0.5) * adblStdevVar[nDim - 1]/adblEpsilon[k],k);
		pAvgALN->SetWeightMax( pow(2.0,0.5) * adblStdevVar[nDim - 1]/adblEpsilon[k],k);
    // now impose the a priori bounds on weights
    if(dblMinWeight[k] > pAvgALN->GetWeightMin(k))
    {
		  pAvgALN->SetWeightMin(dblMinWeight[k],k);
    }
    if(dblMaxWeight[k] < pAvgALN->GetWeightMax(k))
    {
		  pAvgALN->SetWeightMax(dblMaxWeight[k],k);
    }
	}
	(pAvgALN)->SetEpsilon(0,nDim -1); // small output tolerance based on variance error
	fprintf(fpProtocol,"Smoothing epsilon for training the average ALN is set the same as for approximation\n\n");
	fflush(fpProtocol);
	// Tell the training algorithm the way to access the data using fillvector
	dblLearnRate = 0.3;
	nEpochSize = (int)floor(1.199 / dblLearnRate) * nRowsTR;
	pAvgALN->SetDataInfo(nEpochSize,nDim,NULL,NULL);
	dblMinRMSE = 0;
		int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO;
	//*********
	if(bEstimateRMSError)
	{
		// Quickstart: get centroids and weights in neighborhood of good values
		nNumberEpochs = 10;
		if(!pAvgALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // fast change to start, dblLearnRate 1.0 for 2 epochs
		{
			 fprintf(fpProtocol,"Training failed!\n");
		}
	}
	else
	{
		// use the weights and centroids from linear regression
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
	if(bEstimateRMSError)
	{
		delete [] adblLRC; // delete this from free store
		delete [] adblLRW;
		adblLRC = NULL; // set the pointers to NULL
		adblLRW = NULL;
	}
	dblMinRMSE = 0.0; //THIS SHOULD BE REPLACED BY SOMETHING
	dblLearnRate = 0.2;
	nEpochSize = (int)floor(1.199 / dblLearnRate) * nRowsTR; //training average ALN
  nNumberLFNs = 1;  // initialize at 1
	//int nOldNumberLFNs = 1;  // having little change in nNumberLFNs is a stopping criterion
  for(int iteration = 0; iteration < 20; iteration++) // is 20 iterations enough?
	{
	  // Call the training function
		bStopTraining = TRUE;
	  if (!pAvgALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // average
	  {
		  fprintf(fpProtocol,"Average ALN training failed!\n");
      fflush(fpProtocol);
      exit(0);
	  }
		if (bStopTraining == FALSE)
		{
			fprintf(fpProtocol, "This training stopped because all leaf nodes have stopped changing!\n");
			fprintf(fpProtocol, "\nOvertraining of the average ALN completed at iteration %d \n", iteration);
			fflush(fpProtocol);
			break;
		}
		fprintf(fpProtocol,"\nIteration %d of training average ALN, RMSE = %f\n", iteration, dblTrainErr);
		fflush(fpProtocol);
		/* this may not be useful
		// we don't need to worry about overtraining because there is no noise
		if((double) nNumberLFNs < (double)nOldNumberLFNs * 0.8)
		{
			fprintf(fpProtocol,"Stopping training average ALN, number %d of active LFNs has shrunk too much\n", nNumberLFNs);
			fflush(fpProtocol);
			break;  // if active LFN growth has almost completely stopped, stop iterating
		}
		nOldNumberLFNs = nNumberLFNs; */
	}
	/*
	// polish the average ALN
  // keep the final value of nNumberLFNs and hence the same nEpochSize
	dblLearnRate = 0.05;
	if (!pAvgALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate,FALSE, nNotifyMask)) // polish average
	{
		fprintf(fpProtocol,"Polishing average ALN training failed!\n");
    fflush(fpProtocol);
    exit(0);
	}	
	else
	{
		fprintf(fpProtocol,"Polishing average ALN (ie at learning rate 0.05): RMSE = %f\n", dblTrainErr);
	}
	// now we have an average ALN
	if(bEstimateRMSError == FALSE)
	{
		fprintf(fpProtocol,"No variance done\n");
	}
	else
	{
		fprintf(fpProtocol,"Average ALN variance error w. r. t. variance set = %f\n", dblVarianceErr);
	}
	*/
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
    TRfile.Destroy();
    TSfile.Destroy();
    VARfile.Destroy();
    free(adblEpsilon);
		free(anInclude);
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
// This is called by the callback in myALN to fill in a data vector for the ALN.
{
	long nRow;
	nRow = (long)floor(ALNRandFloat() * (double) nRowsTR);
	for(int i = 0; i < nDim; i++)
	{
		adblX[i] = TRfile.GetAt(nRow,i,0);
	}
	if(bTrainingAverage) 
	{
		// To average several ALNs, we create new samples around those in TRfile and average the ALN outputs
		// The average ALN is created with negligible smoothing
		const ALNCONSTRAINT* pConstr;
		int nalns = nALNs;
		double dblValue;
		ALNNODE* pActiveLFN;
		for(int i = 0; i < nDim -1; i++)
		{
			// The idea here is the same as jitter
			// The chosen point triangularly distributed around the initial point
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

