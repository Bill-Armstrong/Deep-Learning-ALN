// ALN Library
// file train_ops.cpp
// Copyright (C) 1995 - 2010 William W. Armstrong.
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


// include classes
#include ".\cmyaln.h" 

//#include "alnfit.h"
#include "alnextern.h"
#include "alnintern.h"

int SplitDtree(DTREE** ppDest, DTREE* pSrc, int nMaxDepth);
static CMyAln* pAvgALN;      // an ALN representing the bagged average
static CMyAln** apALN = NULL;       // an array of pointers to ALNs
void ALNAPI dolinearregression();
void ALNAPI onealnfit();
void ALNAPI approximate();
void ALNAPI validate(CMyAln * pALN);
void fillvector(double * adblX, CMyAln* paln);
double	dblMinRMSE = 0;
double  dblLearnRate = 0.3;
int nNumberEpochs	= (int)floor(6.0 /dblLearnRate);

void ALNAPI dolinearregression() // routine
{
  fprintf(fpProtocol,"\n****** Finding an upper bound on error begins *****\n");
	// begin linear regression on a single ALN -- not a good method in general
	CMyAln * pALN = NULL;
	// Set up the ALN
	pALN = new CMyAln;
	if (!(pALN->Create(nDim, nDim-1)))
	{
		fprintf(fpProtocol,"Stopping: linear regression ALN creation failed!\n");
    fflush(fpProtocol);
		exit(0);
	}
  fflush(fpProtocol);
  // We don't set this tree growable -- no splitting
	// Set constraints on variables for the ALN	
	for( int m = 0; m < nDim-1; m++)
	{
		pALN->SetEpsilon(adblEpsilon[m],m);
    if(adblEpsilon[m] == 0)
    {
      fprintf(fpProtocol,"Stopping: Variable %d appears to be constant. Try removing it.\n",m);
      fflush(fpProtocol);
			exit(0);
    }
		pALN->SetMin(adblMinVar[m] - 0.1 * adblStdevVar[m],m);
		pALN->SetMax(adblMaxVar[m] + 0.1 * adblStdevVar[m],m);
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
  createTrainValfiles();
  bTrainingAverage = FALSE;
	double dblOldRmse = adblStdevVar[nDim - 1]; // this is a high initial value
  dblTolerance = 0.001 * adblStdevVar[nDim - 1]; // doesn't matter as there is only one piece // changed 2009.11.04
	pALN->SetEpsilon(dblTolerance,nDim -1);
 	(pALN->GetRegion(0))->dblSmoothEpsilon = 0.0; // linear regression: value doesn't matter since there are no MAX or MIN nodes
	int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO;
	// Tell the training algorithm the way to access the data using fillvector
  int nEpochSize = nDim * 100; // this is overkill, we just need an upper bound! (But it's fun to see it work.)
	pALN->SetDataInfo(nEpochSize,nDim,NULL,NULL);
	//********* dolinearregression
	// Quickstart: get centroids and weights in neighborhood of good values
	dblMinRMSE =  0; //don't stop early
	dblLearnRate = 0.2; // try to correct most of the error on the first few passes
  nNumberEpochs	= (int)floor(6.0 /dblLearnRate); // doing 30 epochs at a 0.2 learning rate
	if(!pALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // fast change to start, dblLearnRate 1.0 for 2 epochs
	{
		 fprintf(fpProtocol,"Training failed!\n");
	}	
	//**********
	dblLearnRate = 0.1;  //adjusted downwards after each iteration
	nNumberEpochs	= (int)floor(6.0 /dblLearnRate); // linear regression -- no splitting
	dblMinRMSE = 0;
	nNumberLFNs = 1;  // this is always correct here
	for(int iteration = 0; iteration < 10; iteration++)  // any number of iterations is enough; we'll get an upper bound no matter what
	{
		// Call the training function
		if (!pALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // linear regression
		{
			cerr << "Training failed!" << endl;
		}	
    fflush(fpProtocol);
    if(bEstimateRMSError)
		{
			validate(pALN);
			// is the new validation error less than the old one
			// if not, stop the iterations
			if((iteration > 3) && (dblValidationErr >= dblOldRmse)) // sooner or later this will happen change 2009.11.04
			{
				dblValidationErr = dblOldRmse; // return to the previous validation error and quit iterations
				fprintf(fpProtocol,"Stopping linear regression: validation error = %f has increased\n", dblValidationErr);
				break;
			}
			dblOldRmse = dblValidationErr;
			fprintf(fpProtocol,"\nIteration %d of linear regression continues, validation error = %f \n", iteration, dblValidationErr);
		}
		else
		{
			fprintf(fpProtocol, "\n Iteration %d of linear regression continues without validation\n", iteration);
			dblValidationErr = dblOldRmse; // just to have a value, even if only the stdev of output
		}
		// we adopt a constantly decreasing learning rate to increase accuracy
		dblLearnRate = 0.9 * dblLearnRate;
		nNumberEpochs	= (int)floor(6.0 /dblLearnRate);
	}
  // we should now have the best linear regression fit
  // find the weights on the linear piece using an evaluation at the 0 vector (could be anything).
	double * adblX = (double *) malloc((nDim) * sizeof(double));
	ALNNODE* pActiveLFN;
  for(int m=0; m < nDim; m++)
	{
		adblX[m] = 0;
	}
	double dummy = pALN->QuickEval(adblX, &pActiveLFN);
	fprintf(fpProtocol,"Linear regression weights on ALN\n");
	for (int m = 0; m < nDim-1; m++)
	{
    if(nLag[m] == 0)
    {
  		fprintf(fpProtocol,"Weight on %s is %f centroid is %f\n", varname[nInputCol[m]] , ((pActiveLFN)->DATA.LFN.adblW)[m+1], ((pActiveLFN)->DATA.LFN.adblC)[m]); 
		}
    else
    {
  		fprintf(fpProtocol,"Linear regression weight on ALN input %s@lag%d is %f centroid is %f\n", varname[nInputCol[m]], nLag[m] , ((pActiveLFN)->DATA.LFN.adblW)[m+1], ((pActiveLFN)->DATA.LFN.adblC)[m]); 
    }
		adblLRC[m] = ((pActiveLFN)->DATA.LFN.adblC)[m];
		adblLRW[m+1] = ((pActiveLFN)->DATA.LFN.adblW)[m+1];
  }
	fprintf(fpProtocol,"Weight 0 is %f centroid %d is %f\n", ((pActiveLFN)->DATA.LFN.adblW)[0], nDim -1, ((pActiveLFN)->DATA.LFN.adblC)[nDim - 1] ); 

	// now put the negative sum of weighted centroids into the 0 position
	// compress the weighted centroid info into W[0]       
	adblLRW[0] = adblLRC[nDim -1] = ((pActiveLFN)->DATA.LFN.adblC)[nDim -1];
  for (int m = 0; m < nDim - 1; m++)
  {
    adblLRW[0] -= adblLRW[m+1] * adblLRC[m];
  }
	//adblLRW[nDim] = ((pActiveLFN)->DATA.LFN.adblW)[nDim]; // this is just to check that the -1 weight was still there
	adblLRW[nDim] = -1.0;  // the -1 weight may make a difference somewhere, so we set it here
	if(bEstimateRMSError)
	{
		validate(pALN);
		fprintf(fpProtocol,"Linear regression final validation error = %f \n", dblValidationErr);
		dblLinRegErr = dblValidationErr;
		pALN->Destroy();
		// We are finished with that ALN and have destroyed it
		fprintf(fpProtocol,"Linear regression complete\n");
		fflush(fpProtocol);
	}
	else
	{
		fprintf(fpProtocol, "Linear regression validation error not determined\n");
		dblLinRegErr = 0; // the error is not evaluated
	}


}

void ALNAPI onealnfit() // this routine overfits one ALN with no smoothing so the ALNinterpolate the noisy data samples on the training set
{
	// begin fitting a single ALN using new training & validation sets
  // the difference from approximation is that the tolerance is frequently
  // adjusted in view of recent validation error, whereas in approximation
  // the tolerance is set low to begin with and kept constant
	fprintf(fpProtocol,"\n*** Overtraining a single ALN to find low validation error begins***\n");
	fprintf(fpProtocol,"Jitter is never used in this step \n");
	CMyAln * pALN;
	// Set up the ALN
	pALN = new CMyAln;
	if (!(pALN->Create(nDim, nDim-1)))
	{
		fprintf(fpProtocol,"Stopping: noise estimation ALN creation failed!\n");
    fflush(fpProtocol);
    exit(0);
	}
  if (!pALN->SetGrowable(pALN->GetTree()))		
	{
	  fprintf(fpProtocol,"Setting ALN growable failed!\n"); 
    fflush(fpProtocol);
    exit(0);
	}
	// Set constraints on variables for the ALN	
	for( int m = 0; m < nDim-1; m++)
	{
    if(adblEpsilon[m] == 0)
    {
	    fprintf(fpProtocol,"Stopping: Variable %d appears to be constant. Try removing it.\n",m); 
      fflush(fpProtocol);
      exit(0);
    }
		pALN->SetEpsilon(adblEpsilon[m],m);
		pALN->SetMin(adblMinVar[m] - 0.1 * adblStdevVar[m],m);
		pALN->SetMax(adblMaxVar[m] + 0.1 * adblStdevVar[m],m);
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
  createTrainValfiles(); // this means the validation results from linear regression and here can differ
	bTrainingAverage = FALSE;
	dblTolerance = 0.1 * dblValidationErr; // onealnfit: starts from 1/10 the linear regression validation error and gets smaller
	pALN->SetEpsilon(dblTolerance,nDim -1);
	(pALN->GetRegion(0))->dblSmoothEpsilon = 0.0; // for overfitting a single ALN, the smoothing should be zero so as not to interfere.
	int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO;
	fprintf(fpProtocol,"Initial validation error before onealnfit = %fl Tolerance = %f Smoothing = %f \n"
					, dblValidationErr, dblTolerance, (pALN->GetRegion(0))->dblSmoothEpsilon);
	// Tell the training algorithm the way to access the data using fillvector
	int nEpochSize  = nDim * 100;
	pALN->SetDataInfo(nEpochSize,nDim,NULL,NULL);
	//********* initial onealnfit weight and centroid components
	dblMinRMSE =  0; //don't stop early
	dblLearnRate = 0.3; // try to correct all the error on first pass, one more pass to get better values
	nNumberEpochs	= (int)floor(6.0 /dblLearnRate);  // doing just a few epochs at very high learning rate
	// use the weights and centroids from linear regression
	ALNNODE* pActiveLFN;
	pActiveLFN = pALN->GetTree();
	for(int m = 0; m < nDim -1; m++)
	{
		((pActiveLFN)->DATA.LFN.adblC)[m] = adblLRC[m];
		((pActiveLFN)->DATA.LFN.adblW)[m+1] = adblLRW[m+1];
		fprintf(fpProtocol,"Weight on %s is %f centroid is %f\n", varname[nInputCol[m]] , adblLRW[m+1] ,adblLRC[m] ); 

	}
	((pActiveLFN)->DATA.LFN.adblC)[nDim - 1] = adblLRC[nDim - 1];
	((pActiveLFN)->DATA.LFN.adblW)[0] = adblLRW[0];
	((pActiveLFN)->DATA.LFN.adblW)[nDim] = -1.0;
	fprintf(fpProtocol,"Weight 0 is %f centroid %d is %f\n", adblLRW[0], nDim -1, adblLRC[nDim - 1] ); 
	validate(pALN);
	fprintf(fpProtocol,"\nInitial validation error after using LR results = %fl \n", dblValidationErr);
	// We never skip validation for onealnfit()
	splitcontrol(pALN, 500.0);
		// if the validation error of a piece is a given
		// factor (e.g. 1.25) times greater than the
		// training error, splitting of that piece is stopped
		// the greater this number is, the longer splitting continues
		// we favor much more splitting in onealnfit() so we allow the validation error
		// to be *many* times the training error, i.e. we can overtrain to any extent.
	dblLearnRate = 0.2; //this will be diminished as overfitting progresses  TBD!
  nNumberEpochs	= (int)floor(6.0 /dblLearnRate); // it takes longer to stabilize before split at a lower learning rate
	nNumberLFNs = 1;  // starts at 1 and should grow quickly
	double nOldNumberLFNs = 1; // for stopping criterion
	for(int iteration = 0; iteration < 50; iteration++)  // assume 50 iterations is enough
	{
		fprintf(fpProtocol,"Iteration %d of noise estimation continues\n", iteration);
    if(iteration == 5)
    {
      dblLearnRate = 0.1;
			(pALN->GetRegion(0))->dblSmoothEpsilon = 0.0; // Smoothing should never be used when estimating noise!
			fprintf(fpProtocol,"Learning rate changed to %f, Smoothing set to %f\n", dblLearnRate, (pALN->GetRegion(0))->dblSmoothEpsilon);
    }
    if(iteration == 46)// getting near the end of noise estimation, slowing down learning
    {
      dblLearnRate = 0.05;
			nNumberEpochs	= (int)floor(6.0 /dblLearnRate);
			fprintf(fpProtocol,"Learning rate changed to %f\n", dblLearnRate);
    }
		if(iteration == 48)// the two last epochs are for stabilizing
    {
      dblLearnRate = 0.02;
			nNumberEpochs	= (int)floor(6.0 /dblLearnRate);
			fprintf(fpProtocol,"Learning rate changed to %f\n", dblLearnRate);
    }
		// Jitter should never be used for onealnfit, so we set this to false, but leave a TRUE value for approximation
		// call the training function
		if(!pALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // onealnfit
		{
		  fprintf(fpProtocol,"Training failed!\n");
		}	
		
		validate(pALN); //calculates the validation error
		splitcontrol(pALN, 500.0);
		dblTolerance = 0.1 * dblValidationErr; // set new tolerance low compared to validation error -- changed WWA 2009.10.09
		pALN->SetEpsilon(dblTolerance,nDim -1);
		if((iteration >= 4) && ((double) nNumberLFNs < (double)nOldNumberLFNs * 0.8))
		{
			fprintf(fpProtocol,"\nStopping noise estimation: the number %d of active linear pieces has shrunk too much \n", nNumberLFNs );
			fprintf(fpProtocol,"Noise estimation final validation error = %f\n", dblValidationErr);
			fflush(fpProtocol);
			break;
		}
		nOldNumberLFNs = nNumberLFNs;
		fprintf(fpProtocol,"Tolerance reduced to %f \n\n", dblTolerance);
		fflush(fpProtocol);
	}
	fprintf(fpProtocol,"\nNoise estimation training of single ALN: final validation error = %f \n", dblValidationErr);
	double dblTempFraction = pow((nDim + 1.0) / (nDim + 3.0), 0.5);
	dblTolerance = dblValidationErr * dblTempFraction;
	fprintf(fpProtocol, "We set the tolerance for approximation to a dimension-dependent fraction %f of the validation error to limit\n", dblTempFraction);
	fprintf(fpProtocol,"the splitting of linear pieces that fit better than this error level. \n");
	fprintf(fpProtocol,"ESTIMATE OF RMS NOISE IN THE DATA = %f WHICH MAY BE APPLIED TO SIMILAR DATA.\n", dblTolerance);
  fflush(fpProtocol);
  pALN->Destroy();
	// We are finished with that ALN and have destroyed it
}

void ALNAPI approximate() // routine
{
	int nalns = nALNs;  // The number of ALNs over which we average
	fprintf(fpProtocol,"\n**************Approximation with one or more ALNs begins ********\n");
	fprintf(fpProtocol,"Training %d approximation ALNs starts, using noise to limit splitting\n", nalns);
	if(bJitter)
	{
		fprintf(fpProtocol, "Jitter is used during approximation\n");
	}
	else
	{
		fprintf(fpProtocol, "Jitter is not used during approximation\n");
	}
	if(bEstimateRMSError)
	{
		fprintf(fpProtocol, "Validation done during approximation \n");
	}
	else
	{
		fprintf(fpProtocol, "Validation is skipped during approximation \n");
	}
	fprintf(fpProtocol,"Tolerance for training each approximation ALN is %f\n",dblTolerance);
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
    createTrainValfiles();  // we get different training sets at least for each ALN
		bTrainingAverage = FALSE;

    if(bClassify)
    {
      // in classification tasks we skip the linear regression and one ALN steps
      // so we have to put in an initial validation error estimate
      dblValidationErr = 0.288;
    }
    if(bEstimateRMSError == FALSE) // we already have the tolerance if we have estimated the RMS error
    {
      dblTolerance = dblSetTolerance;
    }
    apALN[n]->SetEpsilon(dblTolerance,nDim -1);
		// Tell the training algorithm the way to access the data using fillvector
    int nEpochSize = nDim * 100; // changed 2009.11.06
		apALN[n]->SetDataInfo(nEpochSize,nDim,NULL,NULL);
		// >>>>>>>>>>>>> HERE IS WHERE SMOOTHING IS SET TO 0.5 OF APPROXIMATION TOLERANCE OR TOZERO BY CHANGING WHICH IS COMMENTED OUT   <<<<<<<<<<<<<<<<<<<<<<<
		(apALN[n]->GetRegion(0))->dblSmoothEpsilon = dblSmoothingFraction * dblTolerance;// approximation: we guess that smoothing, if used, should be about half the level of noise
		// for DEEP LEARNING, sMOOTHING SHOULD BE SET TO ZERO
		fprintf(fpProtocol, "The smoothing used for training each approximation ALN and the bagging average is %f\n", dblSmoothingFraction*dblTolerance);
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
			if (!apALN[n]->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // fast change to start, dblLearnRate 1.0 for 2 epochs
			{
				fprintf(fpProtocol, "Training failed!\n");
			}
		}
		dblLearnRate = 0.2;
		nNumberEpochs	= (int)floor( 6.0/dblLearnRate);
		dblMinRMSE =  0.7 * dblTolerance; // stop early only if the error is smaller than the tolerance
    nNumberLFNs = 1;  // initialize at 1
		int nOldNumberLFNs = 0; // for stopping criterion
    fprintf(fpProtocol,"----------  Training approximation ALN %d ------------------\n",n);
		for(int iteration = 0; iteration < 40; iteration++) // is 20 iterations enough?
		{
			fprintf(fpProtocol,"\nIteration %d of approximation with ALN %d\n", iteration,n);
      fflush(fpProtocol);
			// Call the training function
			if (!apALN[n]->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, bJitter, nNotifyMask))
			{
			  fprintf(fpProtocol,"Training failed!\n");
        exit(0);
			}
			splitcontrol(apALN[n],2.0);
      if(bEstimateRMSError == FALSE)
      {
				fprintf(fpProtocol,"Training RMSE = %f\n", dblTrainErr);
				if(dblTrainErr < dblTolerance)
				{
					fprintf(fpProtocol,"Stopping approximation: Training RMSE = %f less than tolerance\n", dblTrainErr);
					fflush(fpProtocol);
					break;
				}
      }
      else
      {
				validate(apALN[n]);
				fprintf(fpProtocol,"Validation error = %f \nTraining RMSE = %f\n", dblValidationErr, dblTrainErr);
				if(dblTrainErr < dblTolerance)
				{
					fprintf(fpProtocol,"Stopping approximation: Training RMSE = %f less than tolerance\n", dblTrainErr);
					fflush(fpProtocol);
					break;
				}
			}
			if(iteration == 15)
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
				splitcontrol(apALN[n], 0.1); // this will stop splitting
				if(iteration == 24) dblLearnRate = 0.05;
				if(iteration == 29) dblLearnRate = 0.04;
				if(iteration == 34) dblLearnRate = 0.03;
				if(iteration == 39) dblLearnRate = 0.02;
				fprintf(fpProtocol,"Number of active LFNs = %d. Tree growth stopped. Next learn rate = %f\n", nNumberLFNs, dblLearnRate);				
			}
			nOldNumberLFNs = nNumberLFNs;
			fflush(fpProtocol);
		} // end of loop of interations over one ALN
  } // end of the loop for n = ALN index
	// we don't destroy the ALNs because they are needed for further work in reporting
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
	pAvgALN = new CMyAln; // NULL initialized ALN
	if (!(pAvgALN->Create(nDim, nDim-1) &&
				pAvgALN->SetGrowable(pAvgALN->GetTree())))
	{
    fprintf(fpProtocol,"Stopping: Average ALN creation failed!\n");
    exit(0);
	}
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
	if(bEstimateRMSError == FALSE)
  {
    dblTolerance = 0.7 * dblSetTolerance/ sqrt((double)nalns);
  }
  else
  {
    dblTolerance = 0.7 * dblValidationErr/ sqrt((double)nalns); 
		                                     // we take the final approximation validation error
		                                     // and correct for the expected increased accuracy of bagging
		                                     // to get the tolerance for training the average ALN
                                         // it seems 0.66 is too small here, leading to overtraining
                                         // while over .75 leads to too few pieces.  This is why we try 0.7
  }
	(pAvgALN)->SetEpsilon(dblTolerance,nDim -1); // small output tolerance based on validation error
	fprintf(fpProtocol,"Tolerance for training the average ALN is %f\n",dblTolerance);
	fprintf(fpProtocol,"Smoothing epsilon for training the average ALN is set the same as for approximation\n\n");
	// Tell the training algorithm the way to access the data using fillvector
	int nEpochSize = nDim * 100;
	pAvgALN->SetDataInfo(nEpochSize,nDim,NULL,NULL);
	dblMinRMSE = 0;
	dblLearnRate = 0.3;
	nNumberEpochs	= (int)floor(6.0 /dblLearnRate);
	int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO;
	//*********
	if(bEstimateRMSError)
	{
		// Quickstart: get centroids and weights in neighborhood of good values
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
	dblMinRMSE =  0.25 * dblTolerance; // stop if the error reaches 1/4 tolerance
	dblLearnRate = 0.2;
	nNumberEpochs	= (int) floor( 6.0/dblLearnRate); //training average ALN
  nNumberLFNs = 1;  // initialize at 1
	int nOldNumberLFNs = 1;  // having little change in nNumberLFNs is a stopping criterion
  for(int iteration = 0; iteration < 20; iteration++) // is 20 iterations enough?
	{
	  // Call the training function
	  if (!pAvgALN->Train(nNumberEpochs, dblMinRMSE, dblLearnRate, FALSE, nNotifyMask)) // average
	  {
		  fprintf(fpProtocol,"Average ALN training failed!\n");
      fflush(fpProtocol);
      exit(0);
	  }	
		fprintf(fpProtocol,"\nIteration %d of training average ALN, RMSE = %f\n", iteration, dblTrainErr);
		splitcontrol(pAvgALN, 10.0); // splitting is stopped if the training error of a piece times a factor 10 is less than tolerance
		// we don't need to worry about overtraining because there is no noise
		if((double) nNumberLFNs < (double)nOldNumberLFNs * 0.8)
		{
			fprintf(fpProtocol,"Stopping training average ALN, number %d of active LFNs has shrunk too much\n", nNumberLFNs);
			fflush(fpProtocol);
			break;  // if active LFN growth has almost completely stopped, stop iterating
		}
		nOldNumberLFNs = nNumberLFNs;
	}
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
		fprintf(fpProtocol,"No validation done\n");
	}
	else
	{
		validate(pAvgALN);
		fprintf(fpProtocol,"Average ALN validation error w. r. t. validation set = %f\n", dblValidationErr);
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

void ALNAPI outputtrainingresults() // routine
{
	fprintf(fpProtocol,"\n**** Analyzing results on the training/validation set begins ***\n");
  // all the ALNs have been trained, now report results
  int i,j,k,n;
	double desired, average, sum; 
	int nalns; // lower case indicates a value on the stack
  nalns = nALNs;
	ALNNODE* pActiveLFN = NULL;
  // test the average of the ALNs against data in the TV set
	double * adblX = (double *) malloc((nDim) * sizeof(double)); 
	double * adblWAcc = (double *) malloc((nDim) * sizeof(double));
	double * adblAbsWAcc = (double *) malloc((nDim) * sizeof(double));
  for(k=0; k < nDim; k++)
  {
    adblWAcc[k] = 0; // zero the accumulator for the average weight
    adblAbsWAcc[k] = 0; // zero the accumulator for the average absolute weight
  }
	double se = 0; // square error accumulator

	int	nClassError = 0;  // for classification problems
	for( j=0; j<nRowsTV; j++)
	{
		sum = 0;

		for (n = 0; n < nalns; n++)
		{
			for( i=0; i < nDim; i++)
			{
				adblX[i] = TVfile.GetAt(j,i,0);
			}
      double dblValue = apALN[n]->QuickEval(adblX, &pActiveLFN);
			sum += dblValue;
      for(int k=0; k < nDim; k++)
      {
        adblWAcc[k] += ((pActiveLFN)->DATA.LFN.adblW)[k+1]; //the adblW vector has the bias in it
                                                        // so the components are shifted
        adblAbsWAcc[k] += fabs(((pActiveLFN)->DATA.LFN.adblW)[k+1]); 
      }
		}
		average = sum / (double) nalns; // this is the result of averaging [0,1]-limited ALNs
		desired = TVfile.GetAt(j,nDim-1,0); // get the desired result	
		se += (average - desired) * (average - desired);
		if(((desired > 0.5) && (average < 0.5)) || ((desired < 0.5) && (average > 0.5))) nClassError++;
  }
	double rmse = sqrt(se / ((double)nRowsTV - 1.0)); // frees se for use below.
  // get the average weight on all variables k
  for(k=0; k < nDim; k++)
  {
    adblWAcc[k] /= (nRowsTV * nalns);
    adblAbsWAcc[k] /= (nRowsTV * nalns);
  }
	fprintf(fpProtocol,"Size of datasets PP TV Test %d  %d  %d \n", nRowsPP, nRowsTV, nRowsALNinputTestFile );
	fprintf(fpProtocol,"Root mean square error of the average over %d ALNS is %f \n", nalns, rmse);
	fprintf(fpProtocol,"Warning: the above result is optimistic, see results on the test set below\n");
	fprintf(fpProtocol,"Importance of each input variable:\n");
	fprintf(fpProtocol,"Abs imp = stdev(input var) * average absolute weight / stdev(output var) \n" );
  fprintf(fpProtocol,"Abs imp is numerical and indicates ups and downs in output when the given input varies.\n");
  fprintf(fpProtocol,"For example a sawtooth function with six teeth would have importance 12.\n");
  fprintf(fpProtocol,"First we have to compute the standard deviation of the output variable.\n");
  //compute the average of the output variable in the TVset
  k = nDim - 1; 
  desired = 0;
  for( j=0; j<nRowsTV; j++) 
  {
    desired += TVfile.GetAt(j,k,0);
  }
  desired /= nRowsTV; // now desired holds the average for variable k

  // compute the standard deviation of the output variable in the TVset
  se = 0;
  double temp;
  for( j=0; j < nRowsTV; j++) 
  {
    temp = TVfile.GetAt(j,k,0);
    se += ( temp - desired) * (temp - desired);
  }
  se /= ((double) nRowsTV - 1.0); // sample variance of the output variable
  double stdevOutput = sqrt(se);
  fprintf(fpProtocol,"\nStandard deviation of the output in the TVfile %f\n", stdevOutput);
  if(fabs(stdevOutput) < 1e-10)
  {
    fprintf(fpProtocol,"\nStopping: The standard deviation of the output on the TV set is near 0.\n");
    fclose(fpProtocol);
    exit(0);
  }
  // we compute the variance of each column of TV
  se = 0;
  for(k = 0; k < nDim - 1; k++) // do each variable k
  {
    //compute the average of variable k in TVset
    desired = 0;
    for( j=0; j<nRowsTV; j++) 
    {
      desired += TVfile.GetAt(j,k,0);
    }
    desired /= nRowsTV; // now desired holds the average for variable k
 
    // compute the standard deviation of variable k in TVset
    se = 0;
    double temp;
    for( j=0; j < nRowsTV; j++) 
    {
      temp = TVfile.GetAt(j,k,0);
      se += ( temp - desired) * (temp - desired);
    }
    se /= ((double) nRowsTV - 1.0); // sample variance of variable k
    dblImportance[k] = sqrt(se) * adblAbsWAcc[k]/stdevOutput;
    if(nLag[k] == 0)
    {
      fprintf(fpProtocol,"Variable %s: stdev = \t%f; avg.wt = \t%f; abs imp = \t%f\n",
      varname[nInputCol[k]], sqrt(se), adblWAcc[k], dblImportance[k]);
    }
    else
    {
      fprintf(fpProtocol,"Variable %s@lag%d: stdev = \t%f; avg.wt = \t%f; abs imp = \t%f\n",
      varname[nInputCol[k]], nLag[k], sqrt(se), adblWAcc[k], dblImportance[k]);
    }
    // we use the product of the variance of k and the average absolute weight as a measure of importance
	}
	if(bClassify)
	{
    fprintf(fpProtocol,"Number of TV file cases misclassified = %d out of %d\n",nClassError, nRowsTV);
    fprintf(fpProtocol,"Percentage of TV file cases misclassified = %f",100.0*(double) nClassError/(double)nRowsTV);
	}
  free(adblX); 
	free(adblWAcc);
	free(adblAbsWAcc);
	fflush(fpProtocol);
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
    ALNinputTestFile.Destroy();
    ALNinputValFile.Destroy();
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



void ALNAPI validate(CMyAln * pALN) // routine
{
	// validate the ALN using the validation file
	double * adblX = (double *) malloc((nDim) * sizeof(double));
	double desired = 0;
	double predict = 0;
  double      se = 0; // square error accumulator
	ALNNODE* pActiveLFN;
  // The following code has the goal of reducing the number of cycles
  // to do validation by sampling only nNumberLFNs * nDim * 128 samples.
  // This can make a big difference in speed when the validation file
  // ALNinputValFile is big.
  // If this sampling uses more than 50% of the validation set,
  // we use the whole validation set.
  // The number of LFNs is obtained from cmyaln at the end of each epoch
  // It must be initialized in the calling routine
  int n128TimesOversampled = nNumberLFNs * nDim * 128;
  int nSamplesUsed = n128TimesOversampled;
  BOOL b128TimesOversampled = TRUE; // indicates which sample size we are using
  if(n128TimesOversampled > 0.5 * nRowsALNinputValFile)
  {
    nSamplesUsed = nRowsALNinputValFile;
    b128TimesOversampled = FALSE;
  }
  dblValidationErr = 0; // makes sure we don't use an old value
  long k;
  if(b128TimesOversampled)
  {
	 for(long j = 0; j < nSamplesUsed; j++)
	 {
      // this could be done without replacement too
      k = (long)( ALNRandFloat()* (double) nRowsALNinputValFile);
      for(int i = 0; i < nDim; i++)
      {
	      adblX[i] = ALNinputValFile.GetAt(k,i,0);
      }
      adblX[nDim - 1] = 0; // not used in evaluation by QuickEval
      desired = ALNinputValFile.GetAt(k,nDim - 1,0); // get the desired result
      predict = pALN->QuickEval(adblX, &pActiveLFN);
      se += (predict - desired) * (predict - desired);
    } // end loop over samples (with replacement!)
  }
  else
  {
	 for(long j = 0; j < nSamplesUsed; j++)
	 {
		for(int i = 0; i < nDim; i++)
		{
		   adblX[i] = ALNinputValFile.GetAt(j,i,0);
		}
		adblX[nDim - 1] = 0; // not used in evaluation by QuickEval
		desired = ALNinputValFile.GetAt(j,nDim - 1,0); // get the desired result
		predict = pALN->QuickEval(adblX, &pActiveLFN);
      se += (predict - desired) * (predict - desired);
    } // end loop over validation set
  }
  dblValidationErr = sqrt(se / ((double)nSamplesUsed));
  fflush(fpProtocol);
  free(adblX);
}

void fillvector(double * adblX, CMyAln* paln) // routine
// This is called by the callback in myALN to fill in a data vector for the ALN.
{
	int nalns = nALNs; // put nALNs onto the stack
	ALNNODE* pActiveLFN;
	long nRow;
	double dblValue;
	nRow = (long)floor(ALNRandFloat() * (double) nRowsTR);
	for(int i = 0; i < nDim; i++)
	{
		adblX[i] = TRfile.GetAt(nRow,i,0);
	}
	if(bTrainingAverage) 
	{
		// resampling: we take a training point with inputs close to
		// the original data point inputs
		// the output value will be calculated at the moved point
		// and the new training point will be used to train the
		// average ALN
		// this is done with negligible smoothing so it makes sense even if
		// there is only one ALN in the average
		const ALNCONSTRAINT* pConstr;
		for(int i = 0; i < nDim -1; i++)
		{
			//adblX[i] += adblStdevVar[i] * (ALNRandFloat() - 0.5)/(double)nRowsTV;  // random variation of inputs
			// WWA changed 2009.10.05 to be the same as jitter except for the output variable
			// just stays close to original samples
			pConstr = paln->GetConstraint(i, 0);
			ASSERT(pConstr != NULL);
			adblX[i] += (ALNRandFloat() - ALNRandFloat()) * pConstr->dblEpsilon;
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

