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
#include <string.h>
#include <alnpp.h>
#include <dtree.h>
#include <datafile.h>
#include <malloc.h>
#include ".\cmyaln.h" 
#include "alnextern.h"
#include "alnintern.h"

// We use dblRespTotal in two ways and the following definition helps.
#define DBLNOISEVARIANCE dblRespTotal

// files used in training operations
CDataFile TRfile; // Training data, changed for different purposes.
long nRowsTR; // the number of rows in the current training set in TRfile

//routines
void ALNAPI createTR_file();
void ALNAPI createNoiseVarianceTool(); // This prepares to create the noise variance samples.
void ALNAPI approximate(); // Actually does training avoiding overtraining.
void ALNAPI outputTrainingResults(); // Shows how well approximate() has done.
void ALNAPI constructDTREE(int nMaxDepth); // Takes the average ALN and turns it into a DTREE for high speed evaluation.
void ALNAPI cleanup(); // Destroys ALNs and other objects when they are no longer needed.
void fillvector(double * adblX, CMyAln* paln); // An alternative way to select a vector for training from a file or an online stream.
double dist(double*, double*); // calculates the distance between domain points. This should be changed to reflect
// the rate of change of noise variance on the various axes so that the same distance implies the same change.

// ALN pointers
CMyAln* pBaseNeuron = NULL; // declares a pointer to an ALN. Now we only need one for this program for testing.

// Some global variables
double dblLimit; // If dblLimit <= 0, noise is estimated and used to stop splitting.  Otherwise it stops splitting.
double dblMinRMSE = 0; // Training is stopped when the mean square training error is smaller than this
double dblLearnRate = 0.2;  // Roughly, 0.2 corrects 20% of the deviation of ALN from desired. Fifteen passes through TRfile corrects most of the error.
int nMaxEpochs; // This controls the number of epochs between splittings of linear pieces. The linear pieces get time to fit better before splitting.
//double dblLimit = -1  ;// A negative value splits pieces based on an F test, otherwise they split if training MSE < dblLimit.
// MYTEST above now part of ALNDATAINFO
BOOL bStopTraining = FALSE; // Set to TRUE and becomes FALSE if any (active) linear piece still needs training
int nNotifyMask = AN_TRAIN | AN_EPOCH | AN_VECTORINFO; // Used with callbacks at different times for reporting on learning progress.
double * adblX = NULL; // This buffer holds an input vector including the desired output as last component.
double* aNoiseSampleTool = NULL; // aNoiseSampleTool helps create noise variance samples based on LFN weights during training.

using namespace std;

void ALNAPI createTR_file() // routine
{
	// This routine uses the TVfile to set up TRfile.
	// The V stands for validation, but we now no longer need a validation set.
	// The TVfile is all of the PreprocessedDataFile which is not used for testing..
	long i;
	int j;
	fprintf(fpProtocol, "Setting up the training data in TRfile\n");
	nRowsTR = TVfile.RowCount();
	TRfile.Create(nRowsTR, nDim);
	// First we fill TRfile from TVfile
	double dblValue;
	for (i = 0; i < nRowsTV; i++)
	{
		for (j = 0; j < nDim; j++)
		{
			dblValue = TVfile.GetAt(i, j, 0);
			TRfile.SetAt(i, j, dblValue, 0);
		}
	}
}

void ALNAPI createNoiseVarianceTool()
{
	fprintf(fpProtocol, "\n ********* Begin Creation of Noise Variance Tool ********\n");
	/*
	The array stores for each sample (X,y) in the training set, the vector
	from that sample to the L1-closest other sample and the difference of those sample values.
	These will be used later together with the weights of an LFN during training
	to construct a noise variance sample related to (X,y) on the LFN. Array aNoiseSampleTool is used
	to avoid overtraining.
	Potential speedup:
	We now go through all samples j to find one of the closest other samples to sample i.
	This could be greatly speeded up by some good nearest neighbor algorithm.
	*/
	// Set up the data
	createTR_file(); // This selects the training data after some samples for testing have been removed. 
	ASSERT(nRowsTR == TRfile.RowCount());
	double* aXcentral = NULL;
	double* aXnearby = NULL;
	aXcentral = (double*)malloc(nDim * sizeof(double));
	aXnearby = (double*)malloc(nDim * sizeof(double));
	aNoiseSampleTool = (double*) malloc(nRowsTR * nDim * sizeof(double));
	long i, j, iTimesnDim; // The TRfile can be huge. 
	int k;  // The dimension is small, eg less than 20 usually.
	double ds, dstemp;
	for (i = 0; i < nRowsTR; i++) // i is the central sample.
	{
		// Get the central sample aXcentral with index i.
		for (k = 0; k < nDim; k++)
		{
			aXcentral[k] = TRfile.GetAt(i, k, 0); // just the domain components
		}
		// We have to go through ALL samples j now to get the closest sample to i.
		ds = DBL_MAX;
		iTimesnDim = i * nDim;
		for (j = 0; j < nRowsTR; j++)
		{
			if (j != i) // sample i is not included
			{
				// Get the domain point for the j-th sample into aXnearby
				for (k = 0; k < nDim; k++)
				{
					aXnearby[k] = TRfile.GetAt(j, k, 0);
				}
				dstemp = dist(aXcentral, aXnearby);
				if(dstemp < ds)
				{
					ds = dstemp;
					// In either case, insert the new sample as closest yet found.
					// The first nDim -1 components are the vector, the last
					// is the value difference.
					for (int kk = 0; kk < nDim; kk++)
					{
						aNoiseSampleTool[iTimesnDim + kk] = aXnearby[kk] - aXcentral[kk];
					}
				}
			}
		} // end of j loop
	} // end loop over i
}

void ALNAPI approximate() // routine
{
	fprintf(fpProtocol, "\n**************Approximation with one or more ALNs begins ********\n");
	fflush(fpProtocol);
	fprintf(fpProtocol,"Training approximation ALN with the goal of avoiding overtraining\n");
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
	int nColumns = TRfile.ColumnCount(); // This is always nDim for training.
	ASSERT(nColumns == nDim);
	// ************ SET UP THE ALN FOR TRAINING **********
	// Set up the sample buffer.
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	// Set up the approximation ALN
	pBaseNeuron = new CMyAln; // NULL initialized ALN
	if (!pBaseNeuron->Create(nDim, nDim-1))
	{
	   fprintf(fpProtocol,"ALN creation failed!\n");
      fflush(fpProtocol);
			exit(0);
	}
	// Now make the tree growable
	if (!pBaseNeuron->SetGrowable(pBaseNeuron->GetTree()))		
	{
	  fprintf(fpProtocol,"Setting ALN growable failed!\n");
    fflush(fpProtocol);
    exit(0);
	}
	// Set constraints on variables for ALN
	// NB The following loop excludes the output variable of the ALN, index nDim -1.
	for (int m = 0; m < nDim - 1; m++)
	{
		pBaseNeuron->SetEpsilon(adblEpsilon[m], m);
		if (adblEpsilon[m] == 0)
		{
			fprintf(fpProtocol, "Stopping: Variable %d appears to be constant. Try removing it.\n", m);
			fflush(fpProtocol);
			exit(0);
		}
		// The minimum value of the domain is a bit smaller than the min of the data points
		// in TVfile, and the maximum is a bit larger.
		pBaseNeuron->SetMin(adblMinVar[m] - 0.1 * adblStdevVar[m], m);
		pBaseNeuron->SetMax(adblMaxVar[m] + 0.1 * adblStdevVar[m], m);

		// A rough value for the range of output (for a uniform dist.) divided by the
		// likely distance between samples in axis m.
		pBaseNeuron->SetWeightMin(-pow(3.0, 0.5) * adblStdevVar[nDim - 1] / adblEpsilon[m], m);
		pBaseNeuron->SetWeightMax(pow(3.0, 0.5) * adblStdevVar[nDim - 1] / adblEpsilon[m], m);

		// Impose the a priori bounds on weights which have been set by the user
		if (dblMinWeight[m] > pBaseNeuron->GetWeightMin(m))
		{
			pBaseNeuron->SetWeightMin(dblMinWeight[m], m);
		}
		if (dblMaxWeight[m] < pBaseNeuron->GetWeightMax(m))
		{
			pBaseNeuron->SetWeightMax(dblMaxWeight[m], m);
		}
	}
	(pBaseNeuron->GetRegion(0))->dblSmoothEpsilon = 0;
	fprintf(fpProtocol, "The smoothing for training approximation is %f\n", 0.0); 
	nMaxEpochs = 20; // The number of passes through the data without splitting LFNs.
	dblMinRMSE = 1e-20; // Stops training when the error is tiny.
	dblLearnRate = 0.2;
	bStopTraining = FALSE; // Set to TRUE at the start of each epoch in alntrain.cpp. 
	// Set FALSE by any piece needing more training. 
  nNumberLFNs = 1;  // initialize at 1
	// Tell the training algorithm the way to access the data using fillvector
	nRowsTR = TRfile.RowCount();
	ASSERT(nRowsTR == nRowsTV);
	const double* adblData = TRfile.GetDataPtr(); // This is where training gets samples.
	// The third parameter in the following could also set to NULL instead of adblData.
	// Then, instead of using FillInputVector(), the program uses fillvector()
	// for setting up the input vectors to the ALN.  fillvector() allows the
	// system to choose training vectors more flexibly (even online with proper programming).
	// The advantage of giving the pointer adblData instead of NULL is that training permutes
	// the order of the samples and goes through all samples exactly once per epoch.
	dblLimit = -1.0; // Negative to split leaf nodes according to an F test;
	// positive to split if training MSE > dblLimit.
	if (dblLimit <= 0)
	{
		fprintf(fpProtocol, "An F test is used to decide whether to split a piece depending on hit count. \n");
		bEstimateNoiseVariance = TRUE;
	}
	else
	{
		fprintf(fpProtocol, "A manually set limit, %f, is used to decide whether to split a piece. \n", dblLimit);
		bEstimateNoiseVariance = FALSE;
	}
	fflush(fpProtocol);
	pBaseNeuron->SetDataInfo(nRowsTR, nDim, adblData, NULL,dblLimit);
	fprintf(fpProtocol,"----------  Training approximation ALN ------------------\n");
	fflush(fpProtocol);
	for(int iteration = 0; iteration < 100; iteration++) 
	{
		fprintf(fpProtocol, "\nIteration %d of %d epochs ", iteration, nMaxEpochs);
		fflush(fpProtocol);
		// TRAIN ALNS WITHOUT OVERTRAINING   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
		if (!pBaseNeuron->Train(nMaxEpochs, dblMinRMSE, dblLearnRate, bJitter, nNotifyMask))
		{
			fprintf(fpProtocol,"Training failed!\n");
			fflush(fpProtocol);
      exit(0);
		}
		if (bStopTraining == TRUE)
		{
			fprintf(fpProtocol, "\nTraining of approximation ALN is complete after iteration %d \n", iteration);
			fprintf(fpProtocol, "All leaf nodes have stopped changing!\n");
			fflush(fpProtocol);
			break;
		}
		fprintf(fpProtocol, "Learning rate is %f\n", dblLearnRate);
		fflush(fpProtocol);
		if (iteration == 90) dblLearnRate = 0.15;
		if (iteration == 95) dblLearnRate = 0.05;
		if (iteration == 99) dblLearnRate = 0.01;
	} // end of loop of training interations over one ALN
	free(adblX);
	// we don't destroy the ALN because it is needed for further work in reporting
}

void ALNAPI outputTrainingResults() // routine
{
	fprintf(fpProtocol, "\n**** Analyzing results of approximation begins ***\n");
	// all the ALNs have been trained, now report results
	int i, j, k;
	double desired, sum;
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
		for (i = 0; i < nDim; i++)
		{
			adblX[i] = TVfile.GetAt(j, i, 0);
		}
		double dblValue = pBaseNeuron->QuickEval(adblX, &pActiveLFN);
		sum += dblValue;
		for (int k = 0; k < nDim; k++)
		{
			adblWAcc[k] += ((pActiveLFN)->DATA.LFN.adblW)[k + 1]; //the adblW vector has the bias in it
																											// so the components are shifted
			adblAbsWAcc[k] += fabs(((pActiveLFN)->DATA.LFN.adblW)[k + 1]);
		} ; 
		desired = TVfile.GetAt(j, nDim - 1, 0); // get the desired result	
		se += (sum - desired) * (sum - desired);
		if (fabs(desired - sum) > 0.5)  nClassError++; // desired must be integer
	}
	double rmse = sqrt(se / ((double)nRowsTV - 1.0)); // frees se for use below.
	// get the average weight on all variables k
	for (k = 0; k < nDim; k++)
	{
		adblWAcc[k] /= nRowsTV ;
		adblAbsWAcc[k] /= nRowsTV;
	}
	fprintf(fpProtocol, "Size of datasets PP TV Test %d  %d  %d \n", nRowsPP, nRowsTV, nRowsTS);
	fprintf(fpProtocol, "Root mean square error of ALN is %f \n", rmse);
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

void ALNAPI constructDTREE(int nMaxDepth) // routine
{
	// ******************  CONSTRUCT A DTREE FOR THE AVERAGE ALN *******************************
  fprintf(fpProtocol,"\n***** Constructing an ALN decision tree from the  ALN *****\n");
	DTREE* pBaseNeuronDTR;

	// Create a single-layer DTREE directly with ConvertDtree.
	// Setting nMaxDepth to a higher value than 1
	// requires a lot of time to generate the DTREE, but
	// allows much faster evaluation. This could turn out to be
	// useful for extremely demanding real-time tasks like
	// controlling nuclear fusion in ITER.
	pBaseNeuronDTR = pBaseNeuron->ConvertDtree(nMaxDepth);
  if(pBaseNeuronDTR == NULL)
	{
		fprintf(fpProtocol,"No DTREE was generated from the ALN. Stopping. \n");
    exit(0);
	}
	else
	{
		WriteDtree(szDTREEFileName,pBaseNeuronDTR);
		fprintf(fpProtocol,"The DTREE of the ALNs %s  was written.\n",szDTREEFileName );
	  fflush(fpProtocol);
  }
}

void ALNAPI cleanup() // routine
{
  if(bTrain)
  {
		// cleanup what was allocated for training in approximation
		pBaseNeuron->Destroy();
    // the TV file is not created for evaluation
		TVfile.Destroy();
		TSfile.Destroy();
    TRfile.Destroy();
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

void fillvector(double * adblX, CMyAln* paln) // routine
// This is called by the callback in CMyAln to fill in a data vector for training.
// If adblData is set by TRfile.GetDataPtr() and used in paln.SetDataInfo(nPoints,nCols,adblData)
// then the callback does not call this routine to get a training vector, it uses FillInputVector(...).
// To call this routine, set up the ALN for training with paln.SetDataInfo(nPoints,nCols,NULL).
{
	// IMPORTANT: This routine can be adapted to create training vectors
	// for real-time applications.
	long nRow;
	nRow = (long)floor(ALNRandFloat() * (double) nRowsTR); // This is where the TRfile is indicated for training.
	for(int i = 0; i < nDim; i++)
	{
		adblX[i] = TRfile.GetAt(nRow,i,0); // Notice that TRfile is fixed.
	}
}

double dist(double* adblA, double* adblB)
{
	// Computes the L1 distance between domain points. Since the domain axes k have different rates of change
	// of the noise variance, this metric gives a bound on the maximum change. In general, there will be
	// rates of change r(k): sum += r(k) * fabs(adblA[k] - adblB[k]);
	double sum = 0.0;
	for (int k = 0; k < nDim - 1; k++) // the sample value at nDim - 1 doesn't matter
	{
		sum += fabs(adblA[k] - adblB[k]);
	}
	return sum;
}
