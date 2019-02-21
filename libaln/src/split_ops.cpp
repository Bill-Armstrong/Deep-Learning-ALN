// ALN Library (libaln)
// file split_ops.cpp
// 
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

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// include files
#include <stdafx.h>
#include <alnpp.h>
#include <datafile.h>


// include classes
#include ".\cmyaln.h"
#include "aln.h"

// We use dblRespTotal in two ways and the following definition helps.
#define DBLNOISEVARIANCE dblRespTotal

extern CDataFile TRfile; // The training file (global).
extern CDataFile NVfile;// The file containing noise variance samples. Its purpose is to prevent overtraining.
extern int nDim; // The dimension of the domain of the function to be learned plus one. Equals the number of ALN inputs .
extern double dblLimit; // The limit for deciding whether to split a piece.
extern BOOL bStopTraining; // This becomes TRUE and stops training when pieces are no longer splitting.
extern BOOL bTrainNV_ALN; // TRUE when a noise variance ALN is being trained.
void splitControl(ALN* pALN, double dblLimit);
void zeroSplitValues(ALN* pALN, ALNNODE* pNode);
void splitUpdateValues(ALN * pALN);
void doSplits(ALN* pALN, ALNNODE* pNode, double dblLimit);
int ALNAPI SplitLFN(ALN* pALN, ALNNODE* pNode);

// We use the first three fields in ALNLFNSPLIT (declared in aln.h)
// in two different ways: for training and between training intervals.
// The following routines use the SPLIT typedef between trainings, at the end of alntrain.cpp.
// Between nMaxEpochs epochs of training, where the hyperplanes in leaf nodes change weights,
// the pieces should be doing the best they can and decisions must be made on whether to split each leaf node. 
// To do that, splitcontrol, takes an average of the square training errors of
// each piece. This is compared to the average of the noise
// variance samples on the piece from NVfile.  If the average square training error is 
// greater than dblLimit times the average of the noise variance samples on the piece, then using
// an F-test, the piece is split because it does not yet fit within the limits of noise.

// Explanation of dblLimit
// dblLimit = 2.59 says that splitting of a linear piece is prevented when the mean square
// training error of a piece becomes less than 2.59 times the average of the noise variance
// samples on it. This value comes from tables of the F-test for d.o.f. > 7 and probability 90%.
// For 90% with 3 d.o.f the value is 5.39, i.e. with fewer d.o.f. training stops sooner
// and the training error will generally be larger than with a lower F-value.
// The values below for adblFconstant35 are interpolated and may not be accurate. 25 is the reciprocal of 75.

static const double adblFconstant90[13]{ 9.00, 5.39, 4.11, 3.45, 3.05, 2.78, 2.59, 2.44, 2.32, 1.79, 1.61, 1.51, 1.40 };
static const double adblFconstant75[13]{ 3.00, 2.36, 2.06, 1.89, 1.78, 1.70, 1.64, 1.59, 1.55, 1.36, 1.28, 1.24, 1.19 };
static const double adblFconstant50[13]{ 1,1,1,1,1,1,1,1,1,1,1,1,1 };
// The following two have not had any beneficial effect; but who knows when they might be useful?
static const double adblFconstant35[13]{ 0.58, 0.65, 0.70, 0.73, 0.75, 0.77, 0.78, 0.79, 0.80, 0.86, 0.88, 0.90, 0.92 };
static const double adblFconstant25[13]{ 0.333, 0.424, 0.485, 0.529, 0.562, 0.588, 0.610, 0.629, 0.645, 0.735, 0.781, 0.806, 0.840 };

void splitControl(ALN* pALN, double dblLimit)  // routine
{
  ASSERT(pALN);
	ASSERT(pALN->pTree);
	// initialize all the SPLIT values to zero
	zeroSplitValues(pALN, pALN->pTree);
	// get square errors of pieces on training set and the noise variance estimates from NVfile
	splitUpdateValues(pALN);
	// With the above statistics, dosplitcontol actually splits eligible pieces.
  doSplits(pALN, pALN->pTree, dblLimit);
  // Resetting the SPLIT components to zero by zeroSplitValues is done in alntrain.
}



// Routines that set some fields to zero

void zeroSplitValues(ALN* pALN, ALNNODE* pNode) // routine
{
	// initializes split counters of all leaf nodes before the next training period
	ASSERT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		zeroSplitValues(pALN, MINMAX_LEFT(pNode));
		zeroSplitValues(pALN, MINMAX_RIGHT(pNode));
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		(pNode->DATA.LFN.pSplit)->nCount = 0;
		(pNode->DATA.LFN.pSplit)->dblSqError = 0;
		(pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE = 0;
	}
}

// Routines that get the training errors and noise variance values.

void splitUpdateValues(ALN * pALN) // routine
{
	// Assign the square errors on the training set and the noise variance
	// sample values to the leaf nodes of the ALN.
	// IMPORTANT: THE SAMPLE DOMAINS IN TRfile AND NVfile MUST BE IN 1-TO-1 CORRESPONDENCE
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	double fromFile = 0;
	double predict = 0;
	long nDimm1 = nDim - 1;

	ALNNODE* pActiveLFN;
	long nrows = TRfile.RowCount();
	for (long i = 0; i < nrows; i++) // TRfile and NVfile have to have a 1-to-1 correspondence of samples
	{
		for (int j = 0; j < nDim; j++)
		{
			adblX[j] = TRfile.GetAt(i, (long) j, 0);
		}
		predict = ALNQuickEval(pALN, adblX, &pActiveLFN); // the current ALN value
		if (LFN_CANSPLIT(pActiveLFN)) // Skip this leaf node if it can't split anyway.
		{
			fromFile = adblX[nDimm1]; //adblX[nDim - 1] is not used in evaluation by ALNQuickEval
			(pActiveLFN->DATA.LFN.pSplit)->nCount++;
			(pActiveLFN->DATA.LFN.pSplit)->dblSqError += (predict - fromFile) * (predict - fromFile);
			fromFile = NVfile.GetAt(i, nDimm1, 0);
			(pActiveLFN->DATA.LFN.pSplit)->DBLNOISEVARIANCE += fromFile; // A noise variance sample value for the same X (same index)
		}
	} // end loop over both files
	free(adblX);
} // END of splitUpdateValues

void doSplits(ALN* pALN, ALNNODE* pNode, double dblLimit) // routine
{
	// This routine visits all the leaf nodes and determines whether or not to split.
	// If dblLimit < 0, it uses an F test with d.o.f. based on the number of samples counted,
	// but if dblLimit >= 0 it uses the actual dblLimit value to compare to the square training error.

	ASSERT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		doSplits(pALN, MINMAX_LEFT(pNode), dblLimit);
		doSplits(pALN, MINMAX_RIGHT(pNode), dblLimit);
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		if (LFN_CANSPLIT(pNode))
		{
			long Count = (pNode->DATA.LFN.pSplit)->nCount;
			if (Count >= nDim) // Don't consider splitting if there are too few samples; do nothing.
			{
				double dblPieceSquareTrainError = (pNode->DATA.LFN.pSplit)->dblSqError; // total square error on the piece
				double dblPieceNoiseVariance = (double)Count; // These two values will be changed to do the F test
				double dblSplitLimit = dblLimit; // if dblLimit is <= 0, otherwise they test training MSE < dblLimit
				if (dblLimit <= 0) // if this is TRUE, we do the F test.
				{
					dblPieceNoiseVariance = (pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE; // total noise variance samples
					int dofIndex; // get the dblSplitLimit corresponding to the degrees of freedom of the F test
					dofIndex = Count - 2;
					if (Count > 10) dofIndex = 8;
					if (Count > 20) dofIndex = 9;
					if (Count > 30) dofIndex = 10;
					if (Count > 40) dofIndex = 11;
					if (Count > 60) dofIndex = 12;
					dblSplitLimit = adblFconstant50[dofIndex]; // One can reject the H0 of a good fit with various percentages
					// 90, 75, 50, 35, 25. E.g. 90% says that if the training error is greater than the dblSplitLimit prescribes
					// it is 90% sure that the fit is bad.  A higher percentage needs less training time.
					// Note that when there are few hits on the piece, the dblSplitLimit is larger and 
					// the criterion for fitting well enough is easier to satisfy.
					if (bTrainNV_ALN) // A special case of doing the F-test with a different NVfile
					{
						// SNV = dblPieceNoiseVariance is the sum of NV sample values over the piece
						// The training SSE is compared to SNV; the counts are equal;
						dblSplitLimit = adblFconstant50[dofIndex];
					}
				}
				if (dblPieceSquareTrainError <= dblPieceNoiseVariance * dblSplitLimit)
				{
					// The piece fits well, stop splitting it. 
					LFN_FLAGS(pNode) &= ~LF_SPLIT;  // this flag setting prevents further splitting 
				}
				else
				{
					// The piece doesn't fit and needs to split; then training must continue.
					SplitLFN(pALN, pNode); // We split *every* leaf node that reaches this point.
					// We start an epoch with bStopTraining == TRUE, but if any leaf node might still split,
					bStopTraining = FALSE; //  we set it to FALSE and continue to another epoch of training.
				}
			} // do nothing if Count < nDim, the piece is underdefined.
		}
	}
}
