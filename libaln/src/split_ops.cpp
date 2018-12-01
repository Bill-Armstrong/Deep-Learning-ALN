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

extern double dblMax;
extern double* adblEpsilon;
extern double dblFlimit;
extern int nDim;
extern double dblSetTolerance;
extern BOOL bEstimateNoiseVariance;
extern BOOL bStopTraining;
extern CDataFile VARfile;
extern CDataFile TRfile;

void dosplitcontrol(ALN* pALN, ALNNODE* pNode, double dblFlimit);
void dozerosplitvalues(ALN* pALN, ALNNODE* pNode);
void spliterrorsetTR(ALN * pALN);
void dodivideTR(ALN* pALN, ALNNODE* pNode);
void splitNoiseSetVAR(ALN * pALN);
void dodivideVAR(ALN* pALN, ALNNODE* pNode);
int ALNAPI SplitLFN(ALN* pALN, ALNNODE* pNode);
// the following routines use the SPLIT typedef between trainings of an ALN


void splitcontrol(ALN* pALN, double dblFlimit)  // routine
{
  ASSERT(pALN);
	ASSERT(pALN->pTree);
	// initialize all the SPLIT values to zero
	dozerosplitvalues(pALN, pALN->pTree);
	// get square errors of pieces on training set
	spliterrorsetTR(pALN);
	// divide the training errors of the pieces by the hit counts and
	// set nCount component of split to zero
	dodivideTR(pALN,pALN->pTree);
	if (dblFlimit > 1.0)
	{
		// get the values in VARfile which estimate the local noise variance
		splitNoiseSetVAR(pALN);
		// divide the sum of local noise estimates on each piece by its count of hits
		// now done in dosplitcontrol along with the F-test
	}
	// splitcontrol also splits with dblLimit <= 1.0 without noise variance
	// With the above statistics, dosplitcontol actually splits eligible pieces.
  dosplitcontrol(pALN, pALN->pTree, dblFlimit);
  // reset the SPLIT components to zero
	//dozerosplitvalues(pALN, pALN->pTree); This is done in alntrain.
}

// We use the first three fields in ALNLFNSPLIT (declared in aln.h)
// in two different ways: first for training and second between training intervals.
// Between trainings, where the hyperplanes in leaf nodes change weights,
// those changes stop and decisions must be made on whether or not to split each leaf node. 
// To do that, this routine, splitcontrol, takes an average of the square errors of
// each piece on the training set hits. This is compared to the average of the noise
// variance samples on the piece from VARfile.  If the average square training error is 
// greater than a specified fraction of the average of the noise variance samples on the piece, then according to
// an F-test with limit dblFlimit, the piece is not split. It does not yet fit within the limits of noise.
// We want to implement a new idea: when a linear piece does not satisfy the criterion
// for fitting within the limits of noise, and a whole training epoch has passed without
// significant progress towards fitting, we do split it so training can continue.
// The hope is that eventually all pieces will satisfy the F-test criterion for fitting
// and the training of the ALN can be stopped.

// Routines that set some fields to zero

void dozerosplitvalues(ALN* pALN, ALNNODE* pNode) // routine
{
	// initializes split counters of all leaf nodes before the next training period
	ASSERT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		dozerosplitvalues(pALN, MINMAX_LEFT(pNode));
		dozerosplitvalues(pALN, MINMAX_RIGHT(pNode));
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

void spliterrorsetTR(ALN * pALN) // routine
{
	// assign the square errors on the training set to the leaf nodes of the ALN
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	double desired = 0;
	double predict = 0;
	double      se = 0; // square error accumulator
	ALNNODE* pActiveLFN;
	for (long j = 0; j < nRowsTR; j++)
	{
		for (int i = 0; i < nDim; i++)
		{
			adblX[i] = TRfile.GetAt(j, i, 0);
		}
		predict = ALNQuickEval(pALN, adblX, &pActiveLFN); // the current ALN value
		if (LFN_CANSPLIT(pActiveLFN)) // Skip this leaf node if it can't split anyway.
		{
			desired = adblX[nDim - 1]; //adblX[nDim - 1] is not used in evaluation by ALNQuickEval
			se = (predict - desired) * (predict - desired);
			(pActiveLFN->DATA.LFN.pSplit)->nCount++;
			(pActiveLFN->DATA.LFN.pSplit)->dblSqError += se;
		}
	} // end loop over TRset
	free(adblX);
} // END of spliterrorsetTR

void dodivideTR(ALN* pALN, ALNNODE* pNode) // routine
{
	// divides the total square errors of the pieces by their hit count
	ASSERT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		dodivideTR(pALN, MINMAX_LEFT(pNode));
		dodivideTR(pALN, MINMAX_RIGHT(pNode));
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		if (LFN_CANSPLIT(pNode)) // skip this leaf node if it has stopped splitting
		{
			long nCountTemp = (pNode->DATA.LFN.pSplit)->nCount;
			if (nCountTemp > 0) // avoid division by 0 as well as having enough hits
			{
				(pNode->DATA.LFN.pSplit)->dblSqError /= nCountTemp;
				(pNode->DATA.LFN.pSplit)->nCount = 0; // after we get the value of the MSE,
				// the nCount is zeroed. nCount is used again for the noise variance samples
				// in VARfile attached to this LFN
			}
		}
	}
}

void splitNoiseSetVAR(ALN * pALN) // routine
{
	// NB  It might be possible to fuse this with spliterrorsetTR
	// but then we couldn't do several noise decompositions
	// and get a lot more noise variance samples in the future
	// assign the noise variance samples to the leaf nodes of the ALN and add them up
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	double desired = 0;
	double value = 0;
	ALNNODE* pActiveLFN;
	double      se = 0; // sample value accumulator for LFN DBLNOISEVARIANCE
	long nNoiseVarianceSamples = VARfile.RowCount();
	for (long i = 0; i < nNoiseVarianceSamples; i++)
	{
		for (int j = 0; j < nDim; j++)
		{
			adblX[j] = VARfile.GetAt(i, j, 0); // the value at nDim - 1 is used only for desired
		}
		if (adblX[nDim - 1] > 0) // noise variance samples can be 0 -- if too few samples - skip
		{
			value = ALNQuickEval(pALN, adblX, &pActiveLFN);  // all that matters is which LFN the X-vector lies on
			{
				(pActiveLFN->DATA.LFN.pSplit)->nCount++; // we have to zero this before this routine is called.
				(pActiveLFN->DATA.LFN.pSplit)->DBLNOISEVARIANCE += adblX[nDim - 1];
					//add the value of the noise variance sample
			}
		}
	} // end loop over VARfile
	free(adblX);
} // END of splitNoiseSetVAR

void dosplitcontrol(ALN* pALN, ALNNODE* pNode, double dblFlimit) // routine
{
	// This routine visits all the leaf nodes and determines whether or not to split.
	// During linear regression, there is no splitting anyway,so dblFlimit doesn't matter.
	// During overtraining the dblFlimit should be zero, which causes a lot of splitting
	// until the pieces all have nDim samples defining them (at least that's the hope!).
	double dblPieceSquareTrainError;
	double dblPieceNoiseVariance;
	ASSERT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		dosplitcontrol(pALN, MINMAX_LEFT(pNode), dblFlimit);
		dosplitcontrol(pALN, MINMAX_RIGHT(pNode), dblFlimit);
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		if (LFN_CANSPLIT(pNode))
		{
			dblPieceSquareTrainError = (pNode->DATA.LFN.pSplit)->dblSqError; // average square error on the piece
			if (((pNode->DATA.LFN.pSplit)->nCount > 0)) // For overtraining, the count is zero and dblFlimit <= 1.0
			{
				dblPieceNoiseVariance = (pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE / (pNode->DATA.LFN.pSplit)->nCount;
			}
			else
			{
				dblPieceNoiseVariance = 1.0;
			}
			// average noise variance on the piece or 1.0 as the case may be
			if (dblPieceSquareTrainError < dblPieceNoiseVariance * dblFlimit)
			{
				// This implements the F-test criterion for stopping training (here stopping splitting).
				// If we get here, this piece is fitting within the noise variance by the F-test.
				// Splitting is also stopped if noise variance is zero (e.g. overtraining) and the training MSE < dblFlimit.
				// Since this piece fits well enough, we stop all future splitting of this leaf node (LFN).
				LFN_FLAGS(pNode) &= ~LF_SPLIT;  // this flag setting prevents further splitting
			}
			else
			{
				// If dblFlimit < 1.0, or F-test shows poor fit, LFN splits & training continues
				SplitLFN(pALN, pNode); // We split *every* leaf node that reaches this point.
				// We start an epoch with bStopTraining == TRUE, but if any leaf node might still split,
				bStopTraining = FALSE; //  we set it to FALSE and continue to another epoch of training.
			}
		}
	}
}
