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
extern BOOL bEstimateRMSError;
extern BOOL bStopTraining;

void dosplitcontrol(ALN* pALN, ALNNODE* pNode, double dblFlimit);
void dodivideTR(ALN* pALN, ALNNODE* pNode);
void dodivideVAR(ALN* pALN, ALNNODE* pNode);
void spliterrorsetTR(ALN * pALN);
void spliterrorsetVAR(ALN * pALN);
void dozerosplitvalues(ALN* pALN, ALNNODE* pNode);
void dozerosplitNOISEVARIANCE(ALN* pALN, ALNNODE* pNode);
int ALNAPI SplitLFN(ALN* pALN, ALNNODE* pNode);
// the following tells us how to use the SPLIT typedef between trainings of an ALN


void splitcontrol(ALN* pALN, double dblFlimit)  // routine
{
  ASSERT(pALN);
	ASSERT(pALN->pTree);
	// initializing all the SPLIT values to zero is done by void dozerosplitvalues(ALN* pALN, ALNNODE* pNode)
	if(dblFlimit >1.0)
	{
		// get square errors of pieces on training set
		spliterrorsetTR(pALN);
		// divide the training errors of the pieces by the hit counts, set counts to zero
		dodivideTR(pALN,pALN->pTree);
		// zero the NOISEVARIANCE component of split
		dozerosplitNOISEVARIANCE(pALN, pALN->pTree); 
		// get the values in VARfile which estimate the local noise variance
		spliterrorsetVAR(pALN);
		// divide the sum of local noise estimates on each piece by its count of hits
		dodivideVAR(pALN,pALN->pTree);
		//} I removed the test above so that we can get the noise samples analyzed
		// this conducts a search for all leaves in the ALN
	}
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



void dozerosplitNOISEVARIANCE(ALN* pALN, ALNNODE* pNode) // routine
{
	// initializes NOISEVARIANCE counters of all leaf nodes before processing of noise variance samples
	// Question: is this routine necessary after zeroing in dodivideTR(...)?
	ASSERT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		dozerosplitNOISEVARIANCE(pALN, MINMAX_LEFT(pNode));
		dozerosplitNOISEVARIANCE(pALN, MINMAX_RIGHT(pNode));
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		(pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE = 0;
	}
}



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
			if (nCountTemp >= nDim) // avoid division by 0 as well as having enough hits
			{
				(pNode->DATA.LFN.pSplit)->dblSqError /= nCountTemp;
				(pNode->DATA.LFN.pSplit)->nCount = 0; // after we get the value of the MSE, this can be zeroed.
				// nCount is used again for the noise variance samples in VARfile attached to this LFN
			}
			else
			{
				// This node can't be allowed to split any more because it has too few samples on it.
				pNode->fNode &= ~LF_SPLIT;
			}
			(pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE = 0; // This will be used in doing noise variance samples.
		}
	}
}


void dodivideVAR(ALN* pALN, ALNNODE* pNode) // routine
{
	// divides the total square error of each piece by its respective hit count
  ASSERT(pNode);
  if (NODE_ISMINMAX(pNode))
  {
    dodivideVAR(pALN, MINMAX_LEFT(pNode));
    dodivideVAR(pALN, MINMAX_RIGHT(pNode));
  }
  else
  {
    ASSERT(NODE_ISLFN(pNode));
		if (LFN_CANSPLIT(pNode)) // Don't do this leaf node if it can't split.
		{
			if((pNode->DATA.LFN.pSplit)->nCount >= nDim) // we need at least nDim samples to do stats
			{
				(pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE /=
											 (pNode->DATA.LFN.pSplit)->nCount;
			}
			else
			{
				// If there are too few or zero noise variance samples, make sure we don't divide
				LFN_FLAGS(pNode) &= ~LF_SPLIT;  // ~ is the unary one's complement operator.;
				       // setting the flag this way prevents this leaf node splitting.
			}
		}
	}
}



void dosplitcontrol(ALN* pALN, ALNNODE* pNode, double dblFlimit) // routine
{
	// This routine visits all the leaf nodes and determines whether or not to split.
	// During linear regression, there is no splitting anyway,so dblFlimit doesn't matter.
	// During overtraining the dblFlimit should be zero, which causes lots splitting
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
		if ( !LFN_CANSPLIT(pNode) || (NODE_RESPCOUNT(pNode) == 0)) // We have to allow pieces to share points
		{
			return;   // no splitting of this leaf node: already it can't split or has too few samples on it.
		}
		dblPieceSquareTrainError = (pNode->DATA.LFN.pSplit)->dblSqError; // average square error on the piece
		dblPieceNoiseVariance = (pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE; // average noise variance on the piece
		if ((dblFlimit > 1.0) && (dblPieceSquareTrainError < dblPieceNoiseVariance * dblFlimit))
		{
			// This implements the F-test criterion for stopping training (here equals stopping splitting).
			// If we get here, this piece is fitting within the noise variance by the F-test.
			// Since this piece fits well enough, we stop all future splitting of this leaf node (LFN).
			LFN_FLAGS(pNode) &= ~LF_SPLIT;  // this flag setting prevents further splitting
		}
		else
		{
			// If dblFlimit is >1, then the F-test has indicated more training is required for this piece.
			// If dblFlimit is 0, there is no stopping splitting unless there are too few points on the piece.
			// Ideally, we want the training error to go to zero.
			SplitLFN(pALN, pNode); // We split this and *every* leaf node that reaches this point.
			// We start off with bStopTraining == TRUE, but if any leaf node still needs more training, we set to FALSE
			bStopTraining = FALSE; // More training is required
		}
	}
}
