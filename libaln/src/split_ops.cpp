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

extern double dblMax;
extern double* adblEpsilon;
extern double dblFlimit;
extern int nDim;
extern double dblSetTolerance;
extern long nRowsTR;
extern long nRowsVAR;
extern CDataFile TRfile;
extern CDataFile VARfile;
extern BOOL bEstimateRMSError;
BOOL bStopTraining = FALSE;


void dosplitcontrol(ALN* pALN, ALNNODE* pNode, double dblFlimit);
void dodivideTR(ALN* pALN, ALNNODE* pNode);
void dodivideVAR(ALN* pALN, ALNNODE* pNode);
void spliterrorsetTR(ALN * pALN);
void spliterrorsetVAR(ALN * pALN);
void dozerosplitvalues(ALN* pALN, ALNNODE* pNode);
void dozerosplitNOISEVARIANCE(ALN* pALN, ALNNODE* pNode);
int ALNAPI SplitLFN(ALN* pALN, ALNNODE* pNode);
// the following tells us how to use the SPLIT typedef between trainings of an ALN
#define DBLNOISEVARIANCE dblRespTotal

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
  dosplitcontrol(pALN, pALN->pTree, dblFlimit);
  // reset the SPLIT components to zero
	//dozerosplitvalues(pALN, pALN->pTree); RESET COUNTERS???
}

// We use the first three fields in ALNLFNSPLIT (declared in aln.h)
// in two different ways: first for training and second between training intervals.
// Between trainings, where the hyperplanes in leaf nodes change weights,
// those changes stop and decisions must be made on whether or not to split each leaf node. 
// To do that, this routine, splitcontrol, takes an average of the square errors of
// each piece on the training set hits. This is compared to the average of the noise
// variance samples from VARfile.  If the average square training error is greater than a specified
// fraction of the average of the noise variance samples on the piece, then according to
// an F-test with limit dblFlimit, the piece is not split. It does not yet fit within the limits of noise.
// We want to implement a new idea: when a linear piece does not satisfy the criterion
// for fitting within the limits of noise, and a whole training epoch has passed without
// significant progress towards splitting, we do split it so training can continue.
// The hope is that eventually all pieces will satisfy the F-test criterion for fitting
// and the training of the ALN can be stopped.
// Another idea is to reduce the dblFlimit in stages, but we just have to see if this can do anything useful.


void dozerosplitvalues(ALN* pALN, ALNNODE* pNode) // routine
{
	// initializes counters of all leaf nodes before the next training period
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
				// this node can't be allowed to split any more
				pNode->fNode &= ~LF_SPLIT;
			}
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
		if (LFN_CANSPLIT(pNode)) // only do this leaf node if it has not stopped splitting
		{
			if((pNode->DATA.LFN.pSplit)->nCount >= nDim) // we need at least nDim samples to do stats
			{
				(pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE /=
											 (pNode->DATA.LFN.pSplit)->nCount;
			}
			else
			{
				// if there are too few or no variance samples, make sure we don't divide
				LFN_FLAGS(pNode) &= ~LF_SPLIT;  // ~ is the unary one's complement operator; this flag prevents splitting
				(pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE = dblMax;
			}
		}
	}
}



void dosplitcontrol(ALN* pALN, ALNNODE* pNode, double dblFlimit) // routine
{
	// This routine visits all the leaf nodes and determines whether or not to split.
	// During linear regression, there is no splitting anyway.
	// During overtraining the dblFlimit is should be zero, which causes lots of splitting.
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
		if ( !LFN_CANSPLIT(pNode) || (NODE_RESPCOUNT(pNode) < pALN->nDim))
		{
			return;   // no splitting of this leaf node: already can't split or too few samples
		}
		dblPieceSquareTrainError = (pNode->DATA.LFN.pSplit)->dblSqError; // average square error on the piece
		dblPieceNoiseVariance = (pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE; // average noise variance on the piece


		if ((dblFlimit > 1.0) && (dblPieceSquareTrainError < dblPieceNoiseVariance * dblFlimit)) // this implements the F-test criterion for stopping training
		{
			// If we get here, this piece is fitting within the noise variance by the F-test.
			// Since this piece fits well enough, we stop all future splitting of this leaf node (LFN).
			LFN_FLAGS(pNode) &= ~LF_SPLIT;  // this flag setting prevents further splitting
		}
		else
		{
			// If dblFlimit is >1, then the F-test has indicated more training is required.
			// If dblFlimit is 0, there is no stopping unless there are too few points on the piece.
			// We want the training error to go to almost zero.
			SplitLFN(pALN, pNode); // we split *every* leaf node that has been training
			// we start off with bStopTraining == FALSE, but if any leaf node needs to continue we set to true
			bStopTraining = FALSE; // if any leaf node gets here, then more training is required
		}
	}
}

void spliterrorsetTR(ALN * pALN) // routine
{
	// assign the square errors on the training set to the leaf nodes of the ALN
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	double desired = 0;
	double predict = 0;
	double      se = 0; // square error accumulator
	ALNNODE* pActiveLFN;
	for (int j = 0; j < nRowsTR; j++)
	{
		for (int i = 0; i < nDim; i++)
		{
			adblX[i] = TRfile.GetAt(j, i, 0);
		}
		desired = adblX[nDim - 1]; // get the desired result
		adblX[nDim - 1] = 0; // not used in evaluation by ALNQuickEval
		predict = ALNQuickEval(pALN, adblX, &pActiveLFN);
	
		//if (LFN_ISINIT(pActiveLFN)) // skip this leaf node if it has stopped training
		{
			se = (predict - desired) * (predict - desired);
			(pActiveLFN->DATA.LFN.pSplit)->nCount++;
			(pActiveLFN->DATA.LFN.pSplit)->dblSqError += se;
		}
	} // end loop over TRset
	free(adblX);
} // END of spliterrorsetTR


void spliterrorsetVAR(ALN * pALN) // routine
{
		// NB  It might be possible to fuse this with spliterrorsetTR but then we couldn't do several noise decompositions
		// and get a lot more noise variance samples in the future
		// assign the noise variance samples to the leaf nodes of the ALN and add them up
		double * adblX = (double *)malloc((nDim) * sizeof(double));
		double desired = 0;
		double value = 0;
		ALNNODE* pActiveLFN;
		double      se = 0; // sample value accumulator for LFN DBLNOISEVARIANCE
		for (int j = 0; j < nRowsVAR; j++)   
		{
			for (int i = 0; i < nDim; i++)
			{
				adblX[i] = VARfile.GetAt(j, i, 0); // the value at nDim - 1 is used only for desired
			}
			// pAln has to be the current approximant! Is this correct?
			value = ALNQuickEval(pALN, adblX, &pActiveLFN);  // all that matters is which LFN the X-vector lies on
			//if (LFN_ISINIT(pActiveLFN)) // skip this leaf node if it has stopped training
			{
				(pActiveLFN->DATA.LFN.pSplit)->nCount++; // we have to zero this before this routine is called.
				(pActiveLFN->DATA.LFN.pSplit)->DBLNOISEVARIANCE += adblX[nDim - 1]; //this is the value of a noise variance sample
			}
		} // end loop over VARfile
		free(adblX);
} // END of spliterrorsetVAR
