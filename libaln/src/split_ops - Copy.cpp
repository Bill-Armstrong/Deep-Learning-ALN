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
BOOL bStopTraining;


void dosplitcontrol(CMyAln* pALN, ALNNODE* pNode, double dblLimit);
void dodivideTR(CMyAln* pALN, ALNNODE* pNode);
void dodivideVAR(CMyAln* pALN, ALNNODE* pNode);
void spliterrorsetTR(CMyAln * pALN);
void spliterrorsetVAR(CMyAln * pALN);
void dozerospliterror(CMyAln* pALN, ALNNODE* pNode);

// the following tells us how to use the SPLIT typedef between trainings of an ALN
#define DBLNOISEVARIANCE dblRespTotal

void splitcontrol(CMyAln* pALN, double dblLimit)  // routine
{
  ASSERT(pALN);
  ASSERT(pALN->GetTree());
	// initialize the SPLIT components to zero
	dozerospliterror(pALN, pALN->GetTree());
	// get square errors of pieces on training set
	spliterrorsetTR(pALN);
	// divide the training errors of the pieces by the hit counts, set counts to zero
	dodivideTR(pALN,pALN->GetTree());
	// get the values in VARfile which estimate the local noise variance
	if(bEstimateRMSError && (dblLimit > 1.0)) // dblLimit is set to 0 for overtraining
	{
		// there is a set of noise variance samples
		spliterrorsetVAR(pALN);
		// divide the sum of local noise estimates on each piece by its count of hits
		dodivideVAR(pALN,pALN->GetTree());
	}
	// this conducts a search for all leaves in the ALN
  dosplitcontrol(pALN, pALN->GetTree(), dblLimit);
  // reset the SPLIT components to zero
	dozerospliterror(pALN, pALN->GetTree());
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


void dozerospliterror(CMyAln* pALN, ALNNODE* pNode) // routine
{
	// initializes counters of all leaf nodes before the next training period
  ASSERT(pNode);
  if (NODE_ISMINMAX(pNode))
  {
    dozerospliterror(pALN, MINMAX_LEFT(pNode));
    dozerospliterror(pALN, MINMAX_RIGHT(pNode));
  }
  else
  {
		ASSERT(NODE_ISLFN(pNode));
		(pNode->DATA.LFN.pSplit)->nCount = 0;
    (pNode->DATA.LFN.pSplit)->dblSqError = 0;
    (pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE = 0;		
	}
}

void dodivideTR(CMyAln* pALN, ALNNODE* pNode) // routine
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
		if (!LFN_ISINIT(pNode)) // skip this leaf node if it has stopped training
		{
			if ((pNode->DATA.LFN.pSplit)->nCount > nDim)  // the piece is overdetermined (if points are distict) 
			{
				(pNode->DATA.LFN.pSplit)->dblSqError /=
					(pNode->DATA.LFN.pSplit)->nCount;
				(pNode->DATA.LFN.pSplit)->nCount = 0;  // inserted WWA 2009.10.06 IMPORTANT ERROR WAS HERE
				// this is zeroed in order to use nCount for spliterrorsetVAR

			}
			else
			{
				// if the count of training samples hitting the piece is nDim or less, the piece is not overdetermined;
				// We don't want to split it so we make the training error of the piece zero
				(pNode->DATA.LFN.pSplit)->dblSqError = 0;
			}
		}
	}
}


void dodivideVAR(CMyAln* pALN, ALNNODE* pNode) // routine
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
		if (!LFN_ISINIT(pNode)) // skip this leaf node if it has stopped training
		{
			if((pNode->DATA.LFN.pSplit)->nCount > nDim) // we need at more than nDim samples to do stats
			{
				(pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE /=
											 (pNode->DATA.LFN.pSplit)->nCount;
			}
			else
			{
				// if there are too few or no variance samples, make sure we don't divide
				(pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE = dblMax;
			}
		}
	}
}



void dosplitcontrol(CMyAln* pALN, ALNNODE* pNode, double dblLimit) // routine
{
	// this routine visits all the leaf nodes and determines whether the
	// training error is below dblLimit times the noise variance
	// During overtraining the dblLimit is less than 1.0. It is set low, which causes lots of splitting
	double dblSqErrorPieceTrain;
	double dblPieceNoiseVariance;
	ASSERT(pNode);
	if (NODE_ISMINMAX(pNode))
	{
		dosplitcontrol(pALN, MINMAX_LEFT(pNode), dblLimit);
		dosplitcontrol(pALN, MINMAX_RIGHT(pNode), dblLimit);
	}
	else
	{
		ASSERT(NODE_ISLFN(pNode));
		if (NODE_ISCONSTANT(pNode))
		{
			return;                     // no adaptation of this constant subtree
		}
		dblSqErrorPieceTrain = (pNode->DATA.LFN.pSplit)->dblSqError; // average square error on TRfile
		if (bEstimateRMSError)
		{
			// if we are evaluating noise variance (which we always do in training until we store a noise variance ALN)
			dblPieceNoiseVariance = (pNode->DATA.LFN.pSplit)->DBLNOISEVARIANCE; 
		}
		else
		{
			dblPieceNoiseVariance = 1.0; // this allows use of dblLimit to stop splitting for regression and overtraining
		}
		if (dblSqErrorPieceTrain < dblPieceNoiseVariance * dblLimit) // this implements the F-test criterion for stopping training
		{
			// since this piece fits well enough, given the noise variance, training of it stops
			// stop all future splitting of this leaf node (LFN)
			LFN_FLAGS(pNode) &= ~LF_SPLIT;  // ~ is the unary one's complement operator
			LFN_FLAGS(pNode) &= NF_CONSTANT; // Make the node constant
		}
		else
		{
			if (bStopTraining == FALSE) bStopTraining = TRUE; // if every leaf has fitted to within noise, training stops
		}
	}
}

void spliterrorsetTR(CMyAln * pALN) // routine
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
		adblX[nDim - 1] = 0; // not used in evaluation by QuickEval
		predict = pALN->QuickEval(adblX, &pActiveLFN);
		if (!LFN_ISINIT(pActiveLFN)) // skip this leaf node if it has stopped training
		{
			se = (predict - desired) * (predict - desired);
			(pActiveLFN->DATA.LFN.pSplit)->nCount++;
			(pActiveLFN->DATA.LFN.pSplit)->dblSqError += se;
		}
	} // end loop over TRset
	free(adblX);
} // END of spliterrorsetTR


void spliterrorsetVAR(CMyAln * pALN) // routine
{

		// assign the square errors on the variance set to the leaf nodes of the ALN
		double * adblX = (double *)malloc((nDim) * sizeof(double));
		double desired = 0;
		double value = 0;
		ALNNODE* pActiveLFN;
		double      se = 0; // square error added to LFN DBLNOISEVARIANCE
		for (int j = 0; j < nRowsVAR; j++)   // this is expensive using the whole variance set, but more accurate
		{
			for (int i = 0; i < nDim; i++)
			{
				adblX[i] = VARfile.GetAt(j, i, 0); // the value at nDim - 1 is used only for desired
			}
			// pAln has to be the current approximant! Is this correct?
			value = pALN->QuickEval(adblX, &pActiveLFN);  // this is the piece of the current approximant that the X-vector lies on.
			if (!LFN_ISINIT(pActiveLFN)) // skip this leaf node if it has stopped training
			{
				//se = (predict - desired) * (predict - desired);
				// now correct for the dimension. the average variance of predict - desired is !+ 2/(nDim+2), so we have to divide se by this
				(pActiveLFN->DATA.LFN.pSplit)->nCount++;
				(pActiveLFN->DATA.LFN.pSplit)->DBLNOISEVARIANCE += adblX[nDim - 1]; //this is the value of a noise variance sample
			}
		} // end loop over VARfile
		free(adblX);
} // END of spliterrorsetVAR
