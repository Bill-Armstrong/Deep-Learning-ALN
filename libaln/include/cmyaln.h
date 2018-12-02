//file cmyaln.h

#ifndef _CMYALN_H
#define _CMYALN_H

#include <stdio.h>

extern "C" int _kbhit();
extern "C" int _getch();
extern FILE* fpProtocol;
extern double dblTrainErr;
extern int nMaxEpochs;
extern int nNumberLFNs;

class CMyAln;
void fillvector(double *, CMyAln *);

class CMyAln : public CAln
{
  // notification overrides
  public:
  virtual BOOL OnTrainStart(TRAININFO* pTrainInfo, void* pvData)
  { 
    //cerr << "Training starts..." << endl; 
		//fprintf(fpProtocol,"Training starts ... \n");
    return TRUE; 
  }
  
	virtual BOOL OnTrainEnd(TRAININFO* pTrainInfo, void* pvData) 
  { 
    //cerr << "Training finished.  RMSE: " << pTrainInfo->dblRMSErr << endl; 
		//fprintf(fpProtocol,"Training finished.  Training set RMSE = %f \n", pTrainInfo->dblRMSErr);
		dblTrainErr = pTrainInfo->dblRMSErr;
		return TRUE;
  }

  virtual BOOL OnEpochStart(EPOCHINFO* pEpochInfo, void* pvData) 
  { 
    //cerr << "Epoch " << pEpochInfo->nEpoch << ":" << flush;
		//fprintf(fpProtocol,"Epoch %d:\n",pEpochInfo->nEpoch);
    return TRUE; 
  }

  virtual BOOL OnEpochEnd(EPOCHINFO* pEpochInfo, void* pvData) 
  {
		if(pEpochInfo->nEpoch == (nMaxEpochs -1))
		{
      nNumberLFNs = pEpochInfo->nActiveLFNs;
		  fprintf(fpProtocol,"Estimated RMSE %f Active/Total LFNs %d/%d\n", pEpochInfo->dblEstRMSErr,
			        pEpochInfo->nActiveLFNs, pEpochInfo->nLFNs);
		}
	  return TRUE;
  }

  virtual BOOL OnVectorInfo(VECTORINFO* pVectorInfo, void* pvData) 
  {
 		fillvector(	pVectorInfo->adblX, this);
    if(_kbhit())
    {
      //cerr << endl << "Checking keyboard input (type 's' to stop training)" << endl;
      return _getch()!='s';
    }
    else
    {
      return TRUE;
    }
  }
};

#endif
