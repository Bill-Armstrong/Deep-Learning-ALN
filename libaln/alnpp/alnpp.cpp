// ALN Library C++ wrapper
// ALNfit Learning Engine for approximation of functions defined by samples.
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

// alnpp.cpp

#ifdef __GNUC__
#include <typeinfo>
#endif

#include <stdlib.h>
#include <memory.h>
#include <alnpp.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#ifndef ASSERT
#define ASSERT ALNASSERT
#endif

/////////////////////////////////////////////////////////////////////////////
// class CAln

CAln::CAln()
{
  m_pALN = NULL;
  memset(&m_datainfo, 0, sizeof(m_datainfo));
  m_nLastError = ALN_GENERIC; // no ALN pointer yet!
}

CAln::~CAln()
{
  Destroy();
  ASSERT(m_pALN == NULL);
}

ALNREGION* CAln::GetRegion(int nRegion)
{
  if (m_pALN == NULL)
    return NULL;

  ASSERT(nRegion >= 0 && nRegion < m_pALN->nRegions);
  return m_pALN->aRegions + nRegion;
}

const ALNREGION* CAln::GetRegion(int nRegion) const
{
  if (m_pALN == NULL)
    return NULL;

  ASSERT(nRegion >= 0 && nRegion < m_pALN->nRegions);
  return m_pALN->aRegions + nRegion;
}

ALNCONSTRAINT* CAln::GetConstraint(int nVar, int nRegion)
{
  if (m_pALN == NULL)
    return NULL;

  ASSERT(nRegion >= 0 && nRegion < m_pALN->nRegions);
  ASSERT(nVar >= 0 && nVar < m_pALN->nDim);
  return m_pALN->aRegions[nRegion].aConstr + nVar;
}

const ALNCONSTRAINT* CAln::GetConstraint(int nVar, int nRegion) const
{
  if (m_pALN == NULL)
    return NULL;

  ASSERT(nRegion >= 0 && nRegion < m_pALN->nRegions);
  ASSERT(nVar >= 0 && nVar < m_pALN->nDim);
  return m_pALN->aRegions[nRegion].aConstr + nVar;
}

double CAln::GetEpsilon(int nVar, int nRegion) const
{
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->dblEpsilon; 
}

void CAln::SetEpsilon(double dblEpsilon, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  pConstr->dblEpsilon = dblEpsilon;
	pConstr->dblSqEpsilon = dblEpsilon * dblEpsilon;
}

double CAln::GetWeightMin(int nVar, int nRegion) const
{ 
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->dblWMin; 
}

void CAln::SetWeightMin(double dblWMin, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);
  ASSERT(nVar != m_pALN->nOutput);  
    // can't change output var weight!
  
  if (nVar != m_pALN->nOutput)
    pConstr->dblWMin = dblWMin; 
}

double CAln::GetWeightMax(int nVar, int nRegion) const
{ 
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->dblWMax; 
}

void CAln::SetWeightMax(double dblWMax, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);
  ASSERT(nVar != m_pALN->nOutput);  
    // can't change output var weight!
  
  if (nVar != m_pALN->nOutput)
    pConstr->dblWMax = dblWMax; 
}  

double CAln::GetMin(int nVar, int nRegion) const
{ 
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->dblMin; 
}

void CAln::SetMin(double dblMin, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);
  
  pConstr->dblMin = dblMin; 
}

double CAln::GetMax(int nVar, int nRegion) const
{ 
  const ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);

  return pConstr->dblMax; 
}

void CAln::SetMax(double dblMax, int nVar, int nRegion)
{ 
  ALNCONSTRAINT* pConstr = GetConstraint(nVar, nRegion);
  ASSERT(pConstr != NULL);
  
  pConstr->dblMax = dblMax; 
}  

ALNNODE* CAln::GetTree()
{
  if (m_pALN == NULL)
    return NULL;

  return m_pALN->pTree;
}

const ALNNODE* CAln::GetTree() const
{
  if (m_pALN == NULL)
    return NULL;

  return m_pALN->pTree;
}

void CAln::SetDataInfo(int nPoints, int nCols, const double* adblData,
                       const VARINFO* aVarInfo /*= NULL*/)
{
  m_datainfo.nCols = nCols;
  m_datainfo.nPoints = nPoints;
  m_datainfo.adblData = adblData;
  m_datainfo.aVarInfo = aVarInfo;
}

BOOL CAln::Create(int nDim, int nOutput)
{
  Destroy();
  ASSERT(m_pALN == NULL);

  m_pALN = ALNCreateALN(nDim, nOutput);
  if (m_pALN != NULL)
  {
    m_nLastError = ALN_NOERROR;
  }
  else
  {
    m_nLastError = ALN_GENERIC;
  }

  return m_pALN != NULL;
}

#ifdef ENABLE_REGIONS
int CAln::AddRegion(int nParentRegion, double dblLearnFactor, 
                    int nConstr, int* anConstr)
{
  m_nLastError = ALN_NOERROR;
  return ALNAddRegion(m_pALN, nParentRegion, dblLearnFactor,
                      nConstr, anConstr);
}
#endif

BOOL CAln::AddLFNs(ALNNODE* pParent, int nParentMinMaxType, 
                   int nLFNs, ALNNODE** apLFNs /*= NULL*/)
{
  m_nLastError = ALNAddLFNs(m_pALN, pParent, nParentMinMaxType, nLFNs, apLFNs);
  return m_nLastError == ALN_NOERROR;
}

BOOL CAln::AddMultiLayer(ALNNODE* pParent, int nParentMinMaxType, 
                         int nLayers, int nMinMaxFanin, int nLFNFanin, 
                         int nFlags)
{
  m_nLastError = ALNAddMultiLayer(m_pALN, pParent, nParentMinMaxType, nLayers, 
                                  nMinMaxFanin, nLFNFanin, nFlags);
  return m_nLastError == ALN_NOERROR;
}

BOOL CAln::AddTreeString(ALNNODE* pParent, const char* pszTreeString, 
                         int& nParsed)
{
  m_nLastError = ALNAddTreeString(m_pALN, pParent, pszTreeString, &nParsed);
  return m_nLastError == ALN_NOERROR;
}

BOOL CAln::SetGrowable(ALNNODE* pNode)
{
  m_nLastError = ALN_NOERROR;
  return ALNSetGrowable(m_pALN, pNode);
}

BOOL CAln::Destroy()
{
  m_nLastError = ALN_NOERROR;
  BOOL b = ALNDestroyALN(m_pALN);
  if (b) m_pALN = NULL;
  return b;
}

// private callback data struct
struct CALLBACKDATA
{
  CAln* pALN;
  void* pvData;
};

// training - uses internal ALNDATAINFO if pData == NULL
BOOL CAln::Train(int nMaxEpochs, double dblMinRMSErr, double dblLearnRate,
                 BOOL bJitter, int nNotifyMask /*= AN_NONE*/, 
                 ALNDATAINFO* pData /*= NULL*/, void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

  CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;

  m_nLastError = ALNTrain(m_pALN, pData, &callback, nMaxEpochs, dblMinRMSErr, dblLearnRate, bJitter);

  return (m_nLastError == ALN_NOERROR || m_nLastError == ALN_USERABORT);
}

double CAln::CalcRMSError(int nNotifyMask /*= AN_NONE*/, 
                          ALNDATAINFO* pData /*= NULL*/, 
                          void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

  CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;
  
  double dbl;
  m_nLastError = ::ALNCalcRMSError(m_pALN, pData, &callback, &dbl);

  if (m_nLastError != ALN_NOERROR)
    dbl = -1.0;

  return dbl;
}

// eval
BOOL CAln::Eval(double* adblResult, int* pnStart /*=NULL*/, 
                int* pnEnd /*= NULL*/, int nNotifyMask /*= AN_NONE*/,
                ALNDATAINFO* pData /*= NULL*/, void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

  CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;

  m_nLastError = ALNEval(m_pALN, pData, &callback, adblResult, pnStart, pnEnd);

  return m_nLastError == ALN_NOERROR;
}

// quick eval
double CAln::QuickEval(const double* adblX, ALNNODE** ppActiveLFN /*= NULL*/)
{
  m_nLastError = ALN_NOERROR;
  return ALNQuickEval(m_pALN, adblX, ppActiveLFN);
}

// get variable monotonicicty, returns -1 on failure
int CAln::VarMono(int nVar)
{
  int nMono;
  m_nLastError = ALNVarMono(m_pALN, nVar, &nMono);
  return nMono;
}

// invert the aln to get an aln for a different output variable WWA
BOOL CAln::Invert(int nVar)
{
	m_nLastError = ALNInvert(m_pALN, nVar);
	return (m_nLastError == ALN_NOERROR);
}

// save ALN to disk file
BOOL CAln::Write(const char* pszFileName)
{
  m_nLastError = ALNWrite(m_pALN, pszFileName);
  return m_nLastError == ALN_NOERROR;
}

// read ALN from disk file... destroys any existing ALN
BOOL CAln::Read(const char* pszFileName)
{
  Destroy();
  ASSERT(m_pALN == NULL);

  m_nLastError = ALNRead(pszFileName, &m_pALN);
  return m_nLastError == ALN_NOERROR;
}

// conversion to dtree
DTREE* CAln::ConvertDtree(int nMaxDepth)
{
  DTREE* pDtree = NULL;
	DTREE** ppDtree = &pDtree;
  m_nLastError = ALNConvertDtree(m_pALN, nMaxDepth, ppDtree);
	pDtree = *ppDtree;
  return pDtree;
}

// confidence intervals
BOOL CAln::CalcConfidence(ALNCONFIDENCE* pConfidence, 
                          int nNotifyMask /*= AN_NONE*/, 
                          ALNDATAINFO* pData /*= NULL*/, 
                          void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

  CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;

  m_nLastError = ALNCalcConfidence(m_pALN, pData, &callback, pConfidence);

  return m_nLastError == ALN_NOERROR;
}

double CAln::ConfidencePLimit(const ALNCONFIDENCE* pConfidence, 
                              double dblSignificance)
{
  double dblPLimit;
  
  if (ALNConfidencePLimit(pConfidence, dblSignificance, &dblPLimit) != ALN_NOERROR)
    return -1.0;
  
  return dblPLimit;
}

double CAln::ConfidenceTLimit(const ALNCONFIDENCE* pConfidence, 
                              double dblInterval)
{
  double dblTLimit;
  
  if (ALNConfidenceTLimit(pConfidence, dblInterval, &dblTLimit) != ALN_NOERROR)
    return -1.0;
  
  return dblTLimit;
}

// lfn analysis
BOOL CAln::LFNAnalysis(void*& pvAnalysis,
                       int& nLFNStats,
                       int nNotifyMask /*= AN_NONE*/, 
                       ALNDATAINFO* pData /*= NULL*/, void* pvData /*= NULL*/)
{
  if (pData == NULL)
    pData = &m_datainfo;

	CALLBACKDATA data;
  data.pALN = this;
  data.pvData = pvData;

  ALNCALLBACKINFO callback;
  callback.nNotifyMask = nNotifyMask;
  callback.pvData = &data;
  callback.pfnNotifyProc = ALNNotifyProc;

  m_nLastError = ALNLFNAnalysis(m_pALN, pData, &callback,
                                pvAnalysis, nLFNStats);

  return m_nLastError == ALN_NOERROR;
}

BOOL CAln::LFNFreeAnalysis(void* pvAnalysis)
{
  return ALNLFNFreeAnalysis(pvAnalysis) == ALN_NOERROR;
}

BOOL CAln::LFNStats(void* pvAnalysis, 
                    int nLFNStat,
                    LFNSTATS& LFNStats,
                    int& nWeights,
                    ALNNODE*& pLFN)
{
  return ALNLFNStats(pvAnalysis, nLFNStat, &LFNStats, 
                     &nWeights, &pLFN) == ALN_NOERROR;
}

BOOL CAln::LFNWeightStats(void* pvAnalysis,
                          int nLFNStat,
                          int nWeightStat,
                          LFNWEIGHTSTATS& LFNWeightStats)
{
  return ALNLFNWeightStats(pvAnalysis, nLFNStat, nWeightStat, 
                           &LFNWeightStats) == ALN_NOERROR;
}

int ALNAPI CAln::ALNNotifyProc(const ALN* pALN, int nCode, void* pParam, 
                               void* pvData)
{
  CALLBACKDATA* pData = (CALLBACKDATA*)pvData;
  CAln* pALNObj = pData->pALN;
  
  ASSERT(pALNObj != NULL && pALN == pALNObj->m_pALN);

  BOOL bContinue = FALSE;

  switch (nCode)
	{
		case AN_TRAINSTART:
			bContinue = pALNObj->OnTrainStart((TRAININFO*)pParam, pData->pvData);
			break;

		case AN_TRAINEND:
		  bContinue = pALNObj->OnTrainEnd((TRAININFO*)pParam, pData->pvData);
			break;

		case AN_EPOCHSTART:
		  bContinue = pALNObj->OnEpochStart((EPOCHINFO*)pParam, pData->pvData);
			break;

		case AN_EPOCHEND:
      bContinue = pALNObj->OnEpochEnd((EPOCHINFO*)pParam, pData->pvData);
			break;

		case AN_ADAPTSTART:
		  bContinue = pALNObj->OnAdaptStart((ADAPTINFO*)pParam, pData->pvData);
      break;

		case AN_ADAPTEND:
		  bContinue = pALNObj->OnAdaptEnd((ADAPTINFO*)pParam, pData->pvData);
      break;

    case AN_LFNADAPTSTART:
      bContinue = pALNObj->OnLFNAdaptStart((LFNADAPTINFO*)pParam, pData->pvData);
      break;
    
    case AN_LFNADAPTEND:
      bContinue = pALNObj->OnLFNAdaptEnd((LFNADAPTINFO*)pParam, pData->pvData);
      break;
    
    case AN_VECTORINFO:
      bContinue = pALNObj->OnVectorInfo((VECTORINFO*)pParam, pData->pvData);
      break;
	}

	return bContinue;
}
