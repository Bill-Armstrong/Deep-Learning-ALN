/*
// dtree.c
// DTREE implementation
// Copyright (C) 1995 Dendronic Decisions Limited
*/

#ifdef DTREEDLL
#define DTRIMP __declspec(dllexport)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <dtree.h>
#include "dtr_priv.h" 

#ifdef _WIN32
#include <windows.h>  // for dynamic TLS routines
#endif

DTRIMP int DTREEAPI GetDtreeVersion()
{                              
  int nVer = ((int)DTREE_VERMINOR) | ((int)DTREE_VERMAJOR << 16);
  return nVer;
}

DTRIMP int DTREEAPI ReadDtree(const char* pszFileName, DTREE** ppDtree)
{      
  FILE* pFile = NULL;
  int nErr = DTR_NOERROR;

  /* set line number to 1 */
  dtree_lineno = 1;

  /* open file */
  if ((pFile = fopen(pszFileName, "r")) == NULL)
  {
    return DTR_FILEERR;
  }           
  
  /* parse the file */
  nErr = ParseDtreeFile(pFile, ppDtree);
  
  /* set the error number */
  dtree_errno = nErr;

  /* close file */
  fclose(pFile);
  return nErr;
}

DTRIMP int DTREEAPI WriteDtree(const char* pszFileName, DTREE* pDtree)
{      
  FILE* pFile = NULL;
  int nErr = DTR_NOERROR;
  
  /* open file */
  if ((pFile = fopen(pszFileName, "w")) == NULL)
  {
    return DTR_FILEERR;
  }           
  
  /* export the file */
  nErr = ExportDtreeFile(pFile, pszFileName, pDtree);
  
  /* set the error number */
  dtree_errno = nErr;

  /* close file */
  fclose(pFile);
  return nErr;
}

DTRIMP int DTREEAPI BinReadDtree(const char* pszFileName, DTREE** ppDtree)
{
  FILE* pFile = NULL;
  int nErr = DTR_NOERROR;

  /* set line number to 1 */
  dtree_lineno = 1;

  /* open file */
  if ((pFile = fopen(pszFileName, "rb")) == NULL)
  {
    return DTR_FILEERR;
  }           
  
  /* parse the file */
  nErr = ReadBinDtreeFile(pFile, ppDtree);
  
  /* set the error number */
  dtree_errno = nErr;

  /* close file */
  fclose(pFile);
  return nErr;
}

DTRIMP int DTREEAPI BinWriteDtree(const char* pszFileName, DTREE* pDtree)
{
  FILE* pFile = NULL;
  int nErr = DTR_NOERROR;
  
  /* open file */
  if ((pFile = fopen(pszFileName, "wb")) == NULL)
  {
    return DTR_FILEERR;
  }           
  
  /* export the file */
  nErr = WriteBinDtreeFile(pFile, pDtree);
  
  /* set the error number */
  dtree_errno = nErr;

  /* close file */
  fclose(pFile);
  return nErr;
}

DTRIMP void DTREEAPI GetDtreeError(int nErrno, char* pBuf, int nMaxBufLen)
{
  GetErrMsg(nErrno, pBuf, nMaxBufLen);
}

DTRIMP int DTREEAPI EvalDtree(DTREE* pDtree, double* adblInput, 
                              double* pdblResult, int* pnLinearIndex)
{                    
  DTREENODE* pNode = NULL;
  int nErr = DTR_NOERROR;
  
  /* find leaf */
  pNode = pDtree->aNodes; 
  while (pNode->nLeaf == 0)
  {
    if (adblInput[DNODE_VARINDEX(pNode)] <= DNODE_THRESHOLD(pNode))
      pNode = pDtree->aNodes + DNODE_LEFTINDEX(pNode);
    else
      pNode = pDtree->aNodes + DNODE_RIGHTINDEX(pNode);
  }
  
  /* eval block */
  nErr = EvalMinMaxTree(pDtree->aBlocks[DNODE_BLOCKINDEX(pNode)].pMinMaxTree,
                        pDtree->aLinearForms, pDtree->nDim, 
                        pDtree->nOutputIndex, 
                        adblInput, pdblResult, pnLinearIndex);
                        
  /* bound output */
  if (*pdblResult < pDtree->aVarDefs[pDtree->nOutputIndex].bound.dblMin)
    *pdblResult = pDtree->aVarDefs[pDtree->nOutputIndex].bound.dblMin;
  else if (*pdblResult > pDtree->aVarDefs[pDtree->nOutputIndex].bound.dblMax)
    *pdblResult = pDtree->aVarDefs[pDtree->nOutputIndex].bound.dblMax;
  
  return nErr;
}
    
DTRIMP int DTREEAPI EvalMinMaxTree(MINMAXNODE* pMMN, LINEARFORM* aLF, 
                                   int nDim, int nOutput, double* adblInput, 
                                   double* pdblResult, int* pnLinearIndex)
{
  if (pMMN->nType == DTREE_LINEAR)
  { 
    int nErr = DTR_NOERROR;
    if ((nErr = EvalLinearForm(aLF + MMN_LFINDEX(pMMN), nDim, nOutput, 
                               adblInput, pdblResult)) != DTR_NOERROR)
      return nErr;
    
    if (pnLinearIndex != NULL)
      *pnLinearIndex = MMN_LFINDEX(pMMN);
    return DTR_NOERROR;
  }
  else if (pMMN->nType == DTREE_MIN || pMMN->nType == DTREE_MAX)
  {   
    /* eval children */
    MINMAXNODE* pList = MMN_CHILDLIST(pMMN);
    while (pList != NULL)
    {
      int lIndex;
      double dbl;
      
      int nErr;
      if ((nErr = EvalMinMaxTree(pList, aLF, nDim, nOutput, 
                                 adblInput, &dbl, &lIndex)) != DTR_NOERROR)
        return nErr;  
      
      if ((pList == MMN_CHILDLIST(pMMN)) ||                     /* first time thru */
          (pMMN->nType == DTREE_MIN && dbl < *pdblResult) ||    /* or this is a min */
          (pMMN->nType == DTREE_MAX && dbl > *pdblResult))      /* or this is a max */
      {
        *pdblResult = dbl;
        if (pnLinearIndex != NULL)
          *pnLinearIndex = lIndex;
      }
      pList = pList->pNext;
    }       
    
    return DTR_NOERROR;
  }
  
  return DTR_GENERIC; /* unknown node type */
}

DTRIMP int DTREEAPI EvalLinearForm(LINEARFORM* pLF, int nDim, int nOutput,
                                   double* adblInput, double* pdblResult)
{        
  register int i;
  if (pLF->adblW[nOutput] == 0)
    return DTR_ZEROOUTPUTWEIGHT;
  /* eval linear form, xi = (w0 + w1x1 + ... + wi-1xi-1 + wi+1 xi+1 + ... + wnxn) / -wi */
  *pdblResult = pLF->dblBias;
  for (i = 0; i < nDim; i++)
  {                         
    if (i == nOutput) continue;    /* skip output var */
    *pdblResult += adblInput[i] * pLF->adblW[i];
  }
  *pdblResult /= -pLF->adblW[nOutput];
  return DTR_NOERROR;
}

#ifdef _WIN32

typedef struct tagDTREETLS
{
  int _dtree_errno;
  int _dtree_lineno;
} DTREETLS;

// thread local storage for dtree_errno and dtree_lineno in WIN32 version
static DWORD g_dwTlsIndex;

static BOOL DtreeAllocTls()
{
  // Allocate a thread-local storage index.
  g_dwTlsIndex = TlsAlloc();

  if (g_dwTlsIndex == TLS_OUT_OF_INDEXES) 
  {
    // The TLS index couldn't be allocated
    return FALSE;
  }

  return TRUE;
}

static DTREETLS* DtreeAllocGlobals()
{
  DTREETLS* ptls = NULL;
  
  ptls = (DTREETLS*)malloc(sizeof(DTREETLS));
  if (ptls != NULL)
  {
    // init values
    ptls->_dtree_errno = ptls->_dtree_lineno = 0;
    TlsSetValue(g_dwTlsIndex, ptls);
  }
  
  return ptls;
}

static void DtreeFreeGlobals()
{
  DTREETLS* ptls = NULL;

  // Ensure that the TLS index was allocated successfully.
  if (g_dwTlsIndex != TLS_OUT_OF_INDEXES) 
  {
    // Get the pointer to the allocated memory.
    ptls = (DTREETLS *)TlsGetValue(g_dwTlsIndex);

    // Test whether memory was ever allocated for this thread.
    if (ptls != NULL) 
    {
      free(ptls);
    }
  }
}

static void DtreeFreeTls()
{
  // Ensure that the TLS index was allocated successfully.
  if (g_dwTlsIndex != TLS_OUT_OF_INDEXES) 
  {
    // free any allocated globals
    DtreeFreeGlobals();
    
    // Free the TLS index.
    TlsFree(g_dwTlsIndex);
  }    
}

#if defined(_DLL) && !defined(NO_DTREE_DLLMAIN)

// DLL entry point
BOOL WINAPI DllMain (HMODULE hMod, DWORD fdwReason, LPVOID lpvReserved)
{
  switch (fdwReason) 
  {
    case DLL_PROCESS_ATTACH:
      // The calling process is attaching a DLL from its address space.
      return DtreeAllocTls();  
      
    case DLL_THREAD_ATTACH:
      // A new thread is being created in the process.
      break;

    case DLL_THREAD_DETACH:
      // A thread is exiting cleanly.
      DtreeFreeGlobals();
      break;

    case DLL_PROCESS_DETACH:
      // The calling process is detaching the DLL from its address space.
      DtreeFreeTls();
      break;
   }
   return(TRUE);
}

static DTREETLS* DtreeGetTls()
{
  // non-DLL version of DtreeGetTls located below...

  DTREETLS* ptls = NULL;
  
  if (g_dwTlsIndex == TLS_OUT_OF_INDEXES)
    return NULL;

  // get thread specific pointer to globals
  ptls = (DTREETLS*)TlsGetValue(g_dwTlsIndex);
  
  if (ptls == NULL)  // allocate memory for this thread?
    ptls = DtreeAllocGlobals();

  return ptls;
}

#else   // !_DLL || NO_DTREE_DLLMAIN

// fake exit point to de-init TLS
static void __cdecl TlsDeInit(void)
{
  DtreeFreeTls();
}

// fake entry point to init TLS
static BOOL __cdecl TlsInit(void)
{
  // register exit function
  atexit(TlsDeInit);

  return DtreeAllocTls();
}

// global var indicates whether TLS initialized
static BOOL bTlsInit = FALSE;

static DTREETLS* DtreeGetTls()
{
  // DLL version of DtreeGetTls is located above...

  DTREETLS* ptls = NULL;

  if (!bTlsInit)  // TLS init here
  {
    bTlsInit = TlsInit();
    if (!bTlsInit)
      return ptls;
  }
  
  if (g_dwTlsIndex == TLS_OUT_OF_INDEXES)
    return NULL;

  // get thread specific pointer to globals
  ptls = (DTREETLS*)TlsGetValue(g_dwTlsIndex);
  
  if (ptls == NULL)  // allocate memory for this thread?
    ptls = DtreeAllocGlobals();

  return ptls;
}

#endif  // _DLL && !NO_DTREE_DLLMAIN

#endif  // _WIN32

/* global dtree error flag */
#ifdef _WIN32
DTRIMP int* __cdecl _dtree_errno(void)
{
  DTREETLS* ptls = DtreeGetTls();
  
  if (ptls) 
  {
    return &(ptls->_dtree_errno);
  }

  return NULL;
}
#else
int dtree_errno = 0;
#endif

/* global file line number counter */
#ifdef _WIN32
DTRIMP int* __cdecl _dtree_lineno(void)
{
  DTREETLS* ptls = DtreeGetTls();
     
  if (ptls) 
  {
    return &(ptls->_dtree_lineno);
  }
  
  return NULL;
}
#else
int dtree_lineno = 0;
#endif
