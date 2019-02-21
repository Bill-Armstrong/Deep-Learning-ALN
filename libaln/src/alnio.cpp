// ALN Library
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

// alnio.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


static int ALNAPI DoALNWrite(FILE* pFile, const ALN* pALN);
static int ALNAPI DoALNRead(FILE* pFile, ALN** ppALN);


// saving ALN to disk file
// returns ALN_* error code, (ALN_NOERROR on success)
ALNIMP int ALNAPI ALNWrite(const ALN* pALN, const char* pszFileName)
{
   // parameter variance
  if (pALN == NULL)
    return ALN_GENERIC;

  if (pszFileName == NULL)
    return ALN_GENERIC;

  // open a file -- binary mode
  FILE* pFile;
  if(fopen_s(&pFile, pszFileName, "wb") != 0)
    return ALN_ERRFILE;

  int nRet = DoALNWrite(pFile, pALN);

  fclose(pFile);  // will not reset errno

  if (nRet != ALN_NOERROR)
  {
    ASSERT(nRet == ALN_ERRFILE);
    int nErr = errno;     // save it
    remove(pszFileName);
    errno = nErr;
  }

  return nRet;
}

// loading ALN from disk file
// pointer to loaded ALN returned in ppALN
// returns ALN_* error code, (ALN_NOERROR on success)
ALNIMP int ALNAPI ALNRead(const char* pszFileName, ALN** ppALN)
{
  // parameter variance
  if (ppALN == NULL)
    return ALN_GENERIC;

  if (pszFileName == NULL)
    return ALN_GENERIC;

  // open a file -- binary mode
  FILE* pFile;
  if(fopen_s(&pFile, pszFileName, "rb") != 0 )
    return ALN_ERRFILE;

  int nRet = DoALNRead(pFile, ppALN);

  fclose(pFile);  // will not reset errno

  return nRet;  
}

/////////////////////////////////////////////////////////////////////////////
// NOTE: byte order is a problem and is machine dependent

// helpers
static int ALNAPI ReadRegion(FILE* pFile, ALN* pALN, ALNREGION* pRegion);
static int ALNAPI WriteRegion(FILE* pFile, const ALNREGION* pRegion);
static int ALNAPI ReadTree(FILE* pFile, ALN* pALN, ALNNODE* pNode);
static int ALNAPI WriteTree(FILE* pFile, const ALN* pALN, const ALNNODE* pNode);

// following macros won't work if you pass a pointer to be written... 
// ie, treat param n as a reference to var being written
#define _WRITE(f, n) ((int)fwrite(&(n), sizeof((n)), 1, f))
#define _READ(f, n) ((int)fread(&(n), sizeof((n)), 1, f))

#define ALNHDR "ALN35"
#define ALNHDRSIZE 5

static int ALNAPI DoALNWrite(FILE* pFile, const ALN* pALN)
{
  ASSERT(pFile);
  ASSERT(pALN);
  
  // header
  if (fwrite(ALNHDR, ALNHDRSIZE, 1, pFile) != 1) return ALN_ERRFILE;
  
  // write version, dim, output, number of regions
  if (_WRITE(pFile, pALN->nVersion) != 1) return ALN_ERRFILE;
  if (_WRITE(pFile, pALN->nDim) != 1) return ALN_ERRFILE;
  if (_WRITE(pFile, pALN->nOutput) != 1) return ALN_ERRFILE;
  if (_WRITE(pFile, pALN->nRegions) != 1) return ALN_ERRFILE;

  // write region array
  for (int i = 0; i < pALN->nRegions; i++)
  {
    int nRet = WriteRegion(pFile, &(pALN->aRegions[i]));
    if (nRet != ALN_NOERROR)
      return nRet;
  }

  // write tree
  return WriteTree(pFile, pALN, pALN->pTree);
}

static int ALNAPI DoALNRead(FILE* pFile, ALN** ppALN)
{
  ASSERT(pFile);
  ASSERT(ppALN);

  *ppALN = NULL;
 
  // header
  char szHdr[ALNHDRSIZE];
  if (fread(szHdr, ALNHDRSIZE, 1, pFile) != 1) return ALN_ERRFILE;
  if (strncmp(szHdr, ALNHDR, ALNHDRSIZE) != 0)
    return ALN_BADFILEFORMAT;

  // alloc aln
  ALN* pALN = (ALN*)malloc(sizeof(ALN));
  if (pALN == NULL) return ALN_OUTOFMEM;
  memset(pALN, 0, sizeof(ALN));
  
  // read version, dim, output, number of regions,
  if (_READ(pFile, pALN->nVersion) != 1) return ALN_ERRFILE;
  if (_READ(pFile, pALN->nDim) != 1) return ALN_ERRFILE;
  if (_READ(pFile, pALN->nOutput) != 1) return ALN_ERRFILE;
  if (_READ(pFile, pALN->nRegions) != 1) return ALN_ERRFILE;

  if (pALN->nVersion > ALNVER || pALN->nDim < 0 ||
      pALN->nOutput < 0 || pALN->nOutput >= pALN->nDim ||
      pALN->nRegions <= 0)
  {
    pALN->nRegions = 0;
    ALNDestroyALN(pALN);
    return ALN_BADFILEFORMAT;
  }

  // skip over confidence params if version is 0x00030007; no other version has these
  if (pALN->nVersion == 0x00030007)
  {
    // P, lower bound, upper bound, Tail05, Interval4P
    double dbl;
    if (_READ(pFile, dbl) != 1) return ALN_ERRFILE;
    if (_READ(pFile, dbl) != 1) return ALN_ERRFILE;
    if (_READ(pFile, dbl) != 1) return ALN_ERRFILE;
    if (_READ(pFile, dbl) != 1) return ALN_ERRFILE;
    if (_READ(pFile, dbl) != 1) return ALN_ERRFILE;
  }
  // else confidence info already zeroed out, or not present
  
  // allocate region array
  pALN->aRegions = (ALNREGION*)malloc(pALN->nRegions * sizeof(ALNREGION));
  if (pALN->aRegions == NULL)
  {
    pALN->nRegions = 0;
    ALNDestroyALN(pALN);
    return ALN_OUTOFMEM;
  }
  memset(pALN->aRegions, 0, pALN->nRegions * sizeof(ALNREGION));

  // read region array
  for (int i = 0; i < pALN->nRegions; i++)
  {
    int nRet = ReadRegion(pFile, pALN, &(pALN->aRegions[i]));
    if (nRet != ALN_NOERROR)
    {
      ALNDestroyALN(pALN);
      return nRet;
    }
  }

  // allocate first node
  pALN->pTree = (ALNNODE*)malloc(sizeof(ALNNODE));
  if (pALN->pTree == NULL)
  {
    ALNDestroyALN(pALN);
    return ALN_OUTOFMEM;
  }

  int nRet = ReadTree(pFile, pALN, pALN->pTree);
  if (nRet != ALN_NOERROR)
  {
    ALNDestroyALN(pALN);
    return nRet;
  }

  // set new version number
  pALN->nVersion = ALNVER;

  // assign to output param
  *ppALN = pALN;

  return ALN_NOERROR;
}

static int ALNAPI WriteRegion(FILE* pFile, const ALNREGION* pRegion)
{
  ASSERT(pFile);
  ASSERT(pRegion);

  // write parent, learn factor
  if (_WRITE(pFile, pRegion->nParentRegion) != 1) return ALN_ERRFILE;
  if (_WRITE(pFile, pRegion->dblLearnFactor) != 1) return ALN_ERRFILE;

  // don't write var map!

  // write number of constraints
  if (_WRITE(pFile, pRegion->nConstr) != 1) return ALN_ERRFILE;

  // write each constraint...
  for(int i = 0; i < pRegion->nConstr; i++)
  {
    if (_WRITE(pFile, pRegion->aConstr[i]) != 1) return ALN_ERRFILE;
  }

  return ALN_NOERROR;
}

static int ALNAPI ReadRegion(FILE* pFile, ALN* pALN, ALNREGION* pRegion)
{
  ASSERT(pFile);
  ASSERT(pALN);
  ASSERT(pRegion);

  memset(pRegion, 0, sizeof(ALNREGION));

  // read parent, learn factor
  if (_READ(pFile, pRegion->nParentRegion) != 1) return ALN_ERRFILE;
  if (_READ(pFile, pRegion->dblLearnFactor) != 1) return ALN_ERRFILE;
  if (pRegion->nParentRegion < -1 || pRegion->nParentRegion >= pALN->nRegions ||
      pRegion->dblLearnFactor < 0)
    return ALN_BADFILEFORMAT;

  // read number of constraints
  if (_READ(pFile, pRegion->nConstr) != 1) return ALN_ERRFILE;
  if (pRegion->nConstr < 0 || pRegion->nConstr > pALN->nDim)
    return ALN_BADFILEFORMAT;

  // alloc var map?
  if (pRegion->nConstr > 0 && pRegion->nConstr < pALN->nDim)
  {
    pRegion->afVarMap = (char*)malloc(MAPBYTECOUNT(pALN->nDim));
      // failure tolerable, since we can always search for constraint
    if (pRegion->afVarMap)
      memset(pRegion->afVarMap, 0, MAPBYTECOUNT(pALN->nDim));
  }
  else pRegion->afVarMap = NULL;

  // alloc constraints
  pRegion->aConstr = (ALNCONSTRAINT*)malloc(pRegion->nConstr * sizeof(ALNCONSTRAINT));
  if (pRegion->aConstr == NULL)
  {
    pRegion->nConstr = 0;
    return ALN_OUTOFMEM;
  }

  // read each constraint...
  for(int i = 0; i < pRegion->nConstr; i++)
  {
    if (_READ(pFile, pRegion->aConstr[i]) != 1) return ALN_ERRFILE;
    if (pRegion->aConstr[i].nVarIndex < 0 || 
        pRegion->aConstr[i].nVarIndex >= pALN->nDim)
      return ALN_BADFILEFORMAT;
  
    if (pRegion->afVarMap)
    {
      if (TESTMAP(pRegion->afVarMap, pRegion->aConstr[i].nVarIndex))
        return ALN_BADFILEFORMAT; // dup var index!

      SETMAP(pRegion->afVarMap, pRegion->aConstr[i].nVarIndex);
    }
  }

  return ALN_NOERROR;
}


static int ALNAPI WriteTree(FILE* pFile, const ALN* pALN, const ALNNODE* pNode)
{
  ASSERT(pFile);
  ASSERT(pNode);

  // write parent region, node flags
  if (_WRITE(pFile, pNode->nParentRegion) != 1) return ALN_ERRFILE;

  // do not write eval flag!  
  int fNode = pNode->fNode & ~NF_EVAL;  
  if (_WRITE(pFile, fNode) != 1) return ALN_ERRFILE;

  if (pNode->fNode & NF_LFN)
  {
    // vector dim
    if (_WRITE(pFile, LFN_VDIM(pNode)) != 1) return ALN_ERRFILE;

    // var map
    if (LFN_VARMAP(pNode))
    {
      char c = 1;
      if (_WRITE(pFile, c) != 1) return ALN_ERRFILE;
      if ((int)fwrite(LFN_VARMAP(pNode), sizeof(char), 
                      MAPBYTECOUNT(pALN->nDim), pFile) 
            != MAPBYTECOUNT(pALN->nDim)) return ALN_ERRFILE;
    }
    else
    {
      char c = 0;
      if (_WRITE(pFile, c) != 1) return ALN_ERRFILE;
    }

    // split
    if (pNode->fNode & LF_SPLIT)
    {
      ASSERT(LFN_CANSPLIT(pNode));
      if (_WRITE(pFile, LFN_SPLIT_COUNT(pNode)) != 1) return ALN_ERRFILE;
      if (_WRITE(pFile, LFN_SPLIT_SQERR(pNode)) != 1) return ALN_ERRFILE;
      if (_WRITE(pFile, LFN_SPLIT_RESPTOTAL(pNode)) != 1) return ALN_ERRFILE;
      if (_WRITE(pFile, LFN_SPLIT_T(pNode)) != 1) return ALN_ERRFILE;  // Added in Version 0x00030009
    }

    // vectors
    if ((int)fwrite(LFN_W(pNode), sizeof(double), LFN_VDIM(pNode) + 1, pFile)
          != (LFN_VDIM(pNode) + 1)) return ALN_ERRFILE;

    if ((int)fwrite(LFN_C(pNode), sizeof(double), LFN_VDIM(pNode), pFile)
          != LFN_VDIM(pNode)) return ALN_ERRFILE;

    if ((int)fwrite(LFN_D(pNode), sizeof(double), LFN_VDIM(pNode), pFile)
          != LFN_VDIM(pNode)) return ALN_ERRFILE;

	}
  else
  {
    ASSERT(pNode->fNode & NF_MINMAX);

    // number of children
    int nChildren = 2;
    if (_WRITE(pFile, nChildren) != 1) return ALN_ERRFILE;

    for(int i = 0; i < MINMAX_NUMCHILDREN(pNode); i++)
    {
      int nRet = WriteTree(pFile, pALN, MINMAX_CHILDREN(pNode)[i]);
      if (nRet != ALN_NOERROR) return nRet;
    }
  }

  return ALN_NOERROR;
}

static int ALNAPI ReadTree(FILE* pFile, ALN* pALN, ALNNODE* pNode)
{
  ASSERT(pFile);
  ASSERT(pALN);
  ASSERT(pNode);

  memset(pNode, 0, sizeof(ALNNODE));

  // read parent region, node flags
  if (_READ(pFile, pNode->nParentRegion) != 1) return ALN_ERRFILE;
  if (_READ(pFile, pNode->fNode) != 1) return ALN_ERRFILE;
  if (pNode->nParentRegion < 0 || pNode->nParentRegion >= pALN->nRegions ||
      (pNode->fNode & (NF_MINMAX | NF_LFN)) == 0)
    return ALN_BADFILEFORMAT;

  if (pNode->fNode & NF_LFN)
  {
    // vector dim
    if (_READ(pFile, LFN_VDIM(pNode)) != 1) return ALN_ERRFILE;
    if (LFN_VDIM(pNode) < 0 || LFN_VDIM(pNode) > pALN->nDim)
      return ALN_BADFILEFORMAT;

    // var map
    char c;
    if (_READ(pFile, c) != 1) return ALN_ERRFILE;
    if (c != 0) // var map exists
    {
      LFN_VARMAP(pNode) = (char*)malloc(MAPBYTECOUNT(pALN->nDim));
      if (LFN_VARMAP(pNode) == NULL)
        return ALN_OUTOFMEM;

      if ((int)fread(LFN_VARMAP(pNode), sizeof(char), 
                     MAPBYTECOUNT(pALN->nDim), pFile) 
            != MAPBYTECOUNT(pALN->nDim)) return ALN_ERRFILE;
    }

    // split
    if (pNode->fNode & LF_SPLIT)
    {
      LFN_SPLIT(pNode) = (ALNLFNSPLIT*)malloc(sizeof(ALNLFNSPLIT));
      if (LFN_SPLIT(pNode) == NULL)
        return ALN_OUTOFMEM;

      if (_READ(pFile, LFN_SPLIT_COUNT(pNode)) != 1) return ALN_ERRFILE;
      if (_READ(pFile, LFN_SPLIT_SQERR(pNode)) != 1) return ALN_ERRFILE;
      if (_READ(pFile, LFN_SPLIT_RESPTOTAL(pNode)) != 1) return ALN_ERRFILE;

      if (pALN->nVersion >= 0x00030009)
      {
        if (_READ(pFile, LFN_SPLIT_T(pNode)) != 1) return ALN_ERRFILE;
      }
      else
        LFN_SPLIT_T(pNode) = 0;
    }

    // alloc and read vectors
    LFN_W(pNode) = (double*)malloc((LFN_VDIM(pNode) + 1) * sizeof(double));
    if (LFN_W(pNode) == NULL)
      return ALN_OUTOFMEM;

    if ((int)fread(LFN_W(pNode), sizeof(double), LFN_VDIM(pNode) + 1, pFile)
          != (LFN_VDIM(pNode) + 1)) return ALN_ERRFILE;

    LFN_C(pNode) = (double*)malloc(LFN_VDIM(pNode) * sizeof(double));
    if (LFN_C(pNode) == NULL)
      return ALN_OUTOFMEM;

    if ((int)fread(LFN_C(pNode), sizeof(double), LFN_VDIM(pNode), pFile)
          != LFN_VDIM(pNode)) return ALN_ERRFILE;

    LFN_D(pNode) = (double*)malloc(LFN_VDIM(pNode) * sizeof(double));
    if (LFN_D(pNode) == NULL)
      return ALN_OUTOFMEM;

    if ((int)fread(LFN_D(pNode), sizeof(double), LFN_VDIM(pNode), pFile)
          != LFN_VDIM(pNode)) return ALN_ERRFILE;
  }
  else
  {
    ASSERT(pNode->fNode & NF_MINMAX);

    // number of children
    int nChildren;
    if (_READ(pFile, nChildren) != 1) return ALN_ERRFILE;
    if (nChildren != 2)
      return ALN_BADFILEFORMAT;

    // init child ptr array
    memset(MINMAX_CHILDREN(pNode), 0, 2 * sizeof(ALNNODE*));

    // read children
    for(int i = 0; i < 2; i++)
    {
      MINMAX_CHILDREN(pNode)[i] = (ALNNODE*)malloc(sizeof(ALNNODE));
      if (MINMAX_CHILDREN(pNode)[i] == NULL)
        return ALN_OUTOFMEM;

      // set parent
      NODE_PARENT(MINMAX_CHILDREN(pNode)[i]) = pNode;

      int nRet = ReadTree(pFile, pALN, MINMAX_CHILDREN(pNode)[i]);
      if (nRet != ALN_NOERROR) return nRet;
    }
  }

  return ALN_NOERROR;
}

