/*
// dtr_mem.c
// DTREE memory management
// Copyright (C) 1995 Dendronic Decisions Limited
*/

#ifdef DTREEDLL
#define DTRIMP __declspec(dllexport)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>

#include <dtree.h>
#include "dtr_priv.h"

DTRIMP DTREE* DTREEAPI CreateDtree()
{
  DTREE* pDtree = NULL;

  if ((pDtree = (DTREE*)malloc(sizeof(DTREE))) == NULL)
    return NULL;

  memset (pDtree, 0, sizeof(DTREE));
  
  return pDtree;
}

DTRIMP void DTREEAPI DestroyDtree(DTREE* pDtree)
{       
  if (pDtree == NULL)
    return;
  
  /* delete arrays */
  DestroyVarDefArray(pDtree->aVarDefs, pDtree->nDim);
  pDtree->nDim = 0; pDtree->nOutputIndex = 0; pDtree->aVarDefs = 0;
  DestroyLinearFormArray(pDtree->aLinearForms, pDtree->nLinearForms);
  pDtree->nLinearForms = 0; pDtree->aLinearForms = NULL;
  DestroyBlockArray(pDtree->aBlocks, pDtree->nBlocks);
  pDtree->nBlocks = 0; pDtree->aBlocks = NULL;
  DestroyDtreeNodeArray(pDtree->aNodes);
  pDtree->nNodes = 0; pDtree->aNodes = NULL;
  free(pDtree);
}

DTRIMP VARDEF* DTREEAPI CreateVarDefArray(int nDim)
{                                              
  VARDEF* aVarDefs = NULL;

  if ((aVarDefs = (VARDEF*)malloc(nDim * sizeof(VARDEF))) == NULL)
    return NULL;
                  
  memset(aVarDefs, 0, nDim * sizeof(VARDEF));
  
  return aVarDefs;
}

DTRIMP void DTREEAPI DestroyVarDefArray(VARDEF* aVarDefs, int nDim)
{
  int i;

  if (aVarDefs == NULL)
    return;
               
  /* delete names */
  for (i = 0; i < nDim; i++)
  {
    if (aVarDefs[i].pszName != NULL)
      free(aVarDefs[i].pszName);
  }
               
  free(aVarDefs);
}

DTRIMP char* DTREEAPI SetVarDefName(VARDEF* pVarDef, const char* pszName)
{
  if (pVarDef == NULL)
    return NULL;
      
  /* delete existing name */
  if (pVarDef->pszName != NULL)   
  {
    free(pVarDef->pszName);
    pVarDef->pszName = NULL;
  }
                            
  if (pszName == NULL)
    return NULL;

  /* copy name */
  pVarDef->pszName = (char*)malloc((strlen(pszName) + 1) * sizeof(char));
  strcpy(pVarDef->pszName, pszName);
  
  return pVarDef->pszName;
}

DTRIMP DTREENODE* DTREEAPI CreateDtreeNodeArray(int nNodes)
{                                              
  DTREENODE* aNodes = NULL;
  
  if ((aNodes = (DTREENODE*)malloc(nNodes * sizeof(DTREENODE))) == NULL)
    return NULL;
                  
  memset(aNodes, 0, nNodes * sizeof(DTREENODE));
  
  return aNodes;
}

DTRIMP void DTREEAPI DestroyDtreeNodeArray(DTREENODE* aNodes)
{
  if (aNodes == NULL)
    return;
               
  free(aNodes);
}

DTRIMP BLOCK* DTREEAPI CreateBlockArray(int nBlocks)
{                                              
  BLOCK* aBlocks = NULL;

  if ((aBlocks = (BLOCK*)malloc(nBlocks * sizeof(BLOCK))) == NULL)
    return NULL;
                  
  memset(aBlocks, 0, nBlocks * sizeof(BLOCK));
  
  return aBlocks;
}

DTRIMP void DTREEAPI DestroyBlockArray(BLOCK* aBlocks, int nBlocks)
{              
  int i;

  if (aBlocks == NULL)
    return; 
    
  /* delete min/max trees */
  for (i = 0; i < nBlocks; i++)
  {
    DestroyMinMaxNode(aBlocks[i].pMinMaxTree);
  }
               
  free(aBlocks);
}

DTRIMP MINMAXNODE* DTREEAPI CreateMinMaxNode()
{
  MINMAXNODE* pMinMaxNode = NULL;
  
  if ((pMinMaxNode = (MINMAXNODE*)malloc(sizeof(MINMAXNODE))) == NULL)
    return NULL;
    
  memset(pMinMaxNode, 0, sizeof(MINMAXNODE));
  
  return pMinMaxNode;
}

DTRIMP MINMAXNODE* DTREEAPI CopyMinMaxNode(MINMAXNODE* pMMN)
{
  MINMAXNODE* pCopy = NULL;    
  
  if (pMMN == NULL)
    return NULL;
    
  if ((pCopy = (MINMAXNODE*)malloc(sizeof(MINMAXNODE))) == NULL)
    return NULL;
  
  memset(pCopy, 0, sizeof(MINMAXNODE));
   
  /* copy */
  pCopy->nType = pMMN->nType;
  if (pCopy->nType == DTREE_LINEAR)  /* leaf node */
  {
    MMN_LFINDEX(pCopy) = MMN_LFINDEX(pMMN);
  }
  else /* internal node */
  {
    /* traverse and copy child list */
    MINMAXNODE* pList = MMN_CHILDLIST(pMMN);
    MINMAXNODE* pListCopy = MMN_CHILDLIST(pCopy);
    while(pList != NULL)
    {
      MINMAXNODE* pChild = CopyMinMaxNode(pList);
      if (pChild == NULL)
      {
        DestroyMinMaxNode(pCopy);  
        return NULL;
      }
      if (pListCopy == NULL)
      {
        MMN_CHILDLIST(pCopy) = pChild;
        pListCopy = pChild;
      }
      else
      {
        pListCopy->pNext = pChild;
        pListCopy = pChild;
      }
      
      pList = pList->pNext;
    }
  }
  
  return pCopy;
}
                      
DTRIMP MINMAXNODE* DTREEAPI AddMinMaxNodeChild(MINMAXNODE* pParent)
{ 
  MINMAXNODE* pChild = NULL;
  
  if (pParent == NULL)
    return NULL;
  
  /* pParent must be a min or max op to add children */
  if (pParent->nType == DTREE_LINEAR)
    return NULL;
    
  if ((pChild = CreateMinMaxNode()) == NULL)
    return NULL;
      
  /* append to END of child list */
  if (MMN_CHILDLIST(pParent) == NULL)
    MMN_CHILDLIST(pParent) = pChild;
  else
  {
    MINMAXNODE* pList = MMN_CHILDLIST(pParent);
    while(pList->pNext != NULL)
      pList = pList->pNext;
    pList->pNext = pChild;  /* fixup sibling pointer */
  }

  return pChild;  
}
                      
DTRIMP void DTREEAPI DestroyMinMaxNode(MINMAXNODE* pMinMaxNode)
{
  if (pMinMaxNode == NULL)
    return;
    
  if (pMinMaxNode->nType != DTREE_LINEAR)
  {
    /* delete child list */
    MINMAXNODE* pList = MMN_CHILDLIST(pMinMaxNode);
    while(pList != NULL)
    {
      MINMAXNODE* pListNext = pList->pNext;
      DestroyMinMaxNode(pList);
      pList = pListNext;
    }
  }
  
  free(pMinMaxNode);
}

DTRIMP LINEARFORM* DTREEAPI CreateLinearFormArray(int nForms, int nDim)
{                                              
  int i;
  size_t sizeVec = sizeof(double) * (size_t)nDim;  
  
  LINEARFORM* aForms = NULL;
  if ((aForms = (LINEARFORM*)malloc(nForms * sizeof(LINEARFORM))) == NULL)
    return NULL;
                  
  memset(aForms, 0, nForms * sizeof(LINEARFORM));
  
  /* alloc weight and centroid vectors */
  for (i = 0; i < nForms; i++)
  {                                  
    /* weights */
    if ((aForms[i].adblW = (double*)malloc(sizeVec)) == NULL)
    {
      DestroyLinearFormArray(aForms, nForms);
      return NULL;
    }
    memset(aForms[i].adblW, 0, sizeVec);
    
    /* centroid */
    if ((aForms[i].adblC = (double*)malloc(sizeVec)) == NULL)
    {
      DestroyLinearFormArray(aForms, nForms);
      return NULL;
    }
    memset(aForms[i].adblC, 0, sizeVec);
  }
  
  return aForms;
}

DTRIMP void DTREEAPI DestroyLinearFormArray(LINEARFORM* aForms, int nForms)
{              
  int i;

  if (aForms == NULL)
    return; 
    
  /* delete weight and centroid vectors */
  for (i = 0; i < nForms; i++)
  {
    free(aForms[i].adblW);
    free(aForms[i].adblC);
  }
               
  free(aForms);
}
