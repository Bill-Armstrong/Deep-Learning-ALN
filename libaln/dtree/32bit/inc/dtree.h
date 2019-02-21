/*
// dtree.h
// DTREE definitions and prototypes - 32 bit version
// Copyright (C) 1995 - 2001 Dendronic Decisions Limited
*/

#ifndef __DTREE_H__
#define __DTREE_H__

#ifdef __cplusplus
extern "C" {
#endif             

/* define __cdecl for non-Microsoft compilers */
#ifndef NODEFINE_CDECL

#if	( !defined(_MSC_VER) && !defined(__cdecl) )
#define __cdecl
#endif

#endif  /* NODEFINE_CDECL */

/* define __stdcall for non-Microsoft compilers */
#ifndef NODEFINE_STDCALL

#if	( !defined(_MSC_VER) && !defined(__stdcall) )
#define __stdcall
#endif

#endif /* NODEFINE_STDCALL */

/*
/////////////////////////////////////////////////////////////////////
// If you are using DTR1032(d).DLL in your application or DLL,
// then you should define the preprocessor symbol DTREEDLL.
*/

#ifdef DTREEDLL
#ifndef DTRIMP
#define DTRIMP __declspec(dllimport)
#endif
#else
#define DTRIMP
#endif

#define DTREEAPI __stdcall
                  
/*                  
/////////////////////////////////////////////////////////////////////
// DTREE version info
*/
#define DTREE_VERMAJOR 0x0001
#define DTREE_VERMINOR 0x0004
                                  
/* returns major version in upper 16 bits, minor in lower 16 bits */
DTRIMP int DTREEAPI GetDtreeVersion();  

/*            
/////////////////////////////////////////////////////////////////////
// Decision Tree Structures
*/

#if defined(_MSC_VER)
#pragma pack(4)                   /* structures aligned to DWORD boundaries */
#endif

/* DTREE min max node types */
#define DTREE_MIN      0          
#define DTREE_MAX      1
#define DTREE_LINEAR   2

typedef struct tagLINEARFORM      /* describes a linear form                */
{                    
  double dblBias;                 /* bias: computed from centroid,weights   */
  double* adblW;                  /* weight vector has nDim elements        */
  double* adblC;                  /* centroid vector has nDim elements      */
} LINEARFORM;                     /* 16 bytes                               */

typedef struct tagMINMAXNODE      /* describes a min/max tree node          */
{
  int nType;                      /* node type                              */
  union
  {       
    int nLFIndex;                 /* linearform index if nType == LINEAR    */
    struct tagMINMAXNODE* pChildList;  /* linked list of children           */
                                  /*   if nType == MIN or nType == MAX      */
  } info;
  struct tagMINMAXNODE* pNext;    /* pointer to sibling                     */
} MINMAXNODE;                     /* 12 bytes                               */

#define MMN_LFINDEX(pMMN) ((pMMN)->info.nLFIndex)
#define MMN_CHILDLIST(pMMN) ((pMMN)->info.pChildList)
 
typedef struct tagDTREENODE       /* describes a decision tree node         */
{ 
  int nLeaf;                      /* non-zero if this is a leaf node        */
  int nParentIndex;               /* index of parent (-1 if no parent)      */
  union
  {
    int nBlockIndex;              /* index of block if this is a leaf       */
    struct
    {                             /* if this is an internal node...         */
      double dblT;                /* threshold                              */
      int nVarIndex;              /* variable index                         */
      int nLeftIndex;             /* index of left child                    */
      int nRightIndex;            /* index of right child                   */
    } node;
  } info;                                         
} DTREENODE;                      /* 28 bytes                               */

#define DNODE_BLOCKINDEX(pNode) ((pNode)->info.nBlockIndex)
#define DNODE_THRESHOLD(pNode) ((pNode)->info.node.dblT)
#define DNODE_VARINDEX(pNode) ((pNode)->info.node.nVarIndex)
#define DNODE_LEFTINDEX(pNode) ((pNode)->info.node.nLeftIndex)
#define DNODE_RIGHTINDEX(pNode) ((pNode)->info.node.nRightIndex)

typedef struct tagBLOCK           /* describes a block structure            */
{
  MINMAXNODE* pMinMaxTree;        /* min/max tree                           */
  int nDtreeIndex;                /* dtree leaf index that refs this        */
} BLOCK;                          /* 8 bytes                                */ 
  
typedef struct tagVARBOUND        /* variable bound structure               */
{
  double dblMin;                  /* variable minimum                       */
  double dblMax;                  /* variable maximum                       */
} VARBOUND;                       /* 16 bytes                               */

typedef struct tagVARDEF          /* variable definition                    */
{          
  VARBOUND bound;                 /* variable bound                         */
  char* pszName;                  /* variable name                          */
} VARDEF;                         /* 20 bytes                               */

#define VARDEF_MIN(pVar) ((pVar)->bound.dblMin)
#define VARDEF_MAX(pVar) ((pVar)->bound.dblMax)
  
typedef struct tagDTREE           /* main decision tree struct              */
{
  int nDim;                       /* dimension of space                     */
  int nOutputIndex;               /* output var index                       */
  VARDEF* aVarDefs;               /* array of variable definitions          */
  int nLinearForms;               /* number linear forms                    */
  LINEARFORM* aLinearForms;       /* array of linear forms                  */
  int nBlocks;                    /* number of blocks                       */
  BLOCK* aBlocks;                 /* array of blocks                        */
  int nNodes;                     /* number of dtree nodes                  */
  DTREENODE* aNodes;              /* array of nodes                         */
} DTREE;                          /* 36 bytes                               */

#if defined(_MSC_VER)
#pragma pack()                    /* restore default structure alignment    */
#endif
                                                                     
/*                                                                     
/////////////////////////////////////////////////////////////////////
// DTREE memory management prototypes
*/

DTRIMP DTREE* DTREEAPI CreateDtree();
DTRIMP void DTREEAPI DestroyDtree(DTREE* pDtree);
          
DTRIMP VARDEF* DTREEAPI CreateVarDefArray(int nDim);          
DTRIMP void DTREEAPI DestroyVarDefArray(VARDEF* aVarDefs, int nDim);
DTRIMP char* DTREEAPI SetVarDefName(VARDEF* pVarDef, const char* pszName);

DTRIMP DTREENODE* DTREEAPI CreateDtreeNodeArray(int nNodes);
DTRIMP void DTREEAPI DestroyDtreeNodeArray(DTREENODE* aNodes);

DTRIMP BLOCK* DTREEAPI CreateBlockArray(int nBlocks);
DTRIMP void DTREEAPI DestroyBlockArray(BLOCK* aBlocks, int nBlocks);

DTRIMP MINMAXNODE* DTREEAPI CreateMinMaxNode();
DTRIMP MINMAXNODE* DTREEAPI CopyMinMaxNode(MINMAXNODE* pMinMaxNode);
DTRIMP MINMAXNODE* DTREEAPI AddMinMaxNodeChild(MINMAXNODE* pParent);
DTRIMP void DTREEAPI DestroyMinMaxNode(MINMAXNODE* pMinMaxNode);

DTRIMP LINEARFORM* DTREEAPI CreateLinearFormArray(int nForms, int nDim);
DTRIMP void DTREEAPI DestroyLinearFormArray(LINEARFORM* aLinearForms, int nForms);

/*
/////////////////////////////////////////////////////////////////////
// Dtree I/O routines
*/

/* global dtree line number counter */
#ifdef _WIN32
DTRIMP int* __cdecl _dtree_lineno(void);
#define dtree_lineno (*_dtree_lineno())
#else 
extern int dtree_lineno;
#endif
                    
/* returns DTR_NOERROR on success 
   places a new DTREE in *ppDtree - use DestroyDtree to Destroy *ppDtree */
DTRIMP int DTREEAPI ReadDtree(const char* pszFileName, DTREE** ppDtree);
DTRIMP int DTREEAPI WriteDtree(const char* pszFileName, DTREE* pDtree);

/*
/////////////////////////////////////////////////////////////////////
// Dtree binary I/O routines
*/
                   
/* returns DTR_NOERROR on success 
   places a new DTREE in *ppDtree - use DestroyDtree to Destroy *ppDtree */
DTRIMP int DTREEAPI BinReadDtree(const char* pszFileName, DTREE** ppDtree);
DTRIMP int DTREEAPI BinWriteDtree(const char* pszFileName, DTREE* pDtree);

/*                   
/////////////////////////////////////////////////////////////////////
// Dtree evaluation routines                   
*/
                                                                        
/* returns DTR_NOERROR on success             
   places a result in *pdblResult
   returns index of linear form that calculated the result in
   plLinearIndex (if not NULL) */
DTRIMP int DTREEAPI EvalDtree(DTREE* pDtree, double* adblInput, 
                              double* pdblResult, int* pnLinearIndex);

/* min/max tree evaluation         
   returns DTR_NOERROR on success         
   places result in pdblResult, 
   places responsible linear index in plLinearIndex (if not NULL) */
DTRIMP int DTREEAPI EvalMinMaxTree(MINMAXNODE* pMMN, LINEARFORM* aLF, int nDim, 
                                   int nOutput, double* adblInput, 
                                   double* pdblResult, int* pnLinearIndex);
                            
/* linear form evaluation         
   returns DTR_NOERROROR on success         
   places result in pdblResult */
DTRIMP int DTREEAPI EvalLinearForm(LINEARFORM* pLF, int nDim, int nOutput, 
                                   double* adblInput, double* pdblResult);

                                 
/*                                 
/////////////////////////////////////////////////////////////////////
// Dtree error handling
*/

/* global dtree error number */
#ifdef _WIN32
DTRIMP int* __cdecl _dtree_errno(void);
#define dtree_errno (*_dtree_errno())
#else
extern int dtree_errno;
#endif

/* retrieves string associated with nErrNo */
DTRIMP void DTREEAPI GetDtreeError(int nErrno, char* pBuf, int nMaxBufLen);

/* error codes */
#define DTR_NOERROR                 0

#define DTR_ERRORBASE               10000
                                 
#define DTR_GENERIC                 (DTR_ERRORBASE + 10)
#define DTR_FILEERR                 (DTR_ERRORBASE + 30)
#define DTR_FILEWRITEERR            (DTR_ERRORBASE + 31)
#define DTR_FILEREADERR             (DTR_ERRORBASE + 32)
#define DTR_FILETYPEERR             (DTR_ERRORBASE + 33)
#define DTR_ENDIANERR               (DTR_ERRORBASE + 34)
#define DTR_MALLOCFAILED            (DTR_ERRORBASE + 50)

#define DTR_BADVERSIONDEF           (DTR_ERRORBASE + 100)
#define DTR_BADVERSIONINT           (DTR_ERRORBASE + 101)
#define DTR_MISSINGVERSIONSEMI      (DTR_ERRORBASE + 102)
#define DTR_UNKNOWNVERSION          (DTR_ERRORBASE + 103)
                                     
#define DTR_BADDIMDEF               (DTR_ERRORBASE + 150)
#define DTR_BADDIMINT               (DTR_ERRORBASE + 151)
#define DTR_BADDIMRANGE             (DTR_ERRORBASE + 152)
#define DTR_MISSINGDIMSEMI          (DTR_ERRORBASE + 153)

#define DTR_BADVARDEF               (DTR_ERRORBASE + 200)
#define DTR_BADVARDEFIDENT          (DTR_ERRORBASE + 201) 
#define DTR_MISSINGVARDEFCOLON      (DTR_ERRORBASE + 202)
#define DTR_MISSINGVARBOUNDSTART    (DTR_ERRORBASE + 203)
#define DTR_MISSINGVARBOUNDCOMMA    (DTR_ERRORBASE + 204)
#define DTR_BADVARBOUND             (DTR_ERRORBASE + 205)
#define DTR_MISSINGVARBOUNDEND      (DTR_ERRORBASE + 206)
#define DTR_NEGATIVEVARBOUNDRANGE   (DTR_ERRORBASE + 207)
#define DTR_MISSINGVARDEFSEMI       (DTR_ERRORBASE + 208)
                                     
#define DTR_BADOUTPUTDEF            (DTR_ERRORBASE + 300)
#define DTR_BADOUTPUTIDENT          (DTR_ERRORBASE + 301)
#define DTR_BADOUTPUTRANGE          (DTR_ERRORBASE + 302)
#define DTR_ZEROOUTPUTWEIGHT        (DTR_ERRORBASE + 303)
#define DTR_MISSINGOUTPUTSEMI       (DTR_ERRORBASE + 304)

#define DTR_UNKNOWNVARIDENT         (DTR_ERRORBASE + 350)
#define DTR_UNDEFINEDREGION         (DTR_ERRORBASE + 351)
#define DTR_UNDEFINEDLINEARFORM     (DTR_ERRORBASE + 352)

#define DTR_BADBLOCKDEF             (DTR_ERRORBASE + 400)
#define DTR_BADBLOCKINT             (DTR_ERRORBASE + 401)
#define DTR_MISSINGBLOCKINDEX       (DTR_ERRORBASE + 402)
#define DTR_DUPBLOCKINDEX           (DTR_ERRORBASE + 403)
#define DTR_MISSINGBLOCKCOLON       (DTR_ERRORBASE + 404)  
#define DTR_BADMINMAXNODE           (DTR_ERRORBASE + 405) 
#define DTR_MISSINGMINMAXLISTSTART  (DTR_ERRORBASE + 406)
#define DTR_EMPTYMINMAXLIST         (DTR_ERRORBASE + 407)
#define DTR_BADMINMAXLISTELEM       (DTR_ERRORBASE + 408)  
#define DTR_BADBLOCKINDEXRANGE      (DTR_ERRORBASE + 409)
#define DTR_MISSINGBLOCKSEMI        (DTR_ERRORBASE + 410)

#define DTR_BADLINEARDEF            (DTR_ERRORBASE + 500)
#define DTR_BADLINEARINT            (DTR_ERRORBASE + 501)
#define DTR_MISSINGLINEARINDEX      (DTR_ERRORBASE + 502)
#define DTR_DUPLINEARINDEX          (DTR_ERRORBASE + 503)
#define DTR_MISSINGLINEARCOLON      (DTR_ERRORBASE + 504)
#define DTR_BADLINEARWEIGHT         (DTR_ERRORBASE + 505)
#define DTR_MISSINGCENTROIDSTART    (DTR_ERRORBASE + 506)
#define DTR_MISSINGCENTROIDEND      (DTR_ERRORBASE + 507)
#define DTR_MISSINGCENTROIDIDENT    (DTR_ERRORBASE + 508)
#define DTR_MISSINGLINEARSIGN       (DTR_ERRORBASE + 509)
#define DTR_DUPVARCENTROID          (DTR_ERRORBASE + 510)
#define DTR_BADLINEARCENTROID       (DTR_ERRORBASE + 511)
#define DTR_BADLINEARINDEXRANGE     (DTR_ERRORBASE + 512)
#define DTR_MISSINGLINEARSEMI       (DTR_ERRORBASE + 513)
                                     
#define DTR_BADDTREEDEF             (DTR_ERRORBASE + 600)                                     
#define DTR_BADDTREEINT             (DTR_ERRORBASE + 601)
#define DTR_MISSINGDTREEINDEX       (DTR_ERRORBASE + 602)
#define DTR_DUPDTREEINDEX           (DTR_ERRORBASE + 603)
#define DTR_MISSINGDTREECOLON       (DTR_ERRORBASE + 604)
#define DTR_MISSINGDTREEVAR         (DTR_ERRORBASE + 605)
#define DTR_MISSINGDTREELE          (DTR_ERRORBASE + 606)
#define DTR_MISSINGDTREETHRESHOLD   (DTR_ERRORBASE + 607)
#define DTR_BADTHRESHOLDRANGE       (DTR_ERRORBASE + 608)
#define DTR_MISSINGDTREEPAREN       (DTR_ERRORBASE + 609)
#define DTR_MISSINGDTREEQUESTION    (DTR_ERRORBASE + 610)
#define DTR_MISSINGDTREECHILDINDEX  (DTR_ERRORBASE + 611)
#define DTR_MISSINGDTREECHILDSEP    (DTR_ERRORBASE + 612)
#define DTR_CIRCULARDTREECHILDREF   (DTR_ERRORBASE + 613)
#define DTR_BADDTREEINDEXRANGE      (DTR_ERRORBASE + 614)
#define DTR_MISSINGDTREEBLOCKINDEX  (DTR_ERRORBASE + 615)
#define DTR_MISSINGDTREESEMI        (DTR_ERRORBASE + 616)
                   
#ifdef __cplusplus
}
#endif

#endif /* __DTREE_H__ */
