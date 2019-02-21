/*
// dtr_io.c
// DTREE I/O routines
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
#include <time.h>

#include <dtree.h>
#include "dtr_priv.h"
  
/*
///////////////////////////////////////////////////////////////////
// lexical & syntactical analysis prototypes                       
*/

/* 
// token types                                                     
*/ 
#define ERR_TOKEN     256
#define IDENTIFIER    257 
#define L_INT       258   /* literals */
#define L_DBL       259
#define R_VERSION     300   /* reserved words */
#define R_OUTPUT      301
#define R_VARDEFS     302
#define R_DTREE       303
#define R_BLOCKS      304
#define R_LINEARFORMS 305
#define R_MIN         306
#define R_MAX         307
#define R_BLOCK       308    
#define S_LE          400   /* symbols */
                          
/*                          
// token structure
*/
#define MAXTOKLEN 128
typedef struct tagTOKEN
{
  long lFilePos;
  union
  {
    long l;
    double dbl;
    char sz[MAXTOKLEN];
  } value;
} TOKEN;        
  
/*
// reserved words
*/
char* _szReserved[] = 
{ 
  "VERSION", "OUTPUT", "VARIABLES",
  "DTREE", "BLOCKS", "LINEARFORMS", 
  "MIN", "MAX", "BLOCK", 
  "", 
};  

int _nReserved[]  = 
{ 
  R_VERSION, R_OUTPUT, R_VARDEFS,
  R_DTREE, R_BLOCKS, R_LINEARFORMS,  
  R_MIN, R_MAX, R_BLOCK,
  0, 
};  
  
/*                         
// parser helpers
*/
int GetVarDef(FILE* pFile, VARDEF* pVarDef);
int GetDtreeNode(FILE* pFile, BLOCK* aBlocks, int nBlocks, 
                 VARDEF* aVarDefs, int nDim, DTREENODE* aNodes, 
                 int nNodes, int nNode);
int GetBlock(FILE* pFile, LINEARFORM* aLinearForms, int nLinearForms, 
             BLOCK* pBlock);
int GetMinMaxTree(FILE* pFile, LINEARFORM* aLinearForms, int nLinearForms,
                  MINMAXNODE** ppMinMaxNode);       
int GetLinearForm(FILE* pFile, int nDim, int nOutput, 
                  LINEARFORM* pLinearForm, VARDEF* aVarDefs);

/* case insensitive string comparison */
#if (defined(_MSC_VER) && !defined(stricmp))
#define stricmp _stricmp     /* Microsoft version */
#endif

/* 
// most UNIX compilers can define the the preprocessor symbol stricmp 
// as strcmpcase in their makefiles...
*/

#ifndef stricmp /* non-ANSI function stricmp must be a macro definition */
#error Please define stricmp to a case-insensitive string comparison routine
#error ... the parameters and return value for stricmp are identical to strcmp.
#endif  /* stricmp */
  
/*
// lexer helpers
*/
void SkipWS(FILE* pFile);                     /* skips white space and comments */
int GetToken(FILE* pFile, TOKEN* pTok);       /* gets next token */
int GetIntegerToken(FILE* pFile, TOKEN* pTok);/* gets an integer token */
int GetDoubleToken(FILE* pFile, TOKEN* pTok); /* gets double token, will convert an integer to double */
void PushbackToken(FILE* pFile, TOKEN* pTok); /* pushes token back onto stream */
int GetSymbol(FILE* pFile, TOKEN* pTok);      /* gets symbol token */
int GetNumber(FILE* pFile, TOKEN* pTok);      /* gets number token */
int GetIdentifier(FILE* pFile, TOKEN* pTok);  /* get identifier token */
int IsReserved(char* psz);                    /* checks if identifier is a reserved word */
int IsComment(FILE* pFile);                   /* eats any comments */

#define ISCSYM(c) ((c) == '_' || isalnum((c)))  /* valid symbol character */
#define ISCSYMF(c) ((c) == '_' || isalpha((c))) /* valid symbol first character */


/*
// export helpers
*/
int WriteMinMaxTree(FILE* pFile, DTREE* pDtree, MINMAXNODE* pMMN, int nIndent);


/* 
// bitmap macros
*/

/* helper macro to calc byte count of bitmap fields for n elements */
/* ... assumes char width is 8 bits */
#define MAPBYTECOUNT(n) (((n) + 7) >> 3)

/* helper macros test, set, and reset bits in map */
/* use n / 8 == n >> 3   and   n % 8 == n & 0x07 */
#define TESTMAP(pMap, n) (((char*)pMap)[(n)>>3] & (0x01<<((n)&0x07)))
#define SETMAP(pMap, n) (((char*)pMap)[(n)>>3] |= (0x01<<((n)&0x07)))
#define RESETMAP(pMap, n) (((char*)pMap)[(n)>>3] &= ~(0x01<<((n)&0x07)))

 
  
/*    
/////////////////////////////////////////////////////////////////////
// main parser
*/
int ParseDtreeFile(FILE* pFile, DTREE** ppDtree)
{
  int nErr = DTR_NOERROR;                     /* error number */
  int i;                                      /* loop counter */
  TOKEN tok;                                  /* token */
  int nVerMajor = 0;                          /* major version */
  int nVerMinor = 0;                          /* minor version */
  int nDim = 0;                               /* dimension */
  int nOutput = 0;                            /* output index */
  VARDEF* aVarDefs = NULL;                    /* array of variable defs */
  int nLinearForms = 0;                       /* number of linear forms */
  LINEARFORM* aLinearForms = NULL;            /* array of linear forms */
  int nBlocks = 0;                            /* number of blocks */
  BLOCK* aBlocks = NULL;                      /* array of blocks */
  int nNodes = 0;                             /* number of nodes */
  DTREENODE* aNodes = NULL;                   /* array of nodes */
  char* aMapSeen = NULL;                      /* map to track parsed indexes */
  
  /* check arguments */
  if (pFile == NULL || ppDtree == NULL)
    return DTR_GENERIC;
  
  /* init dtree pointer */
  *ppDtree = NULL;

  /* parsing helper definitions */
  #define _ON_ERR(nErrCode) { nErr = nErrCode; goto PARSE_ERR; }
  #define _TOKEN(nType, nErrCode) if (GetToken(pFile, &tok) != nType) _ON_ERR(nErrCode)
  #define _INTTOKEN(nErrCode) if (GetIntegerToken(pFile, &tok) != L_INT) _ON_ERR(nErrCode)
  #define _DBLTOKEN(nErrCode) if (GetDoubleToken(pFile, &tok) != L_DBL) _ON_ERR(nErrCode)
  
  /*///////// get version token */
  _TOKEN(R_VERSION, DTR_BADVERSIONDEF);     /* version token */
  _TOKEN('=', DTR_BADVERSIONDEF);           /* assignment */
  _INTTOKEN(DTR_BADVERSIONINT);             /* major version */
  nVerMajor = tok.value.l;
  _TOKEN('.', DTR_BADVERSIONDEF);           /* dot */
  _INTTOKEN(DTR_BADVERSIONINT);             /* minor version */
  nVerMinor = tok.value.l;               
  _TOKEN(';', DTR_MISSINGVERSIONSEMI);      /* end of statement */
  if ((nVerMajor > DTREE_VERMAJOR) || (nVerMinor > DTREE_VERMINOR))
    _ON_ERR(DTR_UNKNOWNVERSION);
                    
  /*//////// get var def list */
  _TOKEN(R_VARDEFS, DTR_BADDIMDEF);         /* vardef token */
  _TOKEN('=', DTR_BADDIMDEF);               /* assignment */
  _TOKEN(L_INT, DTR_BADDIMINT);           /* dimension */
  nDim = tok.value.l;
  _TOKEN(';', DTR_MISSINGDIMSEMI);          /* end of statement */
  if (nDim < 2)
    _ON_ERR(DTR_BADDIMRANGE);
  if ((aVarDefs = CreateVarDefArray(nDim)) == NULL)
    _ON_ERR(DTR_MALLOCFAILED);
  for (i = 0; i < nDim; i++)                /* parse each var */
  {      
    if ((nErr = GetVarDef(pFile, aVarDefs + i)) != DTR_NOERROR)
      _ON_ERR(nErr);
  }
  
  /*///////// get output identifier */
  _TOKEN(R_OUTPUT, DTR_BADOUTPUTDEF);       /* output token */
  _TOKEN('=', DTR_BADOUTPUTDEF);            /* assignment */
  _TOKEN(IDENTIFIER, DTR_BADOUTPUTIDENT);   /* var ident */
  for (i = 0; i < nDim; i++)
  {   
    /* var names are case sensitive */
    if (strcmp(aVarDefs[i].pszName, tok.value.sz) == 0)  
      break;
  }
  _TOKEN(';', DTR_MISSINGOUTPUTSEMI);       /* end of statement */
  if (i == nDim)                        
    _ON_ERR(DTR_UNKNOWNVARIDENT);           /* var ident not found */
  nOutput = i;                          
  if (nOutput < 0 || nOutput >= nDim)
    _ON_ERR(DTR_BADOUTPUTRANGE);            /* out of range */
  
  /*///////// get linear forms */
  _TOKEN(R_LINEARFORMS, DTR_BADLINEARDEF);  /* linearforms token */
  _TOKEN('=', DTR_BADLINEARDEF);            /* assignment */
  _TOKEN(L_INT, DTR_BADLINEARINT);        /* number of linear forms */
  nLinearForms = tok.value.l;
  _TOKEN(';', DTR_MISSINGLINEARSEMI);       /* end of statement */
  if ((aLinearForms = CreateLinearFormArray(nLinearForms, nDim)) == NULL)
    _ON_ERR(DTR_MALLOCFAILED);

  /* allocate index seen map */
  aMapSeen = (char*)malloc(MAPBYTECOUNT(nLinearForms));
  if (aMapSeen == NULL)
    _ON_ERR(DTR_MALLOCFAILED);    
  memset(aMapSeen, 0, MAPBYTECOUNT(nLinearForms));
  for (i = 0; i < nLinearForms; i++)        /* parse each linear form */
  { 
    int nIndex = -1;
    _TOKEN(L_INT, DTR_MISSINGLINEARINDEX);/* linear array index */
    nIndex = tok.value.l;                     
                  
    if (nIndex >= nLinearForms)
      _ON_ERR(DTR_BADLINEARINDEXRANGE);                  
                  
    if (TESTMAP(aMapSeen, nIndex))          /* slot should be empty */
      _ON_ERR(DTR_DUPLINEARINDEX);
    _TOKEN(':', DTR_MISSINGLINEARCOLON);    /* ':' */
    
    /* get form */
    if ((nErr = GetLinearForm(pFile, nDim, nOutput, 
                              aLinearForms + nIndex, aVarDefs)) != DTR_NOERROR)
      _ON_ERR(nErr);
    
    SETMAP(aMapSeen, nIndex);               /* mark as seen */
  } 
  free(aMapSeen);
  aMapSeen = NULL;
  
  /*///////// get blocks */
  _TOKEN(R_BLOCKS, DTR_BADBLOCKDEF);        /* blocks token */
  _TOKEN('=', DTR_BADBLOCKDEF);             /* assignment */
  _TOKEN(L_INT, DTR_BADBLOCKINT);         /* number of blocks */
  nBlocks = tok.value.l;
  _TOKEN(';', DTR_MISSINGBLOCKSEMI);        /* end of statement */
  if ((aBlocks = CreateBlockArray(nBlocks)) == NULL)
    _ON_ERR(DTR_MALLOCFAILED);
  
  /* allocate index seen map */
  aMapSeen = (char*)malloc(MAPBYTECOUNT(nBlocks));
  if (aMapSeen == NULL)
    _ON_ERR(DTR_MALLOCFAILED);    
  memset(aMapSeen, 0, MAPBYTECOUNT(nBlocks));
  for (i = 0; i < nBlocks; i++)             /* for each block */
  {
    int nIndex = -1;
    _TOKEN(L_INT, DTR_MISSINGBLOCKINDEX); /* block array index */
    nIndex = tok.value.l;
      
    if (nIndex >= nBlocks)
      _ON_ERR(DTR_BADBLOCKINDEXRANGE);                        
      
    if (TESTMAP(aMapSeen, nIndex))          /* slot should be empty */
      _ON_ERR(DTR_DUPBLOCKINDEX);
    _TOKEN(':', DTR_MISSINGBLOCKCOLON);     /* ':' */
                                    
    /* get block */
    if ((nErr = GetBlock(pFile, aLinearForms, nLinearForms, 
                         aBlocks + nIndex)) != DTR_NOERROR)
      _ON_ERR(nErr);
    
    SETMAP(aMapSeen, nIndex);               /* mark as seen */
  } 
  free(aMapSeen);
  aMapSeen = NULL;
  
  /*///////// get dtree nodes */
  _TOKEN(R_DTREE, DTR_BADDTREEDEF);         /* dtree token */
  _TOKEN('=', DTR_BADDTREEDEF);             /* assignment */
  _TOKEN(L_INT, DTR_BADDTREEINT);         /* number of nodes */
  nNodes = tok.value.l;
  _TOKEN(';', DTR_MISSINGDTREESEMI);        /* end of statement */
  if ((aNodes = CreateDtreeNodeArray(nNodes)) == NULL)
    _ON_ERR(DTR_MALLOCFAILED);

  /* allocate index seen map */
  aMapSeen = (char*)malloc(MAPBYTECOUNT(nNodes));
  if (aMapSeen == NULL)
    _ON_ERR(DTR_MALLOCFAILED);    
  memset(aMapSeen, 0, MAPBYTECOUNT(nNodes));
  for (i = 0; i < nNodes; i++)              /* parse each node */
  { 
    int nIndex = -1;
    _TOKEN(L_INT, DTR_MISSINGDTREEINDEX); /* node array index */
    nIndex = tok.value.l;        
    
    if (nIndex >= nNodes)
      _ON_ERR(DTR_BADDTREEINDEXRANGE);
    
    if (TESTMAP(aMapSeen, nIndex))          /* slot should be empty */
      _ON_ERR(DTR_DUPDTREEINDEX);
    _TOKEN(':', DTR_MISSINGLINEARCOLON);    /* ':' */
    
    /* get node */
    if ((nErr = GetDtreeNode(pFile, aBlocks, nBlocks, aVarDefs, nDim,
                             aNodes, nNodes, nIndex)) != DTR_NOERROR)
      _ON_ERR(nErr);
    
    SETMAP(aMapSeen, nIndex);               /* mark as seen */
  } 
  free(aMapSeen);
  aMapSeen = NULL;
             
  /*             
  ///////////////////////////////////
  // everything parsed OK!
  */
    
  /* create dtree */
  if ((*ppDtree = CreateDtree()) == NULL)
    _ON_ERR(DTR_MALLOCFAILED);
  
  (*ppDtree)->nDim = nDim;
  (*ppDtree)->nOutputIndex = nOutput;
  (*ppDtree)->aVarDefs = aVarDefs;
  (*ppDtree)->nLinearForms = nLinearForms;
  (*ppDtree)->aLinearForms = aLinearForms;
  (*ppDtree)->nBlocks = nBlocks;
  (*ppDtree)->aBlocks = aBlocks;
  (*ppDtree)->nNodes = nNodes;
  (*ppDtree)->aNodes = aNodes;
  
  goto PARSE_OK;  
        
  /* clear token helper defs */
  #undef _DBLTOKEN
  #undef _INTTOKEN
  #undef _TOKEN
  #undef _ON_ERR
  
PARSE_ERR:    
  DestroyVarDefArray(aVarDefs, nDim);
  DestroyDtreeNodeArray(aNodes);
  DestroyBlockArray(aBlocks, nBlocks);
  DestroyLinearFormArray(aLinearForms, nLinearForms);
  if (aMapSeen != NULL)
    free(aMapSeen);
  
PARSE_OK:    
  dtree_errno = nErr;
  return nErr;
}
               
/*               
/////////////////////////////////////////////////////////////////////
// variable definition parser         
*/
int GetVarDef(FILE* pFile, VARDEF* pVarDef)
{  
  char szName[MAXTOKLEN];
  TOKEN tok;
  
  if (GetToken(pFile, &tok) != IDENTIFIER)    /* identifier */
    return DTR_BADVARDEFIDENT;
  strcpy(szName, tok.value.sz);
  if (GetToken(pFile, &tok) != ':')           /* ':' */
    return DTR_MISSINGVARDEFCOLON;
  if (GetToken(pFile, &tok) != '[')           /* '[' */
    return DTR_MISSINGVARBOUNDSTART;
  if (GetDoubleToken(pFile, &tok) != L_DBL)  /* minimum */
    return DTR_BADVARBOUND;
  pVarDef->bound.dblMin = tok.value.dbl;
  if (GetToken(pFile, &tok) != ',')           /* comma */
    return DTR_MISSINGVARBOUNDCOMMA;
  if (GetDoubleToken(pFile, &tok) != L_DBL)  /* maximum */
    return DTR_BADVARBOUND;
  pVarDef->bound.dblMax = tok.value.dbl;
  if (pVarDef->bound.dblMin >= pVarDef->bound.dblMax) /* check range */
    return DTR_NEGATIVEVARBOUNDRANGE;
  if (GetToken(pFile, &tok) != ']')           /* ']' */
    return DTR_MISSINGVARBOUNDSTART;
  if (GetToken(pFile, &tok) != ';')           /* end of statement */
    return DTR_MISSINGVARDEFSEMI;
  
  SetVarDefName(pVarDef, szName);             /* set name */
  return DTR_NOERROR;
}
        
/*               
/////////////////////////////////////////////////////////////////////
// block parser               
*/
int GetBlock(FILE* pFile, LINEARFORM* aLinearForms, int nLinearForms, 
             BLOCK* pBlock)
{             
  int nErr = DTR_NOERROR;
  TOKEN tok;
  
  /* get min max tree */
  if ((nErr = GetMinMaxTree(pFile, aLinearForms, nLinearForms,
                            &(pBlock->pMinMaxTree))) != DTR_NOERROR)
    return nErr;

  if (GetToken(pFile, &tok) != ';')     /* end of statement */
    return DTR_MISSINGBLOCKSEMI;
          
  return DTR_NOERROR;
}

/*
/////////////////////////////////////////////////////////////////////
// min/max tree parser
*/
int GetMinMaxTree(FILE* pFile, LINEARFORM* aLinearForms, int nLinearForms, 
                  MINMAXNODE** ppMMN)
{             
  TOKEN tok;
  int nType;
  int nErr = DTR_NOERROR;

  /* lookahead to determine if this a min/max node or a linear elem index */
  nType = GetToken(pFile, &tok);
  if (nType != L_INT && nType != R_MIN && nType != R_MAX)
    return DTR_BADMINMAXNODE;               /* bad token type */

  *ppMMN = NULL;    
  if ((*ppMMN = CreateMinMaxNode()) == NULL)
    return DTR_MALLOCFAILED;
  
  if (nType == L_INT)                     /* linear form index */
  {                      
    int nLFIndex = tok.value.l;          
    if(nLFIndex >= nLinearForms || aLinearForms == NULL)
    {
      nErr = DTR_UNDEFINEDLINEARFORM;       /* linear form does not exist */
      goto MINMAXNODE_ERR;
    }
    
    (*ppMMN)->nType = DTREE_LINEAR;
    MMN_LFINDEX((*ppMMN)) = nLFIndex;
    goto MINMAXNODE_OK;
  }
  else                                      /* min/max node */
  {
    MINMAXNODE* pList = NULL;
    (*ppMMN)->nType = (short)((nType == R_MIN) ? DTREE_MIN : DTREE_MAX);
    
    if (GetToken(pFile, &tok) != '(')       /* '(' */
    {
      nErr = DTR_MISSINGMINMAXLISTSTART;
      goto MINMAXNODE_ERR;     
    }
    
    /* read child list */
    nType = GetToken(pFile, &tok);          /* lookahead */
    while (nType != ')')                    /* until we hit close paren */
    {
      MINMAXNODE* pChild = NULL;
      if (nType != L_INT && nType != R_MIN && nType != R_MAX) 
      {                                     /* lookahead revealed next token type is bad */
        nErr = DTR_BADMINMAXLISTELEM;
        goto MINMAXNODE_ERR;
      }
      PushbackToken(pFile, &tok);           /* pushback lookahead token */
      
      /* get child list element */
      if ((nErr = GetMinMaxTree(pFile, aLinearForms, nLinearForms, 
                                &pChild)) != DTR_NOERROR)
        goto MINMAXNODE_ERR;
                             
      /* append to child list */
      if (pList == NULL)
      {
        MMN_CHILDLIST((*ppMMN)) = pChild;
        pList = pChild;
      }
      else
      {
        pList->pNext = pChild;
        pList = pChild;
      }
      
      nType = GetToken(pFile, &tok);        /* get next lookahead token */
      if (nType == ',')                     /* skip if list elem separator */
        nType = GetToken(pFile, &tok);
    }
    
    if (MMN_CHILDLIST((*ppMMN)) == NULL)    /* can't have empty child list */
    {
      nErr = DTR_EMPTYMINMAXLIST;
      goto MINMAXNODE_ERR;      
    }
      
    goto MINMAXNODE_OK;                     /* success */
  }
  
MINMAXNODE_ERR:
  DestroyMinMaxNode(*ppMMN);
  *ppMMN = NULL;
                
MINMAXNODE_OK:
  return nErr;
}
    
/*    
/////////////////////////////////////////////////////////////////////
// linear form parser
*/
int GetLinearForm(FILE* pFile, int nDim, int nOutput, LINEARFORM* pLF, 
                  VARDEF* aVarDefs)
{ 
  int i;
  TOKEN tok;
  int nType;
  int bNeg = 0;                       /* negative flag */
  char* aVarMap = NULL;               /* map of parsed vars */
  int nErr = DTR_NOERROR;
  
  /* allocate var seen map */
  aVarMap = (char*)malloc(MAPBYTECOUNT(nDim));
  if (aVarMap == NULL)
  {
    nErr = DTR_MALLOCFAILED;    
    goto LINEARFORM_ERR;
  }
  memset(aVarMap, 0, MAPBYTECOUNT(nDim));
  
  /* init linear form to have all zero weights, centroid */
  for (i = 0; i < nDim; i++)
  {
    pLF->adblW[i] = 0;
    pLF->adblC[i] = 0;
  }

  /* get token */
  nType = GetToken(pFile, &tok);
  while(nType != ';')
  {
    double dblW;
    double dblC;
    int nVarIndex;
    int nCSign;
                
    if (nType == '-')                  
    {
      bNeg = !bNeg;
      nType = GetToken(pFile, &tok);
    }
                
    /* get weight */   
    if (nType == L_INT)
    {
      tok.value.dbl = tok.value.l;
      nType = L_DBL;
    }
    if (nType != L_DBL)
    {
      nErr = DTR_BADLINEARWEIGHT;
      goto LINEARFORM_ERR;
    }
    dblW = tok.value.dbl;
    if (bNeg) dblW = -dblW;
    
    /* get centroid expr. */
    if (GetToken(pFile, &tok) != '(')
    {
      nErr = DTR_MISSINGCENTROIDSTART;
      goto LINEARFORM_ERR;
    }
    /* get var name */
    if (GetToken(pFile, &tok) != IDENTIFIER)
    {
      nErr = DTR_MISSINGCENTROIDIDENT;
      goto LINEARFORM_ERR;
    }
    for (i = 0; i < nDim; i++)
    {
      /* var names are case sensitive */
      if (strcmp(tok.value.sz, aVarDefs[i].pszName) == 0)
        break;
    }
    if (i == nDim)
    {
      nErr = DTR_UNKNOWNVARIDENT;
      goto LINEARFORM_ERR;
    }
    nVarIndex = i;
    if (TESTMAP(aVarMap, nVarIndex))
    {
      nErr = DTR_DUPVARCENTROID;
      goto LINEARFORM_ERR;
    }    
    
    /* sign */
    nCSign = GetToken(pFile, &tok);
    if (nCSign != '-' && nCSign != '+')
    {
      nErr = DTR_MISSINGLINEARSIGN;
      goto LINEARFORM_ERR;
    }
    if (nCSign == '+') bNeg = 1;
    else bNeg = 0;
    
    /* get centroid */
    if (GetDoubleToken(pFile, &tok) != L_DBL)
    {
      nErr = DTR_BADLINEARCENTROID;
      goto LINEARFORM_ERR;
    }
    dblC = tok.value.dbl;
    if (bNeg) dblC = -dblC;
    
    /* get close paren */
    if (GetToken(pFile, &tok) != ')')
    {
      nErr = DTR_MISSINGCENTROIDEND;
      goto LINEARFORM_ERR;
    }
                 
    /* assign weight and centroid, and mark as seen */
    pLF->adblW[nVarIndex] = dblW;
    pLF->adblC[nVarIndex] = dblC;
    SETMAP(aVarMap, nVarIndex); 
    
    /* lookahead */
    nType = GetToken(pFile, & tok);
    if (nType != ';' && nType != '+' && nType != '-')
    {
      nErr = DTR_MISSINGLINEARSEMI;
      goto LINEARFORM_ERR;
    }
    else if (nType != ';')            /* + or - */
    {
      bNeg = nType == '-';            /* set neg flag for next time thru loop*/
      nType = GetToken(pFile, & tok); /* lookahead again */
    }
  }
  
  if (nOutput >= nDim)                /* verify output index */
    return DTR_BADOUTPUTRANGE;
                                   
  if (pLF->adblW[nOutput] == 0)       /* make sure output weight != 0 */
    return DTR_ZEROOUTPUTWEIGHT;
                     
  /* calc bias weight */
  pLF->dblBias = 0;
  for (i = 0; i < nDim; i++)
    pLF->dblBias -= pLF->adblW[i] * pLF->adblC[i];

LINEARFORM_ERR:

  free(aVarMap);
  return nErr;
}

/*
/////////////////////////////////////////////////////////////////////
// dtree node parser
*/
int GetDtreeNode(FILE* pFile, BLOCK* aBlocks, int nBlocks, VARDEF* aVarDefs,
                 int nDim, DTREENODE* aNodes, int nNodes, int nNode)
{
  int nType;
  int nErr = DTR_NOERROR;
  DTREENODE* pNode = NULL;
  TOKEN tok;
  
  if (nNode < 0 || nNode >= nNodes)
    return DTR_GENERIC;
  pNode = aNodes + nNode;
  
  /* lookahead determines whether this is a leaf or internal node */
  nType = GetToken(pFile, &tok);
  if (nType != R_BLOCK && nType != '(')
    return DTR_BADDTREEDEF;                 /* bad token type */
    
  if (nType == '(')                         /* internal node */
  { 
    int nVarIndex = 0;
    
    if (GetToken(pFile, &tok) != IDENTIFIER)
      return DTR_MISSINGDTREEVAR;
                                                                             
    /* get variable index */
    for (DNODE_VARINDEX(pNode) = 0; 
         DNODE_VARINDEX(pNode) < nDim; 
         DNODE_VARINDEX(pNode)++)
    {       
      /* var names are case sensitive */
      if (strcmp(aVarDefs[DNODE_VARINDEX(pNode)].pszName, tok.value.sz) == 0)
        break;
    }
    if (DNODE_VARINDEX(pNode) == nDim)      /* var ident not found */
      return DTR_UNKNOWNVARIDENT;
    
    if (GetToken(pFile, &tok) != S_LE)      /* '<=' */
      return DTR_MISSINGDTREELE; 
    if (GetDoubleToken(pFile, &tok) != L_DBL)  /* threshold */
      return DTR_MISSINGDTREETHRESHOLD;
    DNODE_THRESHOLD(pNode) = tok.value.dbl;
    if (DNODE_THRESHOLD(pNode) < aVarDefs[DNODE_VARINDEX(pNode)].bound.dblMin ||
        DNODE_THRESHOLD(pNode) > aVarDefs[DNODE_VARINDEX(pNode)].bound.dblMax)
      return DTR_BADTHRESHOLDRANGE;         /* check thresh range */
    
    if (GetToken(pFile, &tok) != ')')       /* ')' */
      return DTR_MISSINGDTREEPAREN;
    if (GetToken(pFile, &tok) != '?')       /* '?' */
      return DTR_MISSINGDTREEQUESTION;
    if (GetToken(pFile, &tok) != L_INT)   /* left child */
      return DTR_MISSINGDTREECHILDINDEX;                                          
    DNODE_LEFTINDEX(pNode) = tok.value.l;      
    if (GetToken(pFile, &tok) != ':')       /* ':' */
      return DTR_MISSINGDTREECHILDSEP;
    if (GetToken(pFile, &tok) != L_INT)   /* right child */
      return DTR_MISSINGDTREECHILDINDEX;                                          
    DNODE_RIGHTINDEX(pNode) = tok.value.l;      
                                          
    pNode->nLeaf = 0;                       /* mark as an internal node */

    if (DNODE_LEFTINDEX(pNode) >= nNodes)   /* out of range */
      return DTR_BADDTREEINDEXRANGE;
    if (DNODE_RIGHTINDEX(pNode) >= nNodes)  /* out of range */
      return DTR_BADDTREEINDEXRANGE;
    if (DNODE_LEFTINDEX(pNode) <= nNode)    /* circ ref */
      return DTR_CIRCULARDTREECHILDREF;
    if (DNODE_RIGHTINDEX(pNode) <= nNode)   /* circ ref */
      return DTR_CIRCULARDTREECHILDREF;
    
    /* fixup parent refs */
    aNodes[DNODE_LEFTINDEX(pNode)].nParentIndex = nNode;
    aNodes[DNODE_RIGHTINDEX(pNode)].nParentIndex = nNode;
  }
  else /* leaf node */
  { 
    if (GetToken(pFile, &tok) != L_INT)
      return DTR_MISSINGDTREEBLOCKINDEX;
  
    DNODE_BLOCKINDEX(pNode) = tok.value.l;
    if(DNODE_BLOCKINDEX(pNode) >= nBlocks)
      return DTR_UNDEFINEDREGION;   
               
    pNode->nLeaf = 1;                       /* mark as a leaf */
    
    /* fixup block back ref */
    aBlocks[DNODE_BLOCKINDEX(pNode)].nDtreeIndex = nNode;
  } 
  
  if (GetToken(pFile, &tok) != ';')         /* end of statement */
    return DTR_MISSINGDTREESEMI;

  return DTR_NOERROR;
}
  
/*  
/////////////////////////////////////////////////////////////////////
// lexers
*/
    
void SkipWS(FILE* pFile) 
{ 
  int c;
  
WS_SEARCH: 
  /* check for eof */
  if (feof(pFile))
    return;
 
  do
  {
    c = fgetc(pFile);
    if (c == '\n') dtree_lineno++;     /* inc linecount */
  } while (isspace(c) && !feof(pFile));
    
  /* check for eof */
  if (feof(pFile))
    return;
  
  /* put back last char */
  ungetc(c, pFile);
  
  /* check for comment */
  if (IsComment(pFile))
  {
    if (feof(pFile))
      return;       /* EOF in comment */
    goto WS_SEARCH; /* skip any more WS */
  }
}    
    
int GetToken(FILE* pFile, TOKEN* pTok)
{ 
  int c, c2;

  if (feof(pFile))
    return EOF;
  
  /* skip any white space */  
  SkipWS(pFile);    

  /* check for eof */
  if (feof(pFile))
    return EOF;    

  /* get token position */
  pTok->lFilePos = ftell(pFile);    
  
  
  /* lookahead at next 2 chars */
  c = fgetc(pFile);  
  c2 = fgetc(pFile);  
  fseek(pFile, pTok->lFilePos, SEEK_SET); /* go back */  
    
  /* get token */
  if (ISCSYMF(c))
    return GetIdentifier(pFile, pTok);
  else if (isdigit(c))
    return GetNumber(pFile, pTok);
  else
  {
    if (c == '-' && isdigit(c2))          /* part of number! */
      return GetNumber(pFile, pTok);
    else
      return GetSymbol(pFile, pTok);
  }
}

void PushbackToken(FILE* pFile, TOKEN* pTok)
{
  fseek(pFile, pTok->lFilePos, SEEK_SET);  
}
 
int GetSymbol(FILE* pFile, TOKEN* pTok)
{
  int c;
 
  if (feof(pFile))
    return EOF;
         
  c = fgetc(pFile);
  pTok->value.l = c;

  /* look ahead to get '<=' symbol */
  if (c == '<' && !feof(pFile))
  {     
    int c2 = fgetc(pFile);
    if (c2 == '=')
      pTok->value.l = S_LE;     /* <= symbol */
    else                  
      ungetc(c2, pFile);  /* not <=, putback */
  }
  
  return (int)pTok->value.l;
}

int GetIntegerToken(FILE* pFile, TOKEN* pTok)
{ 
  int c;                        
  char szTok[MAXTOKLEN];        /* token */
  int nLen = 0;                 /* length of token */
  int nType = ERR_TOKEN;
  
  /* skip any white space */  
  SkipWS(pFile);    
  
  if (feof(pFile))
    return EOF;
  
  /* keep scanning 'til error in integer def */
  szTok[0] = 0;
  while (!feof(pFile) && (nLen < (MAXTOKLEN - 1)))
  {
    c = fgetc(pFile);
    if (!isdigit(c) && ((c != '-') || nLen > 0))
    {
      /* bad char or neg sign in bad place... putback */
      ungetc(c, pFile); break;
    }
    szTok[nLen++] = (char)c;
  } 
  szTok[nLen] = 0;  /* terminate string */

  if (nLen > 0)
  {
    pTok->value.l = atol(szTok);
    nType = L_INT;
  }
          
  return nType;
}

int GetDoubleToken(FILE* pFile, TOKEN* pTok)
{
  int nType = GetToken(pFile, pTok);
  if (nType == L_INT)
  {
    pTok->value.dbl = pTok->value.l;  /* convert int to double */
    nType = L_DBL;
  }
  return nType;
}

int GetNumber(FILE* pFile, TOKEN* pTok)  
{
  int c;                        
  char szTok[MAXTOKLEN];        /* token */
  int nLen = 0;                 /* length of token */
  int nE = 0;                   /* e notation 'e' position */
  int nDots = 0;                /* number of decimal places seen */
  int nType = ERR_TOKEN;
  
  if (feof(pFile))
    return EOF;
  
  /* keep scanning 'til error in number def */
  szTok[0] = 0;
  while (!feof(pFile) && (nLen < (MAXTOKLEN - 1)))
  {
    c = fgetc(pFile);
    if (!isdigit(c))
    {
      if (c == '-')
      {     
        if ((!nE && nLen > 0) || (nE && (nLen != nE + 1)))
        {
          /* neg sign in bad place... putback */
          ungetc(c, pFile);
          break;
        }          
      }
      else if (c == '+')
      {     
        if (!(nE && (nLen == nE + 1)))
        {
          /* plus sign in bad place... putback */
          ungetc(c, pFile);
          break;
        }          
      }
      else if (c == '.')
      {
        nDots++;
        if (nDots > 1)  /* not part of the number, putback */
        {
          ungetc(c, pFile);
          break;
        }
      }
      else if (c == 'e' || c == 'E')
      {  
        if (!nE)
          nE = nLen;
        else            /* not part of the number, putback */
        { 
          ungetc(c, pFile);
          break;
        }
      }
      else              /* not part of the number, putback */
      {
        ungetc(c, pFile);
        break;
      }
    }
    szTok[nLen++] = (char)c;
  } 
  szTok[nLen] = 0;      /* terminate string */

  if (nDots > 0 || nE > 0)
  {           
    pTok->value.dbl = atof(szTok);
    nType = L_DBL;
  }
  else
  {
    pTok->value.l = atol(szTok);
    nType = L_INT;
  }
          
  return nType;
}

int GetIdentifier(FILE* pFile, TOKEN* pTok)
{
  int c;
  char szTok[MAXTOKLEN];        /* token */
  int nLen = 0;                 /* length of token */
  int nType = ERR_TOKEN;
  
  if (feof(pFile))
    return EOF;
  
  /* scanning 'til error in identifier def */
  szTok[0] = 0;
  while (!feof(pFile) && (nLen < (MAXTOKLEN - 1)))
  {
    c = fgetc(pFile);
    if ((nLen == 0 && !ISCSYMF(c)) ||
        (nLen > 0 && !ISCSYM(c)))
    {
      /* putback bad char */
      ungetc(c, pFile);
      break;
    }
    szTok[nLen++] = (char)c;
  }
  szTok[nLen] = 0;  /* terminate string */
  
  nType = IsReserved(szTok);
  if(nType != 0)
  {
    pTok->value.l = nType;
    return nType;
  }
  else
  {
    strcpy(pTok->value.sz, szTok);
    return IDENTIFIER;
  }
}
 
int IsReserved(char* psz)
{ 
  int i;
  for (i = 0; _nReserved[i] != 0; i++)
  {
    if (stricmp(psz, _szReserved[i]) == 0)
    {
      return _nReserved[i];
    }
  }
  
  return 0;
}

int IsComment(FILE* pFile)
{
  int c;
  long lFilePos;
  
  /* check eof */
  if (feof(pFile))
    return 0;       /* should never be EOF here */
    
  /* save filepos */
  lFilePos = ftell(pFile);
         
  /* get char */
  c = fgetc(pFile);
  
  if (c != '/' || feof(pFile))
  {
    ungetc(c, pFile); /* put it back */
    return 0;         /* not a comment */
  }
  
  /* might be a comment, get next char */
  c = fgetc(pFile);
  if (c == '/')       /* C++ style single line comment */
  {
    while (c != '\n' && !feof(pFile))
    {
      c = fgetc(pFile);
      if (c == '\n') dtree_lineno++; /* inc linecount */
    }
    return 1;         /* end of comment or eof in comment */
  }
    
  if (c == '*')       /* C style comment */
  {
  COMMENT_SEARCH:
    if (feof(pFile))
      return 1;       /* end of file in comment */
    do
    {
      c = fgetc(pFile);
      if (c == '\n') dtree_lineno++; /* inc linecount */
    } while (c != '*' && !feof(pFile));
    if (feof(pFile))
      return 1;       /* end of file in comment */
    c = fgetc(pFile);
    if (c == '\n') dtree_lineno++;   /* inc linecount */
    else if (c == '/')
      return 1;       /* end of comment */
      
    goto COMMENT_SEARCH;  /* keep searching */
  }
  
  /* not a comment, push back file pointer */
  fseek(pFile, lFilePos, SEEK_SET); /* go back */
  return 0;
}

/*  
/////////////////////////////////////////////////////////////////////
// Dtree exporter
*/

int ExportDtreeFile(FILE* pFile, const char* pszFileName, DTREE* pDtree)
{
  time_t timeNow;
  int i;        
  
  static char _szHeader[] = "// %s exported on %s\n"
                            "// ALN Decision Tree file format v1.0 (C)1994 Dendronic Decisions Limited\n\n";
                                        
  static char _szVersion[] = "VERSION = 1.0;\n";                                        
  static char _szVarDefs[] = "VARIABLES = %ld;\n";
  static char _szOutput[] = "OUTPUT = %s;\n";
  static char _szLinearForms[] = "LINEARFORMS = %ld;\n";
  static char _szBlocks[] = "BLOCKS = %ld;\n";
  static char _szDtree[] = "DTREE = %ld;\n";
  
  if (pFile == NULL || pDtree == NULL)
    return DTR_GENERIC;
  
  /* write header */
  time(&timeNow);
  fprintf(pFile, _szHeader, pszFileName, asctime(localtime(&timeNow)));
  
  /* write version */
  fprintf(pFile, _szVersion);
  
  /* write variables */
  fprintf(pFile, _szVarDefs, pDtree->nDim);
  if (pDtree->aVarDefs == NULL)
    return DTR_GENERIC;
  for (i = 0; i < pDtree->nDim; i++)
  { 
    static char _szVarDef[] = "%s : [%.19lg, %.19lg];\n";
    VARDEF* pV = pDtree->aVarDefs + i;
    fprintf(pFile, _szVarDef, pV->pszName, pV->bound.dblMin, pV->bound.dblMax);
  }
  
  /* write output */
  if (pDtree->nOutputIndex >= pDtree->nDim)
    return DTR_GENERIC;
  fprintf(pFile, _szOutput, pDtree->aVarDefs[pDtree->nOutputIndex].pszName);
  
  /* write linear forms */
  fprintf(pFile, _szLinearForms, pDtree->nLinearForms);
  if (pDtree->aLinearForms == NULL)
    return DTR_GENERIC;
  for (i = 0; i < pDtree->nLinearForms; i++)
  {      
    int j;
    int nVarsOut;
    static char _szLinearForm[] = "%ld : ";
    static char _szVar[] = "%.19lg (%s %c %.19lg)";
    static char _szSign[] = " %c ";
    LINEARFORM* pLF = pDtree->aLinearForms + i;
    if (pLF == NULL)
      return DTR_GENERIC;
    fprintf(pFile, _szLinearForm, i);
                        
    nVarsOut = 0;     /* no vars output yet */
    for (j = 0; j < pDtree->nDim; j++)
    {                  
      int nWSign, nCSign; 
      double dblW, dblC;
    
      if (pLF->adblW[j] == 0)
        continue;               
      
      dblW = pLF->adblW[j];
      dblC = pLF->adblC[j];

      /* change signs in expression so weight and centroid are always positive */
      if (dblW < 0 && nVarsOut > 0) /* only change sign on weight if not first var */
      {
        nWSign = '-';
        dblW = -dblW;
      }
      else
        nWSign = '+'; 
        
      if (dblC < 0)
      {
        nCSign = '+';
        dblC = -dblC;
      }
      else
        nCSign = '-';
      
      if (nVarsOut > 0)
        fprintf(pFile, _szSign, nWSign);
        
      fprintf(pFile, _szVar, dblW, pDtree->aVarDefs[j].pszName, nCSign, dblC);
      nVarsOut++;
    }
    /* write statement end */
    fputc(';', pFile);
    fputc('\n', pFile);
  }
  
  /* write blocks */
  fprintf(pFile, _szBlocks, pDtree->nBlocks);
  if (pDtree->aBlocks == NULL)
    return DTR_GENERIC;
  for (i = 0; i < pDtree->nBlocks; i++)
  {                                        
    int nErr = DTR_NOERROR;
    static char _szBlock[] = "%ld : ";
    BLOCK* pBlock = pDtree->aBlocks + i;
    fprintf(pFile, _szBlock, i);
    
    /* write minmax node */
    if ((nErr = WriteMinMaxTree(pFile, pDtree, pBlock->pMinMaxTree, 0)) != DTR_NOERROR)
      return nErr;
    
    /* write statement end */
    fputc(';', pFile);
    fputc('\n', pFile);
  }

  /* write dtree */
  fprintf(pFile, _szDtree, pDtree->nNodes);
  if (pDtree->aNodes == NULL)
    return DTR_GENERIC;
  for (i = 0; i < pDtree->nNodes; i++)
  {                                        
    static char _szNode[] = "%ld : ";
    DTREENODE* pNode = pDtree->aNodes + i;
    fprintf(pFile, _szNode, i);
    
    /* is this a leaf ? */
    if (pNode->nLeaf != 0)
    {
      static char _szLeaf[] = "block %ld";
      fprintf(pFile, _szLeaf, DNODE_BLOCKINDEX(pNode));
    }
    else  /* internal node */
    {
      static char _szInternal[] = "(%s <= %.19lg) ? %ld : %ld";
      fprintf(pFile, _szInternal, 
              pDtree->aVarDefs[DNODE_VARINDEX(pNode)].pszName,
              DNODE_THRESHOLD(pNode), DNODE_LEFTINDEX(pNode), 
              DNODE_RIGHTINDEX(pNode));
    }
    
    /* write statement end */
    fputc(';', pFile);
    fputc('\n', pFile);
  }
            
  return DTR_NOERROR;
}

int WriteMinMaxTree(FILE* pFile, DTREE* pDtree, MINMAXNODE* pMMN, int nIndent)
{ 
  if (pFile == NULL || pDtree == NULL || pMMN == NULL)
    return DTR_GENERIC;
                          
  if (pMMN->nType == DTREE_LINEAR)
  {        
    static char _szIndex[] = "%ld";
    
    if (MMN_LFINDEX(pMMN) >= pDtree->nLinearForms || 
        pDtree->aLinearForms == NULL)
      return DTR_UNDEFINEDLINEARFORM;
    
    fprintf(pFile, _szIndex, MMN_LFINDEX(pMMN));
  }
  else if (pMMN->nType == DTREE_MIN || pMMN->nType == DTREE_MAX)
  {   
    int nErr = DTR_NOERROR;
    static char _szMin[] = "MIN";
    static char _szMax[] = "MAX";
    static char _szSep[] = ", ";
    static char _szNL[] = "\n%*s";
    MINMAXNODE* pList = MMN_CHILDLIST(pMMN);                    
    if (pMMN->nType == DTREE_MIN)
      fprintf(pFile, _szMin);
    else
      fprintf(pFile, _szMax);
               
    fputc('(', pFile);
    while(pList)
    { 
      if ((nErr = WriteMinMaxTree(pFile,pDtree,pList,nIndent+4)) != DTR_NOERROR)
        return nErr;
                
      if (pList->pNext != NULL) 
        fprintf(pFile, _szSep);
      
      pList = pList->pNext;
    }     
    fputc(')', pFile);
  }
  else
    return DTR_GENERIC;
                          
  return DTR_NOERROR;
}

