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
// alnaddtreestring.cpp

#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"
#include <ctype.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


/////////////////////////////////////////////////////////////////////////////
// tree parsing

// general form of tree string (notation is YACC like)

// treestring : <null>                  // nothing... leave as a single hyperplane
//            | L <null>                // nothing... leave as a single hyperplane
//            | tree_exp <null>;        // min / max expression
//
// tree_exp   : minmax_type '(' child_list ')';
//
// minmax_type  : MIN | MAX;
//
// child_list : child_exp
//            | child_list ',' child_exp;
//
// child_exp  : L_FANIN 
//            | tree_exp;


// token types

#define ERROR -1
#define T_L 0
#define T_MIN 1
#define T_MAX 2
#define T_L_FANIN 3
#define T_BEGIN_LIST 4
#define T_END_LIST 5
#define T_COMMA 6
#define T_ENDOFSTRING 7

// parse states

#define P_TREESTRING 0
#define P_TREE_EXP 1
#define P_MINMAX_TYPE 2
#define P_CHILD_LIST 3
#define P_CHILD_EXP 4

// parse tree structure
typedef struct tagALNPARSE
{
  int nType;                        // L, MIN, MAX, L_FANIN
  int nChildren;                    // # children for MIN/MAX
  struct tagALNPARSE** apChildren;  // child nodes
} ALNPARSE;

// tokenizers
static const char* MinMaxTokenParse(const char* psz, int* pnTokenType);
static const char* FaninTokenParse(const char* psz, int* pnTokenType, int* pnFanin);
static const char* GetNextTreeToken(const char* psz, int* pnTokenType, int* pnFanin);

// memory mgmt
static void FreeParse(ALNPARSE* pParse);
static ALNPARSE* AllocParse(int nType, int nChildren);
static ALNPARSE** ReAllocParseChildren(ALNPARSE* pParse, int nChildren);

// main parser
static const char* DoParseTreeString(ALNPARSE* pParse, int* pnCurrentToken, const char* psz, int* pnState,
                                     int nFanin);

// maps parse tree onto ALN structure
static int BuildParseTree(ALN* pALN, ALNNODE* pParent, ALNPARSE* pParse);


ALNIMP int ALNAPI ALNAddTreeString(ALN* pALN, ALNNODE* pParent, 
                                   const char* pszTreeString, int* pnParsed)
{
  // param variance
  if (pALN == NULL || pALN->pTree == NULL)
    return ALN_GENERIC;
  if (!NODE_ISLFN(pALN->pTree))
    return ALN_GENERIC;

  // alloc initial parse tree node to match ALN' single LFN
  ALNPARSE* pParse = AllocParse(T_L, 0); 
  if (pParse == NULL)
    return ALN_OUTOFMEM;

  // parse the tree string
  int nParseState = P_TREESTRING;
  int nCurrentToken = ERROR;
  const char* psz = DoParseTreeString(pParse, &nCurrentToken, pszTreeString, &nParseState, 0);

  // count number of chars parsed
  if (pnParsed != NULL)
    *pnParsed = psz - pszTreeString;

  // check final state
  int nReturn = ALN_GENERIC;
  if (nParseState == P_TREESTRING)
  {
    // the whole string should have been parsed by this point
    ASSERT((size_t)(psz - pszTreeString) == strlen(pszTreeString));
    
    // build the corresponding ALN
    nReturn = BuildParseTree(pALN, pALN->pTree, pParse);
      // this may fail, leaving a partially constructed ALN
  }
   
  // free the parse tree
  FreeParse(pParse);
 
  return nReturn;
}

/////////////////////////////////////////////////////////////////////////////
// helpers

// fanin tokenizer... pnTokenType contains token type, returns pointer to next token 
static const char* FaninTokenParse(const char* psz, int* pnTokenType, int* pnFanin)
{
  ASSERT(psz);
  ASSERT(pnTokenType != NULL);
  ASSERT(pnFanin != NULL);

  *pnTokenType = ERROR;
  *pnFanin = -1;

  char* pszEnd;
	char chFirst = *psz;
	
  long l = strtol(psz, &pszEnd, 10);
  psz = pszEnd;   // points to char that stopped scan

	if (l <= 0)
  {
		return psz;   // could not convert to positive int
  }
  else
  {
    *pnTokenType = T_L_FANIN;
    *pnFanin = l;
  	return psz;
  }
}


// min/max tokenizer... pnTokenType contains token type, returns pointer to next token 
static const char* MinMaxTokenParse(const char* psz, int* pnTokenType)
{
  ASSERT(psz);
  ASSERT(pnTokenType != NULL);

  *pnTokenType = ERROR;

  if (*psz != 'M' && *psz != 'm') // first char is M?
    return psz;

  // next char
  psz++;
  if (*psz == 'A' || *psz == 'a') // max ?
  {
    psz++;
    if (*psz == 'X' || *psz == 'x')
    {
      *pnTokenType = T_MAX;
      return psz + 1;
    }
  }
  else if(*psz == 'I' || *psz == 'i') // min ?
  {
    psz++;
    if (*psz == 'N' || *psz == 'n')
    {
      *pnTokenType = T_MIN;
      return psz + 1;
    }
  }

  return psz;
}


// tokenizer... pnTokenType contains token type, returns pointer to next token 
// pnFanin will contain integer fanin if token is L_FANIN
static const char* GetNextTreeToken(const char* psz, int* pnTokenType, int* pnFanin)
{
  ASSERT(psz);
  ASSERT(pnTokenType != NULL);
  ASSERT(pnFanin != NULL);

  *pnTokenType = ERROR;
  *pnFanin = -1;

  // look for tokens
 
  // white space
  while (isspace(*psz))
  {
    psz++;
  }

  // non-white space
  switch(*psz)
  {
    case '\0':
      {
        *pnTokenType = T_ENDOFSTRING;
        return psz;
      }

    case ',':
      {
        *pnTokenType = T_COMMA;
        return psz+1;
      }
    case '(':
      {
        *pnTokenType = T_BEGIN_LIST;
        return psz+1;
      }
    case ')':
      {
        *pnTokenType = T_END_LIST;
        return psz+1;
      }
    case 'L':
    case 'l':
      {
        *pnTokenType = T_L;
        return psz+1;
      }
    case 'M':
    case 'm':
      {
        return MinMaxTokenParse(psz, pnTokenType);
      }
    default:
      {
        if (isdigit(*psz))
          return FaninTokenParse(psz, pnTokenType, pnFanin);
      }
  }
  
  return psz;
}


// memory mgmt

static void FreeParse(ALNPARSE* pParse)
{
  ASSERT(pParse);

  // recurse
  if (pParse->apChildren != NULL)
  {
    ASSERT(pParse->nChildren > 0);
    for (int i = 0; i < pParse->nChildren; i++)
    {
      if (pParse->apChildren[i] != NULL)
        FreeParse(pParse->apChildren[i]);
    }
    free(pParse->apChildren);
  }
  free(pParse);
}

static ALNPARSE* AllocParse(int nType, int nChildren)
{
  ALNPARSE* pParse = (ALNPARSE*)calloc(1, sizeof(ALNPARSE));
  if (pParse == NULL)
    return NULL;

  pParse->nType = nType;
   
  if (nChildren > 0)
  {
    pParse->nChildren = nChildren;
    pParse->apChildren = (ALNPARSE**)calloc(nChildren, sizeof(ALNPARSE*));
    if (pParse->apChildren == NULL)  // failed
    {
      free(pParse);
      return NULL;
    }
  }

  return pParse;
}

static ALNPARSE** ReAllocParseChildren(ALNPARSE* pParse, int nChildren)
{
  ASSERT(pParse);

  // free extra children
  if (nChildren < pParse->nChildren)
  {
    for (int i = nChildren; i < pParse->nChildren; i++)
    {
      if (pParse->apChildren[i] != NULL)
        FreeParse(pParse->apChildren[i]);
    }
  }

  // realloc
  ALNPARSE** apChildren = (ALNPARSE**)realloc(pParse->apChildren, nChildren*sizeof(ALNPARSE**));
  if (apChildren == NULL)  // failed
    return NULL;

  // init new children
  pParse->apChildren = apChildren;
  if (nChildren > pParse->nChildren)
  {
    for (int i = pParse->nChildren; i < nChildren; i++)
    {
      pParse->apChildren[i] = NULL;
    }
  }
  pParse->nChildren = nChildren;

  return pParse->apChildren;
}


// main parser

static const char* DoParseTreeString(ALNPARSE* pParse, int* pnCurrentToken, const char* psz, int* pnState,
                                     int nFanin)
{
  ASSERT(pParse && pnCurrentToken && psz && pnState);

  switch(*pnState)
  {
    case P_TREESTRING:
      {
        // treestring : <null>                  // nothing... leave as a single hyperplane
        //            | L <null>                // nothing... leave as a single hyperplane
        //            | tree_exp <null>;        // min / max expression

        psz = GetNextTreeToken(psz, pnCurrentToken, &nFanin);
        switch(*pnCurrentToken)
        {
          case T_ENDOFSTRING:         // nothing to do
            break;
          case T_L:
            {
              // look for appropriate terminator
              psz = GetNextTreeToken(psz, pnCurrentToken, &nFanin);
              if (*pnCurrentToken != T_ENDOFSTRING)
                *pnState = ERROR;

              break;
            }
          default:
            {
              // assume its a tree expression... parse from current token
              int nSubState = P_TREE_EXP;
              psz = DoParseTreeString(pParse, pnCurrentToken, psz, &nSubState, 0);
              if (nSubState != P_TREE_EXP)
                *pnState = ERROR;
              else
              {
                // look for appropriate terminator
                psz = GetNextTreeToken(psz, pnCurrentToken, &nFanin);
                if (*pnCurrentToken != T_ENDOFSTRING)
                  *pnState = ERROR;
              }

              break;
            }
        }
        return psz;
      }
    case P_TREE_EXP:
      {
        // tree_exp   : minmax_type '(' child_list ')';
        
        // parse minmax type from current token
        int nSubState = P_MINMAX_TYPE;
        psz = DoParseTreeString(pParse, pnCurrentToken, psz, &nSubState, 0);
        if (nSubState != P_MINMAX_TYPE)
          *pnState = ERROR;
        else 
        {
          // get BEGIN_LIST token
          psz = GetNextTreeToken(psz, pnCurrentToken, &nFanin);
          if (*pnCurrentToken != T_BEGIN_LIST)
            *pnState = ERROR;
          else
          {
            // CHILD_LIST
            
            // get next token
            psz = GetNextTreeToken(psz, pnCurrentToken, &nFanin);

            // parse child list from current token
            nSubState = P_CHILD_LIST;
            psz = DoParseTreeString(pParse, pnCurrentToken, psz, &nSubState, nFanin);
            if (nSubState != P_CHILD_LIST)
              *pnState = ERROR;
            else
            {
              // current token should now be END_LIST
              if (*pnCurrentToken != T_END_LIST)
                *pnState = ERROR;
            }
          }
        }
        return psz;
      }
    case P_MINMAX_TYPE:
      {
        // minmax_type  : MIN | MAX;
        switch(*pnCurrentToken)
        {
          case T_MIN:
          case T_MAX:
            pParse->nType = *pnCurrentToken;
            break;

          default:
            *pnState = ERROR; // invalid token
        }
        return psz;
      }
    case P_CHILD_EXP:
      {
        // child_exp  : L_FANIN 
        //            | tree_exp;

        switch(*pnCurrentToken)
        {
          case T_L_FANIN:
            {
              // alloc new children
              ASSERT(nFanin > 0);
              ASSERT(pParse->nType == T_MIN || pParse->nType == T_MAX);
              int nOldChildren = pParse->nChildren;
              int nTotal = nOldChildren + nFanin;
              ALNPARSE** apChildren = ReAllocParseChildren(pParse, nTotal);
              if (apChildren == NULL)
                *pnState = ERROR;
              else
              {
                for (int i = nOldChildren; i < nTotal; i++)
                {
                  apChildren[i] = AllocParse(T_L, 0);
                  if (apChildren[i] == NULL)
                  {
                    *pnState = ERROR;
                  }
                }
              }
              break;
            }
          case T_MIN:
          case T_MAX:
            {
              // alloc new node
              ALNPARSE* pParseChild = AllocParse(*pnCurrentToken, 0);
              if (pParseChild == NULL)
                *pnState = ERROR;
              else
              {
                // parse tree expression from current token
                int nSubState = P_TREE_EXP;
                psz = DoParseTreeString(pParseChild, pnCurrentToken, psz, &nSubState, nFanin);
                if (nSubState != P_TREE_EXP)
                {
                  FreeParse(pParseChild);
                  *pnState = ERROR;
                }
                else  // add to parent
                {
                  ASSERT(pParse->nType == T_MIN || pParse->nType == T_MAX);
                  int nOldChildren = pParse->nChildren;
                  int nTotal = nOldChildren + 1;
                  ALNPARSE** apChildren = ReAllocParseChildren(pParse, nTotal);
                  if (apChildren == NULL)
                  {
                    FreeParse(pParseChild);
                    *pnState = ERROR;
                    break; // out of list loop
                  }
                  else
                  {
                    apChildren[nOldChildren] = pParseChild;
                  }
                }
              }
              break;
            }
          default:
            {
              *pnState = ERROR;
            }
        }
        return psz;
      }
    case P_CHILD_LIST:
      {
        // child_list : child_exp
        //            | child_list ',' child_exp;

        while (*pnState == P_CHILD_LIST)
        {
          // CHILD_EXP

          // parse child exp from current token
          int nSubState = P_CHILD_EXP;
          psz = DoParseTreeString(pParse, pnCurrentToken, psz, &nSubState, nFanin);
          if (nSubState != P_CHILD_EXP)
          {
            *pnState = ERROR;
            break;  // out of list loop
          }

          // get next token                        
          psz = GetNextTreeToken(psz, pnCurrentToken, &nFanin);
          if (*pnCurrentToken == T_COMMA)
          {
            // get next token in list
            psz = GetNextTreeToken(psz, pnCurrentToken, &nFanin);
          }
          else  // not a comma... end of list(?)
          {
            ASSERT(*pnState == P_CHILD_LIST);
            break;  // out of list loop
          }
        }
        return psz;
      }
  }

  // default
  *pnState = ERROR;
  return psz;
}


// build ALN from parse tree
static int BuildParseTree(ALN* pALN, ALNNODE* pParent, ALNPARSE* pParse)
{
  if (pParse->nType == T_L)
    return ALN_NOERROR;  // nothing to do

  // build children of current node
  ALNNODE** apChildren;
  apChildren = (ALNNODE**)calloc(pParse->nChildren, sizeof(ALNNODE*));
  if (apChildren == NULL)
    return ALN_OUTOFMEM;
   
  int nReturn = ALNAddLFNs(pALN, 
                           pParent, 
                           (pParse->nType == T_MIN) ? GF_MIN : GF_MAX,
                           pParse->nChildren, 
                           apChildren);
  
  if (nReturn == ALN_NOERROR)
  {
    // recurse over children
    for (int i = 0; i < pParse->nChildren; i++)
    {
      nReturn = BuildParseTree(pALN, apChildren[i], pParse->apChildren[i]);
      if (nReturn != ALN_NOERROR)
      {
        break;  // out of child loop
      }
    }
  }

  // clean up
  free(apChildren);

  return ALN_NOERROR;
}