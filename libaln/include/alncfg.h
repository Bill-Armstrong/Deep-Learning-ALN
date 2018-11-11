/*
// ALN Library
// ALNfit Learning Engine for approximation of functions defined by samples.
// Copyright (C) 2018 William W. Armstrong.

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// Version 3 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

// For further information contact 
// William W. Armstrong
// 3624 - 108 Street NW
// Edmonton, Alberta, Canada  T6J 1B4

*/

/*
// alncfg.h
*/

#ifndef __ALNCFG_H__
#define __ALNCFG_H__

/*
/////////////////////////////////////////////////////////////////////////////
// calling conventions
*/

/* define inline as static for non-C++ compilers */
#ifndef NODEFINE_INLINE

#if ( !defined(__cplusplus) && !defined(inline) )
#define inline static
#endif

#endif  /* NODEFINE_INLINE */

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

/* define __fastcall for non-Microsoft compilers */
#ifndef NODEFINE_FASTCALL

#if	( !defined(_MSC_VER) && !defined(__fastcall) )
#define __fastcall
#endif

#endif /* NODEFINE_FASTCALL */


#define ALNAPI __stdcall


/*
/////////////////////////////////////////////////////////////////////////////
// DLL function exports - define ALNDLL if linking to libalndll(d).DLL
*/

#if defined(ALNDLL) && !defined(ALNIMP)
#define ALNIMP __declspec(dllimport)
#elif !defined(ALNIMP)
#define ALNIMP
#endif

/*
/////////////////////////////////////////////////////////////////////////////
// library files (Microsoft compilers only)
*/

#ifndef ALN_NOFORCE_LIBS

#ifdef _MSC_VER

#ifdef _DEBUG

  #ifdef ALNDLL
  #pragma comment(lib, "libalndlld.lib")      /* ALNDLL */
  #else

    #ifdef _MT
      #ifdef _DLL
      #pragma comment(lib, "libalndmtd.lib")  /* MT DLL */
      #else  /* !_DLL */
      #pragma comment(lib, "libalnmtd.lib")   /* MT */
      #endif /* _DLL */
    #else
      #pragma comment(lib, "libalnd.lib")
    #endif /* _MT */

  #endif /* ALNDLL */


#else   /* !_DEBUG */

  #ifdef ALNDLL
  #pragma comment(lib, "libalndll.lib")       /* ALNDLL */
  #else

    #ifdef _MT
      #ifdef _DLL
      #pragma comment(lib, "libalndmt.lib")   /* MT DLL */
      #else  /* !_DLL */
      #pragma comment(lib, "libalnmt.lib")    /* MT */
      #endif /* _DLL */
    #else
      #pragma comment(lib, "libaln.lib")
    #endif /* _MT */

  #endif /* ALNDLL */

#endif  /* _DEBUG */

#endif  /* _MSC_VER */

#endif  /* ALN_NOFORCE_LIBS */

/*
///////////////////////////////////////////////////////////////////////////////
*/

#endif  /* __ALNCFG_H__ */
