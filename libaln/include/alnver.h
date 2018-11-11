/*
// ALN Library
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

// alnver.h
*/


#ifndef __ALNVER_H__
#define __ALNVER_H__

/* aln library version */

#ifndef ALNLIBVER

#define ALNLIBVER 0x00020001
// Version 0x00010009 improved speed and DTREE building
// Version 0x00020001 includes better DTREE optimization
#endif  /* ALNLIBVER */


/* aln structure version */

#ifndef ALNVER

// Version 0x00030008->9 allowed trees to be of arbitrary structure not just alternating AND/OR
// Version 0x00030009->10 changed SEy to SEE to avoid confusion.

#define ALNVER 0x00030010

#endif  /* ALNVER */


/*
///////////////////////////////////////////////////////////////////////////////
*/

#endif  /* __ALNVER_H__ */
