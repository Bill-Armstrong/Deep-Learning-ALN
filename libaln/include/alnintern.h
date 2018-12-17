// ALN Library (libaln)
// Copyright (C) 2018 William W. Armstrong.
// file: alnintern.h
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


// Functions used only internally
void splitControl(ALN*, double); // if average variance error of a piece is high, splitting is prevented
void doSplits(ALN*, ALNNODE*, double); //does the recursion of splitcontrol
void zeroSplitValues(ALN*, ALNNODE*);  // sets the square error to zero in each LFN
void splitUpdateValues(ALN*); // accumulates the training square error and number of hits on each linear piece
void splitNoiseSetVAR(ALN*); // accumulates the variance square error and number of hits on each linear piece
void dodivideTR(ALN*, ALNNODE*); // divides the total square training set errors of the pieces by their hit count
void dodivideVAR(ALN*, ALNNODE*); // divides the sum of noise variance samples of the pieces by their respective hit counts

// thread procedures
UINT TakeActionProc(LPVOID pParam);  // separate thread

//global variables used only internally
extern BOOL bTrainingAverage;// Switch to tell fillvector whether get a training vector or compute an average
//extern int nMaxEpochs;
extern int nNumberLFNs;  // used to control the epoch size, which should be proportional to this
extern int nColsAuxVariance;
extern int nColsAuxTest;
extern BOOL bDecimal; // means numbers could have a decimal point
extern BOOL bComma;   // means numbers could have a comma
extern int nDim;
extern int nOutputIndex;
extern long nRowsPP;          // The number of rows in the PreprocessedDataFile
extern long nRowsTV;          // The number of rows in the TVfile (Training & Variance File)
extern double* adblEpsilon;         // An array of doubles holding the tolerances of the inputs.
extern double* adblMinVar;          // Array of minima of the variables
extern double* adblMaxVar;          // Array of maxima of the variables
extern double* adblStdevVar;        // Standard deviations of the variables
extern double  dblTrainErr;         // Set in cmyaln.h at the end of training
extern double  dblLinRegErr;        // The error of linear regression for use in upper-bounding output tolerance
extern double* adblLRW; // stores an ALN weight approximation from linear regression
extern double* adblLRC; // ditto for centroids

// Files used only internally
extern FILE *fpData;                 // The data file which contains all data.
extern FILE *fpProtocol;             // the file to record the results of the experiment. Rename to save.
extern FILE *fpReplacement;
extern CDataFile UNfile;             // copy of the data file, but with missing values replaced by special number
extern CDataFile PreprocessedDataFile; // The preprocessed file created from the Universal file
extern CDataFile TVfile;             // The file used for training approximants.  It is split in two to determine noise variance

