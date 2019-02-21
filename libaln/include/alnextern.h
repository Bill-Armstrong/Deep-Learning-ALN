// ALN Library (libaln)
// Copyright (C) 2018 William W. Armstrong
// file: alnextern.h
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
#include "aln.h"
#include "alnpp.h"
#include ".\cmyaln.h"
#include "datafile.h"
#define ALNAPI __stdcall

// functions used externally
void ALNAPI ALNfitSetup();  // sets up automatic noise measurement, bagging and evaluation
int ALNAPI analyzeInputFile(char * szDataFileName, int * nHeaderLines, long * nrows, int * ncols, BOOL bPrint);	// Analyzes the given input data file 
void ALNAPI analyzeTV();	// Computes the standard deviations of the variables in the TVset.
void ALNAPI getTVfile();	// The TVfile created from the PreprocessedDataFile is read in
void ALNAPI getTSfile();	// The TSfile created from the PreprocessedDataFile is read in 
void ALNAPI doLinearRegression();	// This does a truncated linear regression fit to get an upper bound on noise
void ALNAPI createNoiseVarianceFile(); // Creates samples estimating the noise variance.
//void ALNAPI trainNoiseVarianceALN(); // Trains on the noise variance samples to smooth them.
void ALNAPI approximate();	// This creates one or more approximant ALNs using the weight bounds found above.  These are averaged in bagging later.
void ALNAPI reportFunctions();	// reports on the trained function ALNs with stats and plots
void ALNAPI evaluate();	// Evaluate an existing DTREE on the data file after preprocessing
void ALNAPI cleanup();	// destroys allocated items no longer needed
void ALNAPI outputTrainingResults();	// outputs the results of training
void ALNAPI trainAverage();	// trains an average of several smoothed ALNs (i.e. with fillets)
void ALNAPI constructDTREE(int);	// constructs a DTREE from the average ALN from bagging
int ALNAPI analyzeauxiliaryfile(char * szAuxiliaryFileName, int * pAuxheaderlines, long * pAuxrows, int * pAuxcols, BOOL bPrint);
void ALNAPI MakeAuxNumericalFile(char * szAuxiliaryFileName,int nHeaderLinesAuxiliary,long nAuxRows, int nAuxCols, CDataFile & AuxNumericalFile);
void ALNAPI MakeAuxALNinputFile(const CDataFile & AuxNumericalFile, CDataFile & AuxALNinputFile, long nAuxRows, long * nRowsAuxALNinputFile);

// global variables used externally
extern BOOL bClassify;	// This is TRUE if the user chose a classification problem, and FALSE for regression
extern BOOL bRegress;		// This is FALSE iff the TV file outputs are all 0's and 1's
extern BOOL bTrain;			// This is TRUE for training and FALSE for evaluation.
extern BOOL bJitter;    // This is re-set under the options dialog
extern BOOL bReplaceUndefined; // Makes a new column at right if FALSE, otherwise replaces undefined output values
extern BOOL bDiagnostics;  // For controlling printout of diagnostic files
extern BOOL bTimePrefixes;	// Time prefixes are added at the front of some files to distinguish different runs
extern BOOL bPrint;  // controls printing of the input files, may be changed in options dialog
extern int nALNs;	// the number of ALNs trained in bagging which are later averaged
extern int nDTREEDepth; // level of partitioning of the input space to make a DTREE of several ALNs on the parts
extern double dblEvalRMSError;
extern int nEvalMisclassifications;
extern double dblEvalMisclassificationPercent;
extern void fillvector(double *, CMyAln*); // inputs a training data vector, which can be obtained by averaging several ALNs in bagging
extern CDataFile AuxALNinputFile;
extern CDataFile TSfile; // the test file with data that have not been used in training
extern long nRowsTS;
extern double dblSetTolerance; // This used to be for setting the noise variance level when it was constant.  Now it is not used.
extern BOOL bEstimateNoiseVariance; // if TRUE we create noise variance samples.
extern double  dblVarianceErr;    // Set equal to the rmse in the variance step

// variables used and viewed or set in a dialog
// nDim is always the number of ALN input variables plus one for the output
// this is true even when an output column is missing from the data file.
extern int nColsUniv;         // number of columns in the universal file
extern long nRowsUniv;    // number of rows in the universal file
extern int nHeaderLinesUniv;  // the number of headerlines found in the universal file
extern char hdrline[10][2000]; // up to nine header lines plus a buffer
extern char varname[101][32]; //space for datafile column names up to 31 characters long, NULL terminated
extern int nheaderlineitems[10];  // this counts the number of fields read in each header line
extern int nheaderlinewithvars;		//  this is the one of the last two headerlines with the variable names
																	//  however if neither has the correct number of strings this is 10
extern int nALNinputs;	// this is the number of input variables to the ALN, including the desired output. Maybe not all columns are associated with ALN inputs.
extern int  nInputCol[101];	// the indices of datafile columns serving as inputs to the ALNs (columns may be repeated)
extern int  nLag[101]; // the number of rows back (i.e. higher up) from the current row of the output that a column is sampled
extern double dblMax;  // maximum number allowed for weight
extern double dblMaxWeight[101]; // the maximum weight a given input can have
extern double dblMinWeight[101]; // the minimum weight a given input can have
extern double dblImportance[101]; // the importance value, computed during a run and retained for the next run

// Files used externally
extern char szDataFileName[256];        // The name of the data file.;                   
extern char szScatterFileName[256]; // The name of the result file with the column for DTREE output added
extern char szDTREEFileName[256];   // The name of the DTREE file written/read
extern char szProtocolFileName[256]; // Contains report of program actions; same name as DTREE but with txt extension
extern char szR_or_E_FileName[256]; // Same as original data file, missing values replaced by DTREE or column added
extern char szResultMessage[128];    // Contains a summary of results presented in the GUI
extern FILE *fpOutput; 
extern FILE *fpFileSetupProtocol;    // This protocol file records what happens as the data files are analyzed (can be examined during training by Notepad)

// Parameters used externally
extern int nPercentForTest;   // the percentage of the data file used for testing, between 0 and 50 percent
extern BOOL bEstimateNoiseVariance; // This is TRUE when the TVfile is split into two parts to determine the noise variance in the data. 
extern int nMaxLag; // The maximum lag of any input determined in preprocUniversalFile.
extern double dblFracTest;    // The fraction of the PreprocessedDataFile used for TSfile (default 10%) if no separate test file is provided.
extern CDataFile OutputData;  // The result of evaluation with a column added for the DTREE output
extern int nMessageNumber; // messages used in the Doc and View
extern int nPercentProgress; // used for progress indicator
