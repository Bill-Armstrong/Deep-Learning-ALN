// ALN Library (libaln)
// Copyright (C) 1995 - 2010 William W. Armstrong
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
	// including setting up the necessary files
int ALNAPI analyzeinputfile(char * szDataFileName, int * nHeaderLines, long * nrows, int * ncols, BOOL bPrint);
	// Analyzes the given input data file 
void ALNAPI preprocUniversalFile(); // The universal file is the original dataset and can be preprocessed
void ALNAPI createTVTSfiles();	// The PreprocessedDataFile is used to create TV (training only) and TS files, which are written out..
void ALNAPI analyzeTV();	// Computes the standard deviations of the variables in the TVset.
void ALNAPI getTVfile();	// The TVfile created from the PreprocessedDataFile is read in
void ALNAPI getTSfile();	// The TSfile created from the PreprocessedDataFile is read in 
void ALNAPI createTrainVarianceFiles(int nChooseTR);	// creates training and variance files (several times for bagging)
void ALNAPI dolinearregression();	// This does a truncated linear regression fit to get an upper bound on noise
void ALNAPI approximate();	// This creates the final approximant using the weight bounds found above
void ALNAPI reportFunctions();	// reports on the trained function ALNs with stats and plots
void ALNAPI evaluate();	// Evaluate an existing DTREE on the data file after preprocessing
void ALNAPI cleanup();	// destroys allocated items no longer needed
void ALNAPI outputtrainingresults();	// outputs the results of training
//void ALNAPI validate(CMyAln * pALN);// computes the error of the ALN on samples not used in training
void ALNAPI overtrain(CMyAln * pOT);	// overtrains a single ALN to help obtain an estimate of the level of noise in the data
void ALNAPI trainaverage();	// trains an average of several smoothed ALNs (i.e. with fillets)
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
extern void fillvector(double *, CMyAln*); // inputs a training data vector, perhaps by averaging several ALNs
	// this callback routine allows users to get reports of results during training
extern CDataFile AuxALNinputFile;
extern CDataFile NumericalTestFile;  // These CDataFiles are for preprocessing the test and variance files
extern long nRowsNumericalTestFile;
extern int nColsNumericalTestFile;
extern CDataFile TSfile;
extern long nRowsTSfile;
extern CDataFile NumericalValFile;
extern long nRowsNumericalValFile;
extern int nColsNumericalValFile;
extern CDataFile VarianceFile;
extern long  nRowsVAR;
extern double dblSetTolerance;
extern BOOL bEstimateRMSError; // if TRUE we use a onealnfit to estimate RMS error which needs a variance set
extern double  dblVarianceErr;    // Set equal to the rmse in the variance step
extern double  dblSmoothingFraction; // When smoothing is done at approximation and later, this is the fraction of tolerance used

// variables used in the Pro version and viewed/set in a dialog
// nDim is always the number of input variables plus one for the output
// this is true even when an output column is missing from the data file.
extern int nColsUniv;         // number of columns in the universal file
extern long nRowsUniv;    // number of rows in the universal file
extern int nHeaderLinesUniv;  // the number of headerlines found in the universal file
extern char hdrline[10][2000]; // up to nine header lines plus a buffer
extern char varname[101][32]; //space for datafile column names up to 31 characters long, NULL terminated
extern int nheaderlineitems[10];  // this counts the number of fields read in each header line
extern int nheaderlinewithvars;		//  this is the one of the last two headerlines with the variable names
																	//  however if neither has the correct number of strings this is 10
extern int nALNinputs;	// this is the number of input variables to the ALN, including the desired output
												// the input variables of the ALNs are characterized by column numbers in the data file and associated delays
extern int  nInputCol[101];	// the indeces of datafile columns serving as inputs to the ALNs (columns may be repeated)
extern int  nLag[101]; // the number of rows back from the row of the outputthat a column is sampled
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
extern FILE *fpFileSetupProtocol;    // This protocol file records what happens as the data files are analyzed
// Parameters used externally
extern int nPercentForTest;   // the percentage of the data file used for testing, between 0 and 50 percent
extern BOOL bEstimateRMSError; // This is TRUE when the data file is split into several parts to determine
                               // the RMS noise in the data.  First comes a test set, then a training set,
                               // then a variance set of almost equal size to the test set.  
extern int nMaxLag; // The maximum lag of any input is nMaxLag determined in
                    // preprocUniversalFile. If test, training and variance sets are
                    // specified, then the initial parts of each component must be
                    // removed before use.  
extern double dblFracTest;    // The fraction of the PreprocessedDataFile used for TSfile (default 10%) if no separate test file.
extern CDataFile OutputData;  // The result of evaluation with a column added for the DTREE output
extern int nMessageNumber;
extern int nPercentProgress;
