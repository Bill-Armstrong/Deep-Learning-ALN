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

// Description of the program
//
// ALNfit Deep for MS-Windows is a data analysis program for predicting one column
// in a given data file from values in other columns.  The result is expressed as
// an ALN (adaptive logic network). The logic is not Boolean Logic, but rather
// so-called fuzzy logic using minimum and maximum operators instead of AND and OR.
//
// An ALN computes a piecewise linear continuous function (at least in theory)formed by
// using affine (i.e. non-homogeneous linear) functions in the input layer and two-input
// maximum and minimum operators in all other layers of a neural network.
// The single output is at the root.  All values are real (double precision floating
// point in the computer representation).
//
// The inputs of the ALN are taken from the specified columns of the data file. A column,
// including the output column, may furnish several values with different lags
// (a lag of n means that the data is taken from a row
// n rows closer to the first row of the data file than the output row).
// We don't allow the output column to be used as an input value with lag 0,
// since that is the desired output value.
//
// Caution: The above description of an ALN is not the same as was used in the literature
// prior to 1995. There is a strong connection. If you take one of the above ALNs, representing a function f, and
// ask what is the logical value of the comparison  f >= T for some fixed threshold value T.
// You can use the following transformations to transform the ALN into a pre-1995 ALN with ANDs and ORs
// and Boolean-valued inputs computed by Perceptrons (a Perceptron compares an affine function to a constant):
// MIN(g,h) >= T iff g >= T AND h >= T; MAX(g,h) >= T iff g >= T OR h >= T. Changing to
// ALNs which processed real values instead of Booleans was a significant advance
// because the training algorithm had an obvious way to solve the "credit assignment"
// problem. Namely at a minimum node, the input subtree with the lesser of the two input
// values is the one that should get credit for a good or bad output. Except in
// the rare case of ties, this leads to a single linear function whose weights should adapt,
// and the way to do that is well-known.  One can view today's ALNs as using fuzzy logic,
// hence the name has been kept from the Boolean case.  The linear functions adapt by an iterative,
// least-squares approach, and the tree of MAX and MIN operators adapts by growing to better fit the data.
//
// Classes in classification problems are always represented by 0 and 1 values
// Problems with more than two classes are treated like regression problems with integer outputs,
// or must be reformulated as a sequence of two-class problems.
// Input data files must be tab-separated text files.
// The input data files must have no missing or extra fields,
// and all fields must be numeric (integers or decimal numbers)
// using only decimal points (.) or, optionally, commas (,) to indicate fractional parts.
// Undefined inputs and outputs are represented by 99999 or -9999.
// With undefined values represented in this way, the output "R" file replaces those values.
// For evaluation files, the desired output field may be missing in all rows.

// Training and testing mode
//
// The user determines that either training or evaluation of a file using a
// trained ALN will be done.
// For training, the program normally uses the input data file and creates 
// from it one file for training/variance ("TVfile.txt") and
// one file ("TSfile.txt") for testing. The latter file is a specified
// part of the data (eg a randomly selected file consisting of 10% of the rows of
// the input data file). The user can alternatively specify a testing file directly
// which is evaluated using a pre-existing DTREE. The testing file is 
// held back for making unbiased tests after the result of training is ready.
// The result of training is a single DTREE file of variable depth,
// and this is written out as a text file with extension .dtr. 
// The DTREE file can be viewed with a text-file viewer like MS Notepad.
// It can be reloaded and parsed by programs in the DTREE library part of libaln.
// The DTREE can also be written out in binary form (for faster loading).

// Details of the training procedure
//
// Unlike many other machine learning programs, ALNfit Deep depends on
// obtaining an estimate of the RMS error in the training data.
// Of course, looking at any set of data points where no input vector is
// repeated, we can interpolate that data as accurately as we wish and say that the
// error is 0.  However if we also look at a second  set from the same
// data source, then it should define the same function. The difference of functions
// defined by the two data sets can be used to measure the noise.
// Linear regression is first used to fit the training set to determine its mean square error.
// A good ALN fit should have much smaller error, since an ALN can use more than one linear piece.
// Then we try to fit the training data very closely, which means something like interpolation. This amounts to
// what is called "overtraining", and leads to poor "generalization". However
// we need something like overtraining so values of the function can be compared to 
// values of nearby samples to estimate noise variance in a subregion.
// We arrive at a "noise variance function", that smooths noise variance samples such that subsequently
// when we train ALNs on the same data for a good fit, we stop growing the tree (splitting linear
// pieces into two) when the error is below the level of noise variance. 
//
// We assume that the RMS noise added to the ideal function to get samples is slowly varying oner the input space.
// This gives the possibility to learn fine detail in the unknown function when the noise variance is small.
//
//  If the noise in the samples varies proportionately to the value of the ideal function,
//  then we should probably do the training on the logarithm of the output data points to get constant noise variance. 
//
// The approximation step
//
// Once we know how much RMS error there is in the training data, we can train one or more ALNs.
// We can stop elaborating an ALN by splitting linear pieces into two, when the training square error
// becomes less than the noise variance because to continue would just be fitting the noise in the data.
// The smoothing constant used for fillets is set proportionally.
// Several trainings using the noise variance are done and the results averaged.

// The TVfile is used several (nALNs) times to create several ALNs whose average
// will have good generalization performance.  The noise variance function for averaging is set at the
// level found for training one ALN divided by the number of ALNs in the average.
// This is called bagging (inventor: Leo Breimann).


// Training ALNs
//
// Normally several ALNs, eg seven, are trained to get the final approximant.
// Each ALN learns a piecewise linear approximant ( perhaps with additional smoothing using  
// quadratic fillets) giving the output as a smooth function of the inputs.
// During training, each linear piece changes coefficients (weights) to fit the data
// points for which it is responsible. The fit takes account of the fillets.
// If the error of a piece is still above the level of noise variance
// after enough training has been given to move it into position,
// the piece splits into two, initially equal to it, and those two linear pieces then adapt
// their weights to fit their respective data points.  Splitting continues as required
// to get the RMS training error below the output error tolerance on all pieces.
// Each ALN is trained on the entire TVset. This is a departure from standard practise
// which is being tested.

// We also remove some (10%, or at least 10) samples for testing purposes
// (they must never be used for training, just as the
// training samples must never be used for testing if you want a good fit on as
// yet unseen samples. 

// If you have a low dimensional function, you can ask yourself looking at a graph how many
// linear pieces will be required to get an acceptable fit.  If that number is N, and
// there are C columns in the data including the output column, then the number of samples
// (ie rows) in the file must be at least S = 2 * N * C. You can think of it this way:
// Of the S rows, you take away 10% leaving 0.9 * S. That means for each linear piece, you have 1.8 * C data points.  The minimum you
// need for defining a linear piece is C points.  If there is a lot of noise in the
// data, you might need many times this number of rows. Using machine learning requires
// plenty of samples.


// Analysis of the ALNs on the TV set
//
// The average absolute weights over the TV set are computed and scaled by
// multiplying by the standard deviation of the input concerned and divided by
// the standard deviation of the output.

// The final approximant
//
// Taking the average of several ALNs gives a smooth approximant, which is good.  On the downside
// it may take a relatively long time to compute the output value for a given input vector.
// The computation can be greatly speeded up by a process referred to as resampling,
// The average of the ALNs computed at many input points, and those new 
// artificially generated "samples" are used to train a single ALN at a lesser
// noise variance before.   
// The average ALN is converted to a DTREE for speed. The DTREE is strictly piecewise
// linear (ie it has no quadratic or quartic fillets).  It is also smooth because it is fitting
// a smooth function.  It is also very fast, particularly if one splits the input space
// into blocks where only a few linear pieces are required to compute the output.
// NB currently DTREEs are restricted to one level in this program, ie there is no division of
// the input space into pieces.

// Report of the results of the trial:
//
// A report is generated that gives the RMS error of the DTREE approximant on
// the test set. An output file is generated with the original data file with
// an additional column at the right containing the DTREE prediction of the output.
// This is the "E" file, and rightmost column can be compared to others in a spreadsheet.
// If the data file is missing the output during evaluations without training,
// then a zero column takes the place of the missing output column and the DTREE
// output column is to the right.
// A protocol of useful information is written to a file to document all steps and help debugging
#ifdef ALNDLL
#define ALNIMP __declspec(dllexport)
#endif

#include <aln.h>
#include "alnpriv.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// include files
#include <stdafx.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <alnpp.h>
#include <dtree.h>
#include <datafile.h>
#include ".\cmyaln.h" 
#include "alnextern.h"
#include "alnintern.h"


#define ALNAPI __stdcall

// functions used externally
void ALNAPI ALNfitSetup();
int ALNAPI analyzeInputFile(char * szDataFileName, int * nHeaderLines, long * nrows, int * ncols, BOOL bPrint); // Analyzes the given input data file 
void ALNAPI preprocessDataFile(); // creates UNfile, the original dataset without header, with a column added at right
void ALNAPI createTVTSfiles();      // The PreprocessedDataFile is used to create TV (training only) and TS files, which are written out..
void ALNAPI analyzeTV();            // Computes the standard deviations of the variables in the TVset.
void ALNAPI getTVfile();          // The TVfile created from the PreprocessedDataFile is read in
void ALNAPI getTSfile();          // The TSfile created from the PreprocessedDataFile is read in 
void ALNAPI reportFunctions();      // reports on the trained function ALNs with stats and plots
void fillvector(double *, CMyAln*); // Used to input a data vector under program control

// The main steps in training operations (train_ops.cpp)
void ALNAPI doLinearRegression();   // This does a linear regression fit, finding RMS error and weights
void ALNAPI computeNoiseVariance();  // This makes noise variance samples. These are then used to train an ALN so samples are smoothed.
void ALNAPI approximate();          // This creates the final approximant using the weight bounds found above
void ALNAPI trainAverage();         // Does bagging by averaging several ALNs created using noise variance stopping
void ALNAPI outputTrainingResults();// Prints out the results of training 
void ALNAPI constructDTREE(int);		// The DTREE is a VERY fast way of evaluating any ALN.
void ALNAPI evaluate();			        // Evaluate an existing DTREE on the data file after preprocessing
void ALNAPI cleanup();              // Destroys some objects previously allocated

// global variables used externally
BOOL bClassify = FALSE;      // This is TRUE if the user chose a classification problem, and FALSE for regression
BOOL bRegress = TRUE;        // This is FALSE iff the TV file outputs are all 0's and 1's
BOOL bTrain = TRUE;          // This is TRUE for training and FALSE for evaluation.
BOOL bJitter = FALSE;         // This is re-set under the options dialog
BOOL bReplaceUndefined = FALSE; // Makes a new column at right if FALSE, otherwise replaces undefined output values
BOOL bDiagnostics = FALSE;  // For controlling printout of diagnostic files
BOOL bTimePrefixes = TRUE;
BOOL bPrint = TRUE;  // controls printing of the input files, may be changed in options
int nALNs = 1; // since we train on all of the TVfile now, there is no point to doing a lot of bagging
int nMessageNumber =8;
int nPercentProgress = 0;
int nDTREEDepth = 1;
double dblEvalRMSError = -1;
int nEvalMisclassifications = -1;
double dblEvalMisclassificationPercent = -1;
char hdrline[10][2000]; // global for the data file, and is declared extern in alnextern.h
char varname[101][32]; //space for datafile column names up to 31 characters long, NULL terminated
int nColsUniv;         // number of columns in the universal file
long nRowsUniv = 0;    // number of rows in the universal file
int nHeaderLinesUniv;  // the number of headerlines found in the universal file
int nheaderlineitems[10];  // this counts the number of fields read in each header line
int nheaderlinewithvars; //  this is the one of the last two headerlines with the variable names
int nALNinputs; // this is the number of input variables to the ALN, including the desired output
int  nInputCol[101];   // the indeces of datafile columns serving as inputs to the ALNs (columns may be repeated)
int  nLag[101]; // the number of rows back from the row of the outputthat a column is sampled
double dblMax = 1e12;  // maximum number allowed for weight
double dblMaxWeight[101]; // the maximum weight a given input can have
double dblMinWeight[101]; // the minimum weight a given input can have
double dblImportance[101]; // the importance value, computed during a run and retained for the next run

// Files used externally
char szDataFileName[256];        // The name of the data file.;                   
char szScatterFileName[256]; // The name of the result file with the column for DTREE output added
char szDTREEFileName[256];   // The name of the DTREE file written/read
char szProtocolFileName[256]; // Contains report of program actions; same name as DTREE but with txt extension
char szR_or_E_FileName[256]; // Same as original data file, missing values replaced by DTREE or column added
char szResultMessage[128];    // Contains a summary of results presented in the GUI

// Parameters used externally
int nPercentForTest = 10;     // the percentage of the data file used for testing, between 0 and 50 percent
int nMaxLag;                  // The maximum lag of any input is nMaxLag determined in preprocUniversalFile. 
double dblFracTest = 0.1;     // The fraction of the PreprocessedDataFile used for TSfile (default 10%) if no separate test file.

//FILE *fpOutput = NULL;             // the output data file resulting from evaluation
double dblSetTolerance;

// from alnintern.h


// Functions used only internally
//void splitControl(ALN*, double); // if average variance error of a piece is high, splitting is prevented
//void doSplits(ALN*, ALNNODE*, double); //does the recursion of splitcontrol
//void zeroSplitValues(ALN*, ALNNODE*);  // sets the square error to zero in each LFN
//void splitUpdateValues(ALN*); // accumulates the training square error and number of hits on each linear piece
//void dodivideTR(ALN*, ALNNODE*); // divides the total square training set errors of the pieces by their hit count
//void dodivideVAR(ALN*, ALNNODE*); // divides the total square variance errors of the pieces by their hit count

// thread procedures
UINT TakeActionProc(LPVOID pParam);  // separate thread


//global variables used only internally
BOOL bTrainingAverage = FALSE;// Switch to tell fillvector whether get a training vector or compute an average
int nNumberLFNs;  // used to control the epoch size, which should be proportional to this
char szVarName[100][3];
BOOL bEstimateNoiseVariance = TRUE; // if TRUE we estimate noise variance
BOOL bDecimal = TRUE; // means numbers could have a decimal point
BOOL bComma = TRUE;   // means numbers could have a comma
//double dblSetTolerance; // Value of tolerance set in the options dialog.
// It is not used now. It would be the mean of the noise variance function.
int nDim = 0;
int nOutputIndex = 0;
long nRowsPP   = 0;          // The number of rows in the PreprocessedDataFile
long nRowsTV   = 0;          // The number of rows in the TVfile (Training & Variance File)
double* adblEpsilon;         // An array of doubles holding the tolerances of the inputs.
double* adblMinVar;          // Array of minima of the variables
double* adblMaxVar;          // Array of maxima of the variables
double* adblStdevVar;        // Standard deviations of the variables
double  dblTrainErr;         // Set in cmyaln.h at the end of training
double  dblVarianceErr;    // Set equal to the rmse in the variance step
double  dblLinRegErr;        // The error of linear regression for use in upper-bounding output tolerance




// Files used only internally
FILE *fpFileSetupProtocol = NULL;    // This protocol file records what happens as the data files are analyzed
FILE *fpData = NULL;                 // The data file which contains all data.
FILE *fpProtocol = NULL;             // the file to record the results of the experiment. Rename to save.
FILE *fpOutput = NULL;             // the output data file resulting from evaluation
FILE *fpReplacement;
int nColsNumericalTestFile;
long nRowsTS;
long nRowsNumericalValFile;
int nColsNumericalValFile;
CDataFile UNfile;             // copy of the data file, but with missing values replaced by special number
CDataFile PreprocessedDataFile;             // The preprocessed file created from the Universal file
CDataFile NumericalTestFile;  // These CDataFiles are for preprocessing the test and variance files
CDataFile NumericalValFile;   // The original data file has headers and comments. These are removed in this file.
CDataFile TSfile;							// Samples held back for testing
CDataFile TVfile;             // The file used for training and, if no separate file is given, for variance, with nDim columns and nRowsUniv - nRowsTS rows.
CDataFile OutputData;         // The result of evaluation with a column added for the DTREE output

// the variables in a different way when a comparison between the average errors on training and noise variance sets
// is needed to decide whether or not to allow splitting of a linear piece
typedef struct tagSPLIT      // Used in inhibiting splitting -- must be zeroed before and after use
{
  int nCount;                // Number of hits
  double dblSqErrorTrain;    // Squared error of a piece after training                      
  double dblRespTotal;      // Used during training and noise variance estimation
  double dblT_NotUsed;			 // Used in the SDK only
} SPLIT;

// now actually define the storage for weight and centroid values from linear regression
double* adblLRW; // weight components
double* adblLRC; // centroid components

using namespace std;

void ALNAPI ALNfitSetup() // routine
{
  // we can now close the FileSetupProtocol file
	if(bPrint) // if !bPrint it was already closed in ALNfitDeepView.cpp (l. 1526)
  {
    fclose(fpFileSetupProtocol);
  }
  // open the Protocol file which keeps a protocol of all operations
	if ((fpProtocol=fopen(szProtocolFileName,"w")) != NULL)
  {
    fprintf(fpProtocol,"***********************************************************************\n");
    fprintf(fpProtocol,"ALNfit Deep Automatic Regression and Classification Program\n");
    fprintf(fpProtocol,"ALNfit Deep is an open-source program -- see LGPL license in source files\n");
    fprintf(fpProtocol,"********************** ALNfit setting up the run **********************\n");
  }
  else
  {
    exit(0);
  }  
  // print the program header info
   	//seed the random number generator
	CAln::SRand((unsigned int) time(NULL));
	// print out whether we are training or not and classifying or doing regression
	if(bClassify)
	{
		fprintf(fpProtocol,"The task is classification\n");
	}
	else
	{
		fprintf(fpProtocol,"The task is regression\n");
	}
	if(bTrain)
	{
		fprintf(fpProtocol,"A DTREE is to be created by training using the data file\n");
	}
	else
	{
		fprintf(fpProtocol,"An existing DTREE will be used to evaluate the data file\n");
	}

	// ************** ANALYZE THE INPUT DATA FILE ***********************
  int nheaderlines, ncols;
  long nrows;
	nheaderlines = 0;
  nrows = 0;
  ncols = 0;
 	if(!analyzeInputFile(szDataFileName, &nheaderlines, &nrows, &ncols,TRUE))
	{
    fprintf(fpProtocol,"Input file analysis failed\n");
    fflush(fpProtocol);
		exit(0);
	}
	if(nheaderlines > 9)
	{
    fprintf(fpProtocol,"Stopping: You can't have more than nine headerlines in your data file!\n");
    fflush(fpProtocol);
		exit(0);
	}
  nHeaderLinesUniv = nheaderlines;
	nColsUniv = ncols;
	nRowsUniv = nrows; // this means just data rows, not including any header line
  fprintf(fpProtocol,"The number of columns in the data file is %d\n", nColsUniv); 
  fprintf(fpProtocol,"The number of header rows in the data file is %d\n", nHeaderLinesUniv); 
  fprintf(fpProtocol,"The number of data rows in the data file is %d\n", nRowsUniv); 
	/* At this point, we have the number of  rows and columns in the UNfile, and
	the number of headerlines in the data file that have to be eliminated or used
	as column names. n*/

	// ***************** SET UP THE FILES *************************
	// create the necessary files: UNfile PreprocessedDataFile, TVfile, TSfile
  // the UN file just copies the input file,EXCLUDING the header(if present)
  UNfile.Create(nRowsUniv,nColsUniv + 1); // we add an extra column for an evaluated output
	for(int k = 0; k < nColsUniv; k++)
	{
		for(int nrow = 0; nrow < nRowsUniv; nrow++)
		{
			UNfile.SetAt(nrow,k,0,0);
		}
	}
	// initialize the extra column to 0
	for(int nrow = 0; nrow < nRowsUniv; nrow++)
  {
	  UNfile.SetAt(nrow,nColsUniv,0,0);
  }
	fprintf(fpProtocol,"Preprocessing splits data file into TVfile and TSfile (for testing). \n");
	preprocessDataFile();
	createTVTSfiles();
}

int ALNAPI analyzeInputFile(char * szDataFileName, int * pheaderlines, long * prows, int * pcols, BOOL bPrint) // routine
// input: szDataFileName file, outputs: number of header lines, variable names, number of data rows and columns
// this routine accepts files with up to 101 columns
{
	long linecount = 0;
	int  ncount = 0, itemcount = 0, nheaderlines = 0, ncols = 0;
	char item[101][128]; // for reading in string items (e.g. up t0 101 variable names)
	float fitem[101];    // for reading in data lines containing floats
	for (itemcount = 0; itemcount < 101; itemcount++)
	{
		fitem[itemcount] = 0; // zero the array for floats in datalines
	}
	//alnextern.h: static char hdrline[10][2000]; // this will store the whole header
	// max line length in data file is 2000 bytes including null terminator
	ifstream ifs(szDataFileName);
	if (ifs)
	{
		if (bPrint)fprintf(fpFileSetupProtocol, "Opening data file %s for analysis succeeded!\n", szDataFileName);
		if (bPrint)fflush(fpFileSetupProtocol);
	}
	else
	{
		if (bPrint)fprintf(fpFileSetupProtocol, "Stopping. Opening data file %s for analysis failed!\n", szDataFileName);
		if (bPrint)fflush(fpFileSetupProtocol);
		exit(0);
	}
	ncount = 2000; // maximum input line length
	nheaderlines = 0; // current header line number
	linecount = 0; // current data row number, which finally indicates the number of data rows
	itemcount = 0; // current item in the row
	ncols = 0;
	BOOL bStillHeader = TRUE; // we allow empty lines while still in the header, so we need this flag
	// a line of the hdrline array is used as a buffer (it can store up to 10 lines)
	while (ifs.getline(hdrline[nheaderlines], ncount, '\n') && (ncount > 0)) // read one line
	{
		// comments in input files can be indicated by // at the start of the line
		// comment lines are ignored except when output is generated from the input file
		if ((hdrline[nheaderlines][0] == '/') && (hdrline[nheaderlines][1] == '/')) continue;
		// we allow empty lines in the header area, not in the data rows
		if (!bStillHeader && hdrline[nheaderlines][0] == '\0') continue; // if getline reads an empty line in the data part, the while is terminated
		// the above are either comments or blank header lines
		linecount++; // assume it's a data line, but be prepared to reverse this decision and make it a header line
		itemcount = 0;
		if (bStillHeader && (hdrline[nheaderlines][0] == '\0'))
		{
			if (bPrint)fprintf(fpFileSetupProtocol, "The line read from the file is empty\n");
			if (bPrint)fprintf(fpFileSetupProtocol, "The line is assumed to be a header line. Analysis continues.\n");
			if (bPrint)fflush(fpFileSetupProtocol);
			nheaderlineitems[nheaderlines] = 0;
			nheaderlines++;
			linecount--;  // reverse the above increment
			bDecimal = TRUE;
			bComma = TRUE;
			continue;
		}
		// at this point, we know the line is neither a comment, nor empty
		// look at the periods, commas and numbers 0-9 in the line 
		int charcount = 0;
		while (hdrline[nheaderlines][charcount] != '\0') // process up to the null terminator
		{
			if ((hdrline[nheaderlines][charcount] != ' ') && (hdrline[nheaderlines][charcount] != '\t'))
			{
				if (hdrline[nheaderlines][charcount] == '.')
				{
					bComma = FALSE; // this can't be a number line with commas
				}
				else if (hdrline[nheaderlines][charcount] == ',')
				{
					bDecimal = FALSE; // this can't be a number line with decimal points
				}
				else if ((hdrline[nheaderlines][charcount] < '0' || hdrline[nheaderlines][charcount] > '9') &&
					(hdrline[nheaderlines][charcount] != '-') && (hdrline[nheaderlines][charcount] != 'E')
					&& (hdrline[nheaderlines][charcount] != 'e'))
				{
					// in this case, the line does not consist of numbers only
					// numbers like 1,234,567.89 and 1.234.567,89 are not allowed in data input
					// and this is treated like a header line
					bComma = FALSE;
					bDecimal = FALSE;
				}
			}
			charcount++;
		}
		if (bComma == FALSE && bDecimal == FALSE)
		{
			// this is a header line
			if (bPrint)fprintf(fpFileSetupProtocol, "The line read from the file is: \n %s \n", hdrline[nheaderlines]);
			if (bPrint)fprintf(fpFileSetupProtocol, "The line is assumed to be a header line. Analysis continues.\n");
			if (bPrint)fflush(fpFileSetupProtocol);
			itemcount = sscanf(hdrline[nheaderlines], "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
				item + 0, item + 1, item + 2, item + 3, item + 4, item + 5, item + 6, item + 7, item + 8, item + 9,
				item + 10, item + 11, item + 12, item + 13, item + 14, item + 15, item + 16, item + 17, item + 18, item + 19,
				item + 20, item + 21, item + 22, item + 23, item + 24, item + 25, item + 26, item + 27, item + 28, item + 29,
				item + 30, item + 31, item + 32, item + 33, item + 34, item + 35, item + 36, item + 37, item + 38, item + 39,
				item + 40, item + 41, item + 42, item + 43, item + 44, item + 45, item + 46, item + 47, item + 48, item + 49,
				item + 50, item + 51, item + 52, item + 53, item + 54, item + 55, item + 56, item + 57, item + 58, item + 59,
				item + 60, item + 61, item + 62, item + 63, item + 64, item + 65, item + 66, item + 67, item + 68, item + 69,
				item + 70, item + 71, item + 72, item + 73, item + 74, item + 75, item + 76, item + 77, item + 78, item + 79,
				item + 80, item + 81, item + 82, item + 83, item + 84, item + 85, item + 86, item + 87, item + 88, item + 89,
				item + 90, item + 91, item + 92, item + 93, item + 94, item + 95, item + 96, item + 97, item + 98, item + 99, item + 100);
			if ((itemcount == EOF) || (nheaderlines > 8))
			{
				if (bPrint)fprintf(fpFileSetupProtocol, "Too many (>9) headerlines or EOF encountered reading the header line\n");
				if (bPrint)fflush(fpFileSetupProtocol);
				exit(0);
			}
			else
			{
				nheaderlineitems[nheaderlines] = itemcount;
				if (bPrint)fprintf(fpFileSetupProtocol, "There were %d items read in the header line\n", itemcount);
				if (bPrint)fflush(fpFileSetupProtocol);
			}
			nheaderlines++;
			linecount--;
			bDecimal = TRUE;
			bComma = TRUE;
			continue;
		}
		// when we reach this point, we should have seen the last of the header lines
		ASSERT(linecount >= 1);
		bStillHeader = FALSE;
		// hdrline[nheaderlines] is used as a buffer below
		if (bDecimal == FALSE) // but bComma is still TRUE because it's not a header line
		{
			charcount = 0;
			while (hdrline[nheaderlines][charcount] != '\0')
			{
				if (hdrline[nheaderlines][charcount] == ',')
				{
					hdrline[nheaderlines][charcount] = '.'; // convert the number to decimal points
				}
				charcount++;
			}
		}
		// now read the numbers in the above line
		itemcount = sscanf(hdrline[nheaderlines], "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
			fitem + 0, fitem + 1, fitem + 2, fitem + 3, fitem + 4, fitem + 5, fitem + 6, fitem + 7, fitem + 8, fitem + 9,
			fitem + 10, fitem + 11, fitem + 12, fitem + 13, fitem + 14, fitem + 15, fitem + 16, fitem + 17, fitem + 18, fitem + 19,
			fitem + 20, fitem + 21, fitem + 22, fitem + 23, fitem + 24, fitem + 25, fitem + 26, fitem + 27, fitem + 28, fitem + 29,
			fitem + 30, fitem + 31, fitem + 32, fitem + 33, fitem + 34, fitem + 35, fitem + 36, fitem + 37, fitem + 38, fitem + 39,
			fitem + 40, fitem + 41, fitem + 42, fitem + 43, fitem + 44, fitem + 45, fitem + 46, fitem + 47, fitem + 48, fitem + 49,
			fitem + 50, fitem + 51, fitem + 52, fitem + 53, fitem + 54, fitem + 55, fitem + 56, fitem + 57, fitem + 58, fitem + 59,
			fitem + 60, fitem + 61, fitem + 62, fitem + 63, fitem + 64, fitem + 65, fitem + 66, fitem + 67, fitem + 68, fitem + 69,
			fitem + 70, fitem + 71, fitem + 72, fitem + 73, fitem + 74, fitem + 75, fitem + 76, fitem + 77, fitem + 78, fitem + 79,
			fitem + 80, fitem + 81, fitem + 82, fitem + 83, fitem + 84, fitem + 85, fitem + 86, fitem + 87, fitem + 88, fitem + 89,
			fitem + 90, fitem + 91, fitem + 92, fitem + 93, fitem + 94, fitem + 95, fitem + 96, fitem + 97, fitem + 98, fitem + 99, fitem + 100);

		if (itemcount > 101)
		{
			if (bPrint)fprintf(fpFileSetupProtocol, "Error: too many input columns in this line.\n");
			if (bPrint)fprintf(fpFileSetupProtocol, "Maximum is 100 ALN inputs (plus one more for the desired output) for this program.\n");
			if (bPrint)fflush(fpFileSetupProtocol);
			return 0;
		}
		if ((itemcount == 1) && (linecount == 1))
		{
			if (bPrint)fprintf(fpFileSetupProtocol, "The line read from the file is: \n %s \n", hdrline[nheaderlines]);
			if (bPrint)fprintf(fpFileSetupProtocol, "You have only one column, appropriate for a single time series.\n");
			if (bPrint)fprintf(fpFileSetupProtocol, "You have to use delayed values taken from this column to train.\n");
			if (bPrint)fflush(fpFileSetupProtocol);
		}
		if (linecount == 1)
		{
			if (bPrint)fprintf(fpFileSetupProtocol, "First data line to help checking for correctness:\n%s\n", hdrline[nheaderlines]);
			ncols = itemcount; // the number of columns is set in the first data line
			if (bPrint)fprintf(fpFileSetupProtocol, "Count of items in first data line = %d\n", itemcount);
			if (bPrint)fflush(fpFileSetupProtocol);
		}
		else
		{
			if (itemcount != ncols)
			{
				if (bPrint)fprintf(fpFileSetupProtocol, "Inconsistent number of items in data line %d ", linecount);
				if (bPrint)fprintf(fpFileSetupProtocol, "namely %d \n", itemcount);
				if (bPrint)fprintf(fpFileSetupProtocol, "Please correct and make sure each input data file line \n");
				if (bPrint)fprintf(fpFileSetupProtocol, "has the same number of items in it.\n");
				if (bPrint)fflush(fpFileSetupProtocol);
				return 0;
			}
		}
	} // end of while (ifs.getline(hdrline[nheaderlines], ncount, '\n') && (ncount > 0))
	if (bDecimal == FALSE && bComma == FALSE)
	{
		if (bPrint)fprintf(fpFileSetupProtocol, "Stopping: It appears that both commas and decimal points are used in number lines");
		if (bPrint)fflush(fpFileSetupProtocol);
		exit(0);
	}
	*pheaderlines = nheaderlines; // transmit to the caller
	*prows = linecount;
	*pcols = ncols;
	if (bPrint)fflush(fpFileSetupProtocol);
	return 1;
}

void ALNAPI preprocessDataFile()  // routine input: szDataFileName file outputs: UNfile, PreprocessedDataFile
{
	if ((fpData = fopen(szDataFileName, "r")) != NULL)
	{
		fprintf(fpProtocol, "Opening file %s succeeded.\n", szDataFileName);
	}
	else
	{
		fprintf(fpProtocol, "Opening file %s failed. Stop.\n", szDataFileName);
	}
	// recode all the entries to a different array of doubles
	// both original and recoded values for the row are declared here
	char  Dummyline[2000];
	float ALNvar = 0; // limited to 100 input variables plus one output variable
	// copy the entire data file into UNfile EXCEPT the header
	// NB datafiles can only hold doubles!
	for (int i = 0; i < nHeaderLinesUniv; i++)
	{
		// ignore it if it is a comment line and get the next one
		// we are not trying to accumulate headerlines here, just get to the last non-comment header line
		fgets(Dummyline, 2000, fpData);
		while ((Dummyline[0] == '/') && (Dummyline[1] == '/')) fgets(Dummyline, 2000, fpData);
	}
	// we have all the header lines, but there may be other comment lines
	// those can be captured below
	for (int i = 0; i < nRowsUniv; i++)
	{
		fgets(Dummyline, 2000, fpData);
		// ignore any comment lines
		while ((Dummyline[0] == '/') && (Dummyline[1] == '/')) fgets(Dummyline, 2000, fpData);
		// when we get a non-comment line, we keep it, assuming it is data
		// we convert the commas to decimal points if necessary so we can scan as a number
		if (bDecimal == FALSE) // but bComma is still TRUE because it's not a header line
		{
			int charcount = 0;
			while (Dummyline[charcount] != 0)
			{
				if (Dummyline[charcount] == ',')
				{
					Dummyline[charcount] = '.'; // change to decimal point
				}
				charcount++;
			}
		}
		double item[101];

		int itemcount = sscanf(Dummyline, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			item + 0, item + 1, item + 2, item + 3, item + 4, item + 5, item + 6, item + 7, item + 8, item + 9,
			item + 10, item + 11, item + 12, item + 13, item + 14, item + 15, item + 16, item + 17, item + 18, item + 19,
			item + 20, item + 21, item + 22, item + 23, item + 24, item + 25, item + 26, item + 27, item + 28, item + 29,
			item + 30, item + 31, item + 32, item + 33, item + 34, item + 35, item + 36, item + 37, item + 38, item + 39,
			item + 40, item + 41, item + 42, item + 43, item + 44, item + 45, item + 46, item + 47, item + 48, item + 49,
			item + 50, item + 51, item + 52, item + 53, item + 54, item + 55, item + 56, item + 57, item + 58, item + 59,
			item + 60, item + 61, item + 62, item + 63, item + 64, item + 65, item + 66, item + 67, item + 68, item + 69,
			item + 70, item + 71, item + 72, item + 73, item + 74, item + 75, item + 76, item + 77, item + 78, item + 79,
			item + 80, item + 81, item + 82, item + 83, item + 84, item + 85, item + 86, item + 87, item + 88, item + 89,
			item + 90, item + 91, item + 92, item + 93, item + 94, item + 95, item + 96, item + 97, item + 98, item + 99, item + 100);
		ASSERT(itemcount == nColsUniv);
		for (int ncoldata = 0; ncoldata < nColsUniv; ncoldata++)
		{
			// the UNfile consists of all the data lines with header and comments removed
			// it may still contain missing values 99999 -9999
			UNfile.SetAt(i, ncoldata, item[ncoldata], 0);
		}
	}
	//if (bPrint && bDiagnostics) UNfile.Write("DiagnoseUNfile.txt");
	//if (bPrint && bDiagnostics) fprintf(fpProtocol, "DiagnoseUNfile.txt written\n");
	// now we create a file with ALN inputs as lines.
	// this has to use several lines of the UNfile because of lags
	// we have to compute the maximum lag
	nMaxLag = 0;
	for (int ninput = 0; ninput < nALNinputs - 1; ninput++) // the output column is 0 and won't be the max
	{
		if (nLag[ninput] > nMaxLag) nMaxLag = nLag[ninput];
	}
	PreprocessedDataFile.Create(nRowsUniv, nALNinputs); // this file may contain some rows with missing values
	// at the end in the reverse order encountered
	// we first zero it out
	for (int k = 0; k < nRowsUniv; k++)
	{
		for (int ninput = 0; ninput < nALNinputs; ninput++)
		{
			PreprocessedDataFile.SetAt(k, ninput, 0, 0);
		}
	}
	double dblValue = 0;
	int k = 0; // this is the row in the file PreprocessedDataFile. If there are undefined items
						 // in row i of the UNfile, k might not advance in synch with i.
	int back = nRowsUniv - 1;
	for (int i = nMaxLag; i < nRowsUniv; i++)
	{
		BOOL bRowGood = TRUE;
		for (int ninput = 0; ninput < nALNinputs; ninput++)
		{
			// THIS IS WHERE THE TRAINING DATA ARE DEFINED
			dblValue = UNfile.GetAt(i - nLag[ninput], nInputCol[ninput]);
			if (dblValue == 99999 || dblValue == -9999) bRowGood = FALSE;
			PreprocessedDataFile.SetAt(k, ninput, dblValue, 0);
		}
		if (bRowGood)
		{
			// moving on to the next row
			k++;
		}
		else
		{
			// case of missing values: erase the bad row in the PreprocessedDataFile,
			// add the row with missing values at the end of the PreprocessedDataFile going backwards WHY??
			for (int ninput = 0; ninput < nALNinputs; ninput++)
			{
				PreprocessedDataFile.SetAt(k, ninput, 0, 0); // erase the line just filled with some missing values
				dblValue = UNfile.GetAt(i - nLag[ninput], nInputCol[ninput]);
				PreprocessedDataFile.SetAt(back, ninput, dblValue, 0); // fill the last available line with the line just erased
			}
			back--;
		}
	}
	ASSERT(back == k + nMaxLag - 1); // k and back have both moved too far!
	// between the good rows and the rows with undefined, there are nMaxLag rows of garbage
	fprintf(fpProtocol, "PreprocessedDataFile created with %d rows\n", k);
	// write it out
	nRowsPP = k; // this is the number of rows in the date with headers, comments, 
	// commas in numbers and missing values removed 
	//if (bPrint && bDiagnostics) PreprocessedDataFile.Write("DiagnosePreprocessedDataFile.txt");
	//if (bPrint && bDiagnostics) fprintf(fpProtocol, "DiagnosePreprocessedDataFile.txt written\n");
	fclose(fpData);
	if (bTrain)
	{
		fprintf(fpProtocol, "Warning: if the ideal function to be learned is complicated,\n");
		fprintf(fpProtocol, "Or if there is a lot of noise in the data, there must be MANY times\n");
		fprintf(fpProtocol, "more rows in the data file as there are ALN inputs!\n");
		nDim = nALNinputs; // nDim = nALNinputs for training since an output column must be present
		nOutputIndex = nDim - 1;

		// allocate arrays in the free store
		adblLRC = new double[nDim]; // C has all nDim centroids that match weights of index one greater except for the last
		adblLRW = new double[nDim + 1]; // W has the output component which is -1, and uses the 0 component to speed things up
		nRowsTS = (int)(dblFracTest * nRowsPP); // use e.g 10% of data file for test ??? have the rows with missing vals been removed?
		if (nRowsPP < 100)
		{
			fprintf(fpProtocol, "Stopping: There must be at least 100 rows in the data file\n");
			exit(0);
		}
		if ((nRowsPP - nRowsTS) <= nALNinputs)
		{
			fprintf(fpProtocol, "Stopping: There must be at least as many rows in the data file (after removing a given percentage for testing) as there are ALN inputs.\n");
			exit(0);
		}
		nRowsTV = nRowsPP - nRowsTS;
	}
	else // we are evaluating the PreprocessedDataFile
	{
		nRowsTS = nRowsPP;
		nRowsTV = 0;
		// the number of ALNs doesn't matter for evaluation since only a DTREE
		// is used for evaluation
	}
	fprintf(fpProtocol, "The number of rows in the test set is %d\n", nRowsTS);
	if (bTrain) // nRowsTV can't be zero here
	{
		fprintf(fpProtocol, "The number of ALNs to be averaged in bagging is %d\n", nALNs);
		fprintf(fpProtocol, "The dimension of the problem (inputs + one desired output) is %d\n", nDim);
		fprintf(fpProtocol, "The output variable is %s\n", varname[nInputCol[nOutputIndex]]);
	}

	// allocate space for the minima, maxima and standard deviations of the variables
	adblMinVar = (double *)malloc(nALNinputs * sizeof(double));
	adblMaxVar = (double *)malloc(nALNinputs * sizeof(double));
	adblStdevVar = (double *)malloc(nALNinputs * sizeof(double));
} // end preprocUniversalFile

void ALNAPI createTVTSfiles()  // routine
{
  // Taking the PreprocessedDataFile, this creates
  // a test file TSfile of nRowsTS rows
  // and a file TVfile, a Training/noise_Variance file
  // formed by the remaining rows.
  // The number of samples in TSfile
  // nRowsTS is determined 
  // elsewhere either by a fraction of the
  // data set size (default 10%) or by the size of the
  // size of the test file returned by the Options dialog
	long i = 0, k = 0;
	// Training: choose from the PreprocessedDataFile randomly and without replacement nRowsTS rows
	// and put them into a file which is written to disk
	// Evaluation: of the PreprocessedDataFile is the TSfile.

  // First we set up the Test file in all cases
  // Now we produce the TV file (only the the case of training)
 	// get the remaining rows of the PreprocessedDataFile and produce the TVfile
	// which is written to disk
	if(bTrain) // TVfile is created only when training
	{
		if (nRowsTS == 0) // nothing goes into the TSfile
		{
			// in this case, we can take all of the preprocessed data file for the TV file
			TVfile.Create(nRowsPP, nALNinputs);
			nRowsTV = nRowsPP;
			double dblVal;
			for (i = 0; i < nRowsPP; i++)
			{
				for (int j = 0; j < nALNinputs; j++)
				{

					dblVal = PreprocessedDataFile.GetAt(i, j, 0);
					TVfile.SetAt(i, j, dblVal, 0);
				}
			}
		}
		else // nRowsTS > 0
		{
			int* anInTest = (int*)malloc(nRowsPP * sizeof(int));
			int nRandomRow;
			for(int j = 0; j < nRowsPP; j++)
			{
				anInTest[j] = 1; // everything starts off in the TV file
			}
			for (k = 0; k < nRowsTS; k++)
			{
				do
				{
					nRandomRow = (long)floor(ALNRandFloat() * (double)nRowsPP);
				} while (anInTest[nRandomRow] == 0); // keep doing this until you get a row that is still 1
				ASSERT(anInTest[nRandomRow] == 1);
				anInTest[nRandomRow] = 0; // then set that to 0
			}
			// anInTest now has nRowsTS 0 values and the rest 1's
			// TVfile is all of the PPfile except what went to the TSfile
			nRowsTV = nRowsPP - nRowsTS;
			TVfile.Create(nRowsTV, nALNinputs);
			TSfile.Create(nRowsTS, nALNinputs);
			double dblVal;
			long TVrows = 0;
			long TSrows = 0;
			for (i = 0; i < nRowsPP; i++)
			{
				if (anInTest[i] == 1)
				{ // case anInTest[i] == 1 and the row goes into TVfile
					for (int j = 0; j < nALNinputs; j++)
					{
						dblVal = PreprocessedDataFile.GetAt(i, j, 0);
						TVfile.SetAt(TVrows, j, dblVal, 0);
					}
					TVrows++;
				}
				else
				{ // case anInTest[i] == 0 and the row goes into TSfile
					for (int j = 0; j < nALNinputs; j++)
					{
						dblVal = PreprocessedDataFile.GetAt(i, j, 0);
						TSfile.SetAt(TSrows, j, dblVal, 0);
					}
					TSrows++;
				}
			} // end of writing TVfile and TSfile
			ASSERT((TVrows == nRowsTV) && (TSrows == nRowsTS));
		}
		if (bPrint && bDiagnostics)
		{
			//TVfile.Write("DiagnoseTVfile.txt");
			//fprintf(fpFileSetupProtocol, "DiagnoseTVfile.txt written with %d rows.\n", nRowsTV);
			//fflush(fpFileSetupProtocol);
			//TSfile.Write("DiagnoseTSfile.txt");
			//fprintf(fpFileSetupProtocol, "DiagnoseTSfile.txt written with %d rows.\n", nRowsTS);
			//fflush(fpFileSetupProtocol);
		}
	} // end if(bTrain)
	else
	{
		// This is an evaluation.
		nRowsTS = nRowsPP;
		TSfile.Create(nRowsTS, nALNinputs);
		double dblVal;
		for (int i = 0; i < nRowsPP; i++)
		{
			for (int j = 0; j < nALNinputs; j++)
			{
				dblVal = PreprocessedDataFile.GetAt(i, j, 0);
				TSfile.SetAt(i, j, dblVal, 0);
			}
		}
	} // end of writingTSfile
	if (bPrint && bDiagnostics)
	{
		//TSfile.Write("DiagnoseTSfile.txt");
		//fprintf(fpFileSetupProtocol, "DiagnoseTSfile.txt written with %d rows.\n", nRowsTS);
		//fflush(fpFileSetupProtocol);
	}
	fflush(fpProtocol);
} // end of createTVTSfiles

void ALNAPI analyzeTV() // routine
{
  fprintf(fpProtocol,"\n************** Analysis of TV file begins ********\n");
  fflush(fpProtocol);
	// This is done for training to set up the properties of the TVset
	ASSERT(TVfile.RowCount() >= nRowsTV);
  fprintf(fpProtocol,"Dimension nDim is %d, nRowsTV is %d\n", nDim, nRowsTV);
  fflush(fpProtocol);
	// Suppose there are nRowsTV input data points in an nDim - 1 dimensional cube.
	// Consider a uniform distribution on [-sqrt(3) stdev, +sqrt(3) stdev] in each axis
	// where sqrt(3) ~ 1.732. We start with using stdev units on each input axis.
	// The volume is ,3.464^(nDim-1)and each point receives its share of the unit volume.
	// The side of an nDim - 1 dimensional box with that per-point volume? 
	double dblSide = 3.464 * pow(1.0/((double) nRowsTV), 1.0/((double) nDim - 1.0 ));  // volume per TV point, side of above cube is the nDim -1 root.
	fprintf(fpProtocol,"Side of cube per point in the TVset (in stdev units) = %lf\n" ,dblSide );
	adblEpsilon = (double *) malloc((nDim) * sizeof(double));
  // we compute the variance of each column of TV,
	// and the maximum, minimum and standard deviation of each variable
	for(int k = 0; k < nDim; k++)
	{
		// initialize the min and max variables
		adblMinVar[k] = adblMaxVar[k] = TVfile.GetAt(0,k,0);
	}
  bRegress = FALSE; // if it can't be a classification problem, it must be regression
  for(int k = 0; k < nDim; k++) // do each variable k
  {
    //compute the max and min and average of variable k in TVfile
    double desired = 0,se  = 0, value = 0;

    for(int j = 0; j < nRowsTV; j++) 
    {
			value =  TVfile.GetAt(j,k,0);
      desired += value;
			if(value > adblMaxVar[k])
			{
				adblMaxVar[k] = value;
			}
			if(value < adblMinVar[k])
			{
				adblMinVar[k] = value;
			}
      // we want to see whether the output variable forces a regression since it is not an integer or an integer but out of range
      if((k == nDim - 1) && (fabs(floor(value + 0.5) - value)  > 1e-10))
      {
        bRegress = TRUE;
      }
      if((k == nDim - 1) && (value >3.5 || value < -3.)) // we can only classify into classes numbered -3, -2, ...2, 3
      {
        bRegress = TRUE;
      }
    }
    desired /= nRowsTV; // now desired holds the average for variable k
    // compute the standard deviation of variable k in TVset
    se = 0;
    double temp;
    for(int j=0; j < nRowsTV; j++) 
    {
      temp = TVfile.GetAt(j,k,0);
      se += ( temp - desired) * (temp - desired);
    }
    se /= ((double) nRowsTV - 1.0); // sample variance of variable k
		adblStdevVar[k] = sqrt(se);
		if(k < nDim - 1)
		{
			// If the per-point square box above is changed from stdev units to the input units,
			// we get a rectangular box of side adblEpsilon[k] that could contain each point
			// if there were a rectangular grid.
      adblEpsilon[k] = dblSide * adblStdevVar[k];  // changed WWA 2009.10.05
      if(nLag[k] == 0)
      {
        fprintf(fpProtocol, "Stdev of %s = %lf . Epsilon = %lf \n",varname[nInputCol[k]], adblStdevVar[k],adblEpsilon[k]);
      }
      else
      {
        fprintf(fpProtocol, "Stdev of %s@lag%d: = %lf . Epsilon = %lf \n",varname[nInputCol[k]], nLag[k], adblStdevVar[k],adblEpsilon[k]);
      }
    }
		else
		{
      fprintf(fpProtocol, "Stdev of output variable %s = %lf \n",varname[nInputCol[nOutputIndex]],adblStdevVar[k]);
			fprintf(fpProtocol, "The Epsilons above are sides of boxes per point in units of the particular input\n");
		}
  }
  fflush(fpProtocol);
}

/*
void getTVfile() // routine -- not currently used, but this and getTSfile() may be useful for diagnosis
{
	fprintf(fpProtocol,"Opening TAB separated Training/noise_Variance data file.\n");
  if (!TVfile.Read("TVfile.txt"))
  {
		fprintf(fpProtocol,"Reading TVfile failed!\n");
  }
  else
  {
		fprintf(fpProtocol,"Reading TV file succeeded!\n");
  }
  ASSERT(TVfile.ColumnCount() == nDim);
  fprintf( fpProtocol,"TVfile has %d full rows.\n", nRowsTV);
	fflush(fpProtocol);

}

void getTSfile()
{
	fprintf(fpProtocol,"Opening TAB separated test data file.\n");
  if (!TSfile.Read("TSfile.txt"))
  {
		fprintf(fpProtocol,"Reading TSfile failed!\n");
  }
  else
  {
		fprintf(fpProtocol,"Reading TSfile succeeded!\n");
  }
  ASSERT(!bTrain ||(TSfile.ColumnCount() == nDim) && ((TSfile.ColumnCount() == nDim)||(TSfile.ColumnCount() == nDim-1)));
	fprintf(fpProtocol,"TSfile has %d rows.\n",TSfile.RowCount());
  ASSERT(TSfile.RowCount() == nRowsTS);
	fflush(fpProtocol);
}
*/

void ALNAPI evaluate() // routine
{
	// loads the named DTREE file, determines nDim from it and accordingly
	// sets up the Output data file
	long lVersion;          /* DTREE library version */
	DTREE* pDtree;          /* pointer to DTREE structure */
	int nErrCode;           /* library function return value */
	char szErrMsg[256];     /* error message */
	double dblMax, dblMin;  /* temporaries to hold bounds */
	double dblOutput;       /* result of Dtree evaluation */
	int k;                  /* loop counters */
	fprintf(fpProtocol, "\n*********** Opening the DTREE file ****************\n");
	/* get dtree version */
	lVersion = GetDtreeVersion();
	printf("DTREE library v%d.%d\n", lVersion >> 16, lVersion & 0x0000FFFF);

	// load the DTREE file from the earlier run of ALNfit
	fprintf(fpProtocol, "Opening DTREE file %s\n", szDTREEFileName);
	nErrCode = ReadDtree(szDTREEFileName, &pDtree);

	/* check error return */
	if (nErrCode == DTR_NOERROR)
	{
		fprintf(fpProtocol, "DTREE succesfully parsed!\n");
	}
	else
	{
		GetDtreeError(nErrCode, szErrMsg, sizeof(szErrMsg));
		fprintf(fpProtocol, "\nError (%d): %s\n", dtree_lineno, szErrMsg);
		fflush(fpProtocol);
		return;
	}
	/* succesfully loaded DTREE */

	// now look at the data file

	int ncols = nALNinputs; // during evaluation this may be different from nDim
	int nrows = nRowsTS;
	/* compare the output variable index and nDim -- use ncols from analyzeinputfile*/
	if ((nALNinputs != pDtree->nOutputIndex + 1) && (nALNinputs != pDtree->nOutputIndex))
	{
		// the evaluation data file can't have a correct number of columns
		fprintf(fpProtocol, "The number of columns in the data file, with or without outputs,\n");
		fprintf(fpProtocol, "doesn't match the DTREE number of columns.\n");
		fflush(fpProtocol);
		sprintf(szResultMessage, "Data file has the wrong number of colums to be used with this DTREE");
		return; // this case still has to be worked on
	}
	else
	{
		//we have the output column in the data file or else we shall replace it with zeros below
		nDim = pDtree->nOutputIndex + 1;
		nOutputIndex = nDim - 1;
	}
	fprintf(fpProtocol, "Dimension of the DTREE is %d \n", nDim);
	// initialize the random number generator using the time
	srand((unsigned)time(NULL));
	// generate one random vector to evaluate
	// the output component is also random
	// and will be replaced by the DTREE evaluation
	for (k = 0; k < nDim - 1; k++)
	{
		/* get range of all input variables */
		dblMax = pDtree->aVarDefs[k].bound.dblMax;
		dblMin = pDtree->aVarDefs[k].bound.dblMin;
		if (nLag[k] == 0)
		{
			fprintf(fpProtocol, "Variable %s dblMin = %f dblMax = %f \n", varname[nInputCol[k]], dblMin, dblMax);
		}
		else
		{
			fprintf(fpProtocol, "Variable %s@lag%d: dblMin = %f dblMax = %f \n", varname[nInputCol[k]], nLag[k], dblMin, dblMax);
		}
	}
	if (OutputData.Create(nrows, nDim + 1) != NULL)
	{
		fprintf(fpProtocol, "\n******** Testing the DTREE on the test set **********\n");
		fprintf(fpProtocol, "Creating internal data file OutputData succeeded!\n");
		fprintf(fpProtocol, "It has %d rows and %d columns.\n", nrows, nDim + 1);
	}
	else
	{
		fprintf(fpProtocol, "Creating internal data file OutputData failed!\n");
		fflush(fpProtocol);
		exit(0);
	}


	// copy the test file into the output data file
	ASSERT((ncols == nDim) || (ncols == nDim - 1));
	double value;
	for (int j = 0; j < nrows; j++)
	{
		for (k = 0; k < ncols; k++)
		{
			value = TSfile.GetAt(j, k, 0);
			OutputData.SetAt(j, k, value, 0);
		}
	}
	int nCountMisclassifications = 0;
	double dblSE = 0, dblMAE = 0, dblMAXE = 0; // three error accumulators  N.B. this local dblSE has nothing to do with the SmoothingEpsilon abbreviation used elsewhere. 
	double * adblX = (double *)malloc((nDim) * sizeof(double));
	for (int j = 0; j < nRowsTS; j++)
	{
		/* get DTREE evaluation (dblOutput) */
		for (k = 0; k < nDim - 1; k++)
		{
			/* get inputs*/
			adblX[k] = TSfile.GetAt(j, k, 0);
		}
		adblX[nDim - 1] = 0; // File value not used 
		if ((nErrCode = EvalDtree(pDtree, adblX, &dblOutput, NULL)) != DTR_NOERROR)
		{
			GetDtreeError(nErrCode, szErrMsg, sizeof(szErrMsg));
			fprintf(fpProtocol, "\nError (%d): %s\n", dtree_lineno, szErrMsg);
		}
		// put the result into the output data file(round to integer if classification)
		if (bClassify)
		{
			OutputData.SetAt(j, nDim, floor(dblOutput + 0.5), 0);
		}
		else
		{
			OutputData.SetAt(j, nDim, dblOutput, 0);
		}
		if (ncols < nDim) // we put zeros into the output data file in the column of missing outputs
		{
			OutputData.SetAt(j, nDim - 1, 0, 0);
		}
		else
		{
			// we compute the square error the absolute error and the maximum of the latter
			double dblAbsErr;
			dblAbsErr = fabs(dblOutput - TSfile.GetAt(j, nDim - 1, 0));
			dblMAE += dblAbsErr;
			dblSE += dblAbsErr * dblAbsErr;
			if (dblAbsErr > dblMAXE) dblMAXE = dblAbsErr;
		}

		// if we have the correct outputs (even for classifications, we allow slight errors)
		if (bClassify && ncols == nDim)
		{
			if (floor(dblOutput + 0.5) != floor(OutputData.GetAt(j, nDim - 1, 0) + 0.5)) nCountMisclassifications++;
		}
	}
	fprintf(fpProtocol, "The following is based on a test set with %d samples\n", nRowsTS);
	fprintf(fpProtocol, "Please be careful interpreting the results for small numbers of samples!\n");
	if (bClassify && ncols == nDim)
	{
		fprintf(fpProtocol, "The number of misclassifications on the test set is %d\n",
			nCountMisclassifications);
		nEvalMisclassifications = nCountMisclassifications;
		dblEvalMisclassificationPercent = 100.0*(double)nCountMisclassifications / (double)nRowsTS;
		sprintf(szResultMessage, "Samples misclassified = %5.2f percent of test set", 100.0*(double)nCountMisclassifications / (double)nRowsTS);
		fprintf(fpProtocol, szResultMessage);
		fprintf(fpProtocol, "\n");
	}
	else
	{
		dblEvalRMSError = sqrt(dblSE / (double)nRowsTS); // put into global value
		sprintf(szResultMessage, "RMS deviation of DTREE output from desired\n(or from 0 if output column not present) is %f ",
			dblEvalRMSError);
		fprintf(fpProtocol, szResultMessage);
		fprintf(fpProtocol, "\n");
		fprintf(fpProtocol, "Mean absolute deviation of DTREE output from desired is %f\n",
			dblMAE / (double)nRowsTS);
		fprintf(fpProtocol, "Maximum absolute deviation of DTREE output from desired is %f\n", dblMAXE);
	}
	// write out the evaluation file
	if ((fpOutput = fopen(szScatterFileName, "w")) != NULL)
	{
		fprintf(fpProtocol, "\nPlease examine the output file named:\n   %s\n", szScatterFileName);
		fprintf(fpProtocol, "The rightmost column of the output file is the ALN prediction.\n");
		fflush(fpProtocol);
	}
	else
	{
		fprintf(fpProtocol, "Opening file %s failed. Stopping.\n", szScatterFileName);
		exit(0);
	}
	// print out the header lines as they were read
	if (nHeaderLinesUniv > 0)
	{
		for (int linenum = 0; linenum < nHeaderLinesUniv; linenum++)
		{
			fprintf(fpOutput, "%s\n", hdrline[linenum]);
		}
	}
	else // there are no header lines, just output the substitute variable names
	{
		fprintf(fpOutput, "Generated column headings\n");
		for (int k = 0; k < nDim; k++)
		{
			fprintf(fpOutput, "%s", varname[k]);
		}
		fprintf(fpOutput, "\n");
	}
	fflush(fpOutput);
	// now remember that the header shows the columns, but not the ALN inputs
	// which are columns combined with lags!
	fprintf(fpOutput, "The data file columns (and lags) used for ALN inputs and output are:\n");
	for (k = 0; k < nDim; k++)
	{
		if (nLag[k] == 0)
		{
			fprintf(fpOutput, "%s\t", varname[nInputCol[k]]);
		}
		else
		{
			fprintf(fpOutput, "%s@lag%d\t", varname[nInputCol[k]], nLag[k]);
		}
	}
	fprintf(fpOutput, "ALNfitDeep\n");
	fflush(fpOutput);
	// print out the data lines from OutputData
	double dblValue;
	char szValue[32];
	int charcount;
	for (int linenum = 0; linenum < nRowsTS; linenum++)
	{
		for (k = 0; k < nDim + 1; k++)
		{

			dblValue = OutputData.GetAt(linenum, k, 0);
			sprintf(szValue, "%f", dblValue);
			charcount = 0;
			if (!bDecimal)
			{
				while (szValue[charcount] != 0)
				{
					if (szValue[charcount] == '.') szValue[charcount] = ',';
					charcount++;
				}
			}
			fprintf(fpOutput, "%s", szValue);
			if (k < nDim)
			{
				fprintf(fpOutput, "\t");
			}
			else
			{
				fprintf(fpOutput, "\n");
			}
		}
	}
	fflush(fpOutput);
	// Now we do the replacement of missing values in the data file
	// if bReplaceUndefined is true.
	// filename is the sames as the .dtr file with an R appended before .txt.
	// We use the UNfile for this.  It has the whole data of the original data file
	// except for the header lines, which we have stored in the aline array.
	// The rule used for replacement appears at the top of the output. The
	// rules used will be shown with the first used at the top.
	// Blank lines and comments are suppressed in the data part of the file.

	//if(bReplaceUndefined) we always do this, even if creating an extra column
	//{
	int nMaxLag = 0;
	for (int ninput = 0; ninput < nALNinputs - 1; ninput++) // the output has lag 0 which doesn't affect the max
	{
		if (nLag[ninput] > nMaxLag) nMaxLag = nLag[ninput];
	}
	for (int i = nMaxLag; i < nRowsUniv; i++)
	{
		BOOL bRowGood = TRUE; // this means the DTREE can be used to compute an output value and the value in the file is undefined
		// see if the output value in the file is already defined
		ASSERT(nLag[nALNinputs - 1] == 0); // a non-zero lag on the output variable is not allowed
		adblX[nALNinputs - 1] = UNfile.GetAt(i, nInputCol[nALNinputs - 1]);
		// we don't want to generate the value using the DTREE if 
		// we are replacing values in a copy of the input data file, and the value to be replaced is defined
		// or if we can't compute a value because ther is an undefined input value
		if (bReplaceUndefined && (fabs(adblX[nALNinputs - 1] - 99999) > 0.5) && (fabs(adblX[nALNinputs - 1] + 9999) > 0.5))
		{
			bRowGood = FALSE; // the output is defined, so replacement is not needed
		}
		else
		{
			// see if any input to the DTREE is undefined 
			for (int ninput = 0; ninput < nALNinputs - 1; ninput++)
			{
				adblX[ninput] = UNfile.GetAt(i - nLag[ninput], nInputCol[ninput]);
				if ((fabs(adblX[ninput] - 99999) < 0.5) || (fabs(adblX[ninput] + 9999) < 0.5)) bRowGood = FALSE;
			}
		}
		if (bRowGood == TRUE)
		{
			// in this case we have the information necessary to compute a substitute for
			// the missing output value
			if ((nErrCode = EvalDtree(pDtree, adblX, &dblOutput, NULL)) != DTR_NOERROR)
			{
				GetDtreeError(nErrCode, szErrMsg, sizeof(szErrMsg));
				fprintf(fpProtocol, "\nError (%d): %s\n", dtree_lineno, szErrMsg);
			}
			else
			{
				// put the generated value in the right place
				if (bReplaceUndefined)
				{
					UNfile.SetAt(i - nLag[nALNinputs - 1], nInputCol[nALNinputs - 1], dblOutput, 0);
				}
				else
				{
					UNfile.SetAt(i - nLag[nALNinputs - 1], nColsUniv, dblOutput, 0);
				}
			}
		}
	}
	//UNfile.Write(szR_or_E_FileName); // changed to include header
	if ((fpReplacement = fopen(szR_or_E_FileName, "w")) != NULL)
	{
		fprintf(fpProtocol, "\nPlease examine the R or E file named:\n   %s\n", szR_or_E_FileName);
		fprintf(fpProtocol, "Where possible, missing data file output values have been computed.\n");
		fprintf(fpProtocol, "The rule at the top shows the input variables and output column used.\n");
		fprintf(fpProtocol, "Using this file as an input data file, you can do further replacements.\n");
		fflush(fpProtocol);
	}
	else
	{
		fprintf(fpProtocol, "Opening file %s failed. Stopping.\n", szR_or_E_FileName);
		exit(0);
	}
	// copy all the existing rules at the top of the file in the old header
	int nNewRuleNumber = 1;
	for (int linenum = 0; linenum < nHeaderLinesUniv; linenum++)
	{
		if (hdrline[linenum][0] == 'R' &&
			hdrline[linenum][1] == 'u' &&
			hdrline[linenum][2] == 'l' &&
			hdrline[linenum][3] == 'e')
		{
			fprintf(fpReplacement, "%s\n", hdrline[linenum]);
			nNewRuleNumber = linenum + 2;
		}
	}
	if (bReplaceUndefined)
	{
		// now we insert the new replacement rule in the next line:
		fprintf(fpReplacement, "Rule %d: ", nNewRuleNumber);
		for (k = 0; k < nALNinputs - 1; k++)
		{
			if (nLag[k] == 0)
			{
				fprintf(fpReplacement, "%s\t", varname[nInputCol[k]]);
			}
			else
			{
				fprintf(fpReplacement, "%s@lag%d\t", varname[nInputCol[k]], nLag[k]);
			}
		}
		fprintf(fpReplacement, " -> ");
		fprintf(fpReplacement, "%s\n", varname[nInputCol[nALNinputs - 1]]);
	}
	else
	{
		fprintf(fpReplacement, "The rightmost column contains computed values for %s\n", varname[nInputCol[nALNinputs - 1]]);
	}
	// then we put in the non-rule lines of the old header
	for (int linenum = nNewRuleNumber - 1; linenum < nHeaderLinesUniv; linenum++)
	{
		fprintf(fpReplacement, "%s\n", hdrline[linenum]);
	}
	// finally we output the data lines
	// if this is not replacement, then we have to output an extra column from the UNfile
	int nColsPrinted = nColsUniv;
	if (!bReplaceUndefined) nColsPrinted++;
	for (int i = 0; i < nRowsUniv; i++)
	{
		// this has to output commas if required
		for (int j = 0; j < nColsPrinted; j++)
		{
			dblValue = UNfile.GetAt(i, j);
			sprintf(szValue, "%f", dblValue);
			charcount = 0;
			if (!bDecimal)
			{
				while (szValue[charcount] != 0)
				{
					if (szValue[charcount] == '.') szValue[charcount] = ',';
					charcount++;
				}
			}
			fprintf(fpReplacement, "%s", szValue);
			if (j < nColsPrinted - 1)
			{
				fprintf(fpReplacement, "\t");
			}
			else
			{
				fprintf(fpReplacement, "\n");
			}
		}
	}
	fclose(fpReplacement);
	// cleanup
	DestroyDtree(pDtree);
	free(adblX);
	UNfile.Destroy();
	OutputData.Destroy();
} // end evaluate
