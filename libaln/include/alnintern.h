// ALN Library (libaln)
// Copyright (C) 1995 - 2010 William W. Armstrong.
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
void splitcontrol(CMyAln*, double); // if average variance error of a piece is high, splitting is prevented
void dosplitcontrol(CMyAln*, ALNNODE*, double); //does the recursion of splitcontrol
void dozerospliterror(CMyAln*, ALNNODE*);  // sets the square error to zero in each LFN
void spliterrorsetTR(CMyAln*); // accumulates the training square error and number of hits on each linear piece
void spliterrorsetVAR(CMyAln*); // accumulates the variance square error and number of hits on each linear piece
void dodivideTR(CMyAln*, ALNNODE*); // divides the total square training set errors of the pieces by their hit count
void dodivideVAR(CMyAln*, ALNNODE*); // divides the total square variance errors of the pieces by their hit count

// thread procedures
UINT TakeActionProc(LPVOID pParam);  // separate thread


//global variables used only internally

extern CMyAln** apALN;       // an array of pointers to ALNs
//extern CMyAln** apValALN;    // an array of pointers to ALNs to be used on the variance set
//extern static CMyAln* pAvgALN;      // an ALN representing the bagged average
extern BOOL bTrainingAverage;// Switch to tell fillvector whether get a training vector or compute an average
extern int nNumberEpochs;
extern int nNumberLFNs;  // used to control the epoch size, which should be proportional to this
//extern static char szVarName[100][3];
//int nColsUniv = 0;
//long nRowsUniv = 0;
extern int nColsAuxVariance;
extern int nColsAuxTest;
extern long nRowsTR; // size of training file
extern long nRowsVAR; // size of variance file
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
extern int* anInclude;


// Files used only internally
extern FILE *fpData;                 // The data file which contains all data.
extern FILE *fpProtocol;             // the file to record the results of the experiment. Rename to save.
//FILE *fpOutput = NULL;             // the output data file resulting from evaluation (now in ALNfit.cpp)
extern FILE *fpReplacement;
extern CDataFile UNfile;             // copy of the data file, but with missing values replaced by special number
extern CDataFile PreprocessedDataFile;             // The preprocessed file created from the Universal file
extern CDataFile TVfile;             // The file used for training and, if no separate file is given, for variance, with nDim columns and nRowsUniv - nRowsTSfile rows.
extern CDataFile TRfile;             // Training file.  This file is setup separately for each ALN to implement bagging

extern double* adblLRW; // stores an ALN weight approximation from linear regression
extern double* adblLRC; // ditto for centroids

// this typedef allows a cast within ALNfit Pro of ALNLFNSPLIT in aln.h which uses
// the variables in a different way than in the Dendronic Learning Engine SDK
// when a comparison between the average errors on training and variance sets
// is needed to decide whether or not to allow splitting of the hyperplane
/*typedef struct tagSPLIT      // Used in inhibiting splitting -- must be zeroed before and after use
{
  int nCount;                // Number of hits
  double dblSqErrorTrain;    // Squared error of a piece during training                      
  double dblSqErrorVal;      // Squared error of a piece during variance
  double dblT_NotUsed;			 // Used in the SDK only
} SPLIT;
*/