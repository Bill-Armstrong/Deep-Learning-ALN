// ALN Library sample
// A sample giving the same values as LINEST in MS Excel (TM).
// ALNfit Learning Engine for approximation of functions defined by samples.
// Copyright (C) 2018 William W. Armstrong 
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

// realestate.cpp
// This program does not involve learning. Its purpose is to compare statistics
// produced by the ALNfit Learning Engine to the values output by Microsoft (R)
// Excel in a realestate example in the Excel help file for linear regression
// (index term: LINEST). An ALN with a single linear piece is constructed
// having weights given in that example.  The output values are very close.

#ifdef __GNUC__
#include <typeinfo>
#endif
#include <iostream>
#include <stdio.h>
#include <alnpp.h>
#include <datafile.h>
#include <iomanip>

static char szInfo[] = "ALN Library RealEstate sample\n"
                       "Copyright (C)  2018 William W. Armstrong\n"
                       "Licensed under LGPL\n\n";
using namespace std;
class CMyAln : public CAln
{
// notification overrides
public:
  virtual BOOL OnTrainStart(TRAININFO* pTrainInfo, void* pvData)
    { 
      cout << "Training starts..." << endl; 
      return TRUE; 
    }

  virtual BOOL OnTrainEnd(TRAININFO* pTrainInfo, void* pvData) 
  { 
    cout << "Training finished.  RMSE: " << pTrainInfo->dblRMSErr << endl; 
    return TRUE; 
  }
  
  virtual BOOL OnEpochStart(EPOCHINFO* pEpochInfo, void* pvData) 
  { 
    cout << "Epoch " << pEpochInfo->nEpoch << ":" << flush; 
    return TRUE; 
  }
  
  virtual BOOL OnEpochEnd(EPOCHINFO* pEpochInfo, void* pvData) 
  {
    cout << "  Estimated RMSE: " << pEpochInfo->dblEstRMSErr << endl;
    return TRUE;
  }

};

int main()
{
  cout << szInfo;

  // open up data file
  cout << "Opening data file realestate.dat... ";
  
  CDataFile file;
  if (!file.Read("realestate.dat"))
  {
    cout << "failed!" << endl;
    return 1;
  }
  else
  {
    cout << "succeeded!" << endl;
  }

	

 	ASSERT(file.ColumnCount() == 5);
  ASSERT(file.RowCount() == 11);

  // seed ALN random number generator
  CAln::SRand(57);

  // create ALN 
  cout << "Creating 5 dimensional ALN, output index 4... ";
  CAln aln;
  if (!aln.Create(5, 4))
  {
    cout << "failed!" << endl;
    return 1;
  }
  else
  {
    cout << "succeeded!" << endl;
  }
	// set weights directly
	double adblW[] = { 52317.8305,
										 27.6413874,
										 12529.7682,
										 2553.21066,
										 -234.23716,
										 -1.0 };
	for(int ntest = 0; ntest < 2; ntest++)
	{
		if(ntest == 1)
		{
			// before the training test, zero the weights except for the output weight
			for(int j = 0; j < 5; j++)adblW[j] = 0.0;
		}
		ALNNODE* pLFN = aln.GetTree();
		ASSERT(NODE_ISLFN(pLFN));
		double* adblLFNW = LFN_W(pLFN);
		for (int i = 0; i < 6; i++)
		{
			adblLFNW[i] = adblW[i];
		}
	    
		// set data info
		int nPoints = file.RowCount();
		int nCols = file.ColumnCount();
		const double* adblData = file.GetDataPtr();
		aln.SetDataInfo(nPoints, nCols, adblData);
		if(ntest == 1)
		{
			BOOL bJitter = FALSE;
			double dblLearnRate = 0.1;  // small learning rate
			double dblMinRMSE = 0.00001;
			int nMaxEpochs = 10000;
			int nNotifyMask = AN_TRAIN|AN_EPOCH ;
			if (!aln.Train(nMaxEpochs, dblMinRMSE, dblLearnRate, bJitter, nNotifyMask))
			{
				cout << "Training failed!" << endl;
				return 1;
			}
		}
		// perform lfn analysis
		void* pvAnalysis = NULL;
		int nLFNStats = 0;
		aln.LFNAnalysis(pvAnalysis,nLFNStats,AN_NONE,NULL,NULL);

		// output results
		if(ntest == 0)
		{
			printf("\nTest of stats with pre-set ALN weights\n");
		}
		else
		{
			printf("\nTest of stats with trained ALN weights\n");
		}
		printf("LFN       R2        SEE       DF        F         \nFp        nW        RSS       ESS\n");
		//We don't print stats on the bias weight which always has value 1.0
		for (int i = 0; i < nLFNStats; i++)
		{
			LFNSTATS LFNStats;
			int nWeightStats = 0;
			ALNNODE* pLFN = NULL;
			aln.LFNStats(pvAnalysis, i, LFNStats, nWeightStats, pLFN);
			//removed pointer to LFN in printf below
			printf("%d%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%d\t%10.6f\t%10.6f\n",
							i,
							LFNStats.dblR2, 
							LFNStats.dblSEE, 
							LFNStats.dblDF, 
							LFNStats.dblF, 
							LFNStats.dblFp, 
							nWeightStats, 
							LFNStats.dblRSS, 
							LFNStats.dblESS);

			printf("\n\tx\tW\tSEw\tT\tTp\n");
			for (int j = 0; j < nWeightStats-1; j++)
			{
				LFNWEIGHTSTATS LFNWeightStats;
	      
				aln.LFNWeightStats(pvAnalysis, i, j, LFNWeightStats);
				printf("\t%d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n", j, 
							 LFNWeightStats.dblW, 
							 LFNWeightStats.dblSEw, 
							 LFNWeightStats.dblT, 
							 LFNWeightStats.dblTp);
			}
			printf("The output weight (index 4) is always -1, so line number 4 \n");
			printf("which is omitted contains the statistics of the bias weight instead.\n");
		}
		// free analysis data 
		aln.LFNFreeAnalysis(pvAnalysis);
	}

	char c;
	cin >> c;
  return 0;
}
