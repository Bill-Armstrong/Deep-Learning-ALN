##About ALNfitDeep:

ALNfitDeep is open-source machine learning software for MS Windows that **everyone** can use because it is COMPLETELY AUTOMATED.  You just say which tab-separated file of numbers contains your data. If the right-hand column is the one to be learned, you just click on Start and the ALNfitDeep program does the rest. All parameters are set AUTOMATICALLY. The architecture grows AUTOMATICALLY to fit the problem until you have a solution. Overtraining is avoided by estimating the noise level in the data and not growing the architecture when it tends to learn noise. The DEEP LEARNING algorithm is extremely fast in both learning and execution.

###The program can do:
 
* non-linear, non-parametric regression
* filling in of missing values in data files
* multiple time-series analysis
* high-dimensional problems
* imposition of monotonicity constraints on the function to be learned
* imposition of bounds on partial derivatives of the function to be learned
* identification of input variables that can be removed to improve results
* classification into two classes (and more under some conditions)

*Most importantly, the ideas embodied in the software are key for future development of **DEEP LEARNING**.*

The download from GitHub is a bit more than 10 MB:
 **git clone https://github.com/Bill-Armstrong/Deep-Learning-ALN.git**

The program needs two open source libraries: Boost and Eigen.  It expects to find C:\Boost\boost_1_60_0 and C:\Eigen\eigen-eigen-07105f7124f9 or later versions.  Boost provides some special functions and Eigen provides the singular value decomposition (SVD) both for use in stats. You can find them for (free) download at

http://www.boost.org/
http://eigen.tuxfamily.org/index.php?title=Main_Page

###What to do when you have downloaded

This open-source project builds under Visual Studio 2015 Community (free)

https://www.visualstudio.com/en-us/products/visual-studio-community-vs.aspx

It makes use of C++ and Microsoft Foundation Classes (MFC).  The solution (.sln) file is in the folder libaln/win32. Batch build the libaln and alnpp libraries first (10 items), then the ALNfitDeep and realestate items (6 items). The executable ALNfitDeep.exe appears in the libaln/samples/ALNfitDeep/Release directory. You should have a new folder to contain a copy of the program and your data files.

###To begin

Follow the first paragraph of "Examples for ALNfitDeep.md" for a run of ALNfitDeep.exe. Help gives brief suggestions. 

Best wishes

Bill Armstrong
