# README

ALNfitDeep is open-source machine learning software for Windows. This new version's main improvement is that it can handle a varying level of noise at different places in the domain of the unknown function. Please see the .png image in this directory for an example of how ALN learning with noise variance stopping can bring out details in the signal even under noisy conditions.  You can now experiment with a release version of the software. There is a Help button on the interface with instructions. Just browse to your own file (tab-separated) with the desired output column at the right (or change the inputs and output as required).

# Experiments:
Using the ALNfitDeep software you can change the line dblLimit = -1.0 to dblLimit = 10.0 in the approximation routine. The former does noise variance stopping and the latter overtrains. You can generate a test file using the  following description. It is instructive to increase the file size by a factor of 10 to see what happens.

## VeryNoisyGaussianSin2000.txt

* 2000 samples (OR 20000 samples)

* 1 input: a sequence from 0 to 999.5; start with A1=0, A2=A1+0.5 (OR 0.05 if you want 20000 samples)

* output: (copy into B1 and C1 the expressions below)

=(1000-A1)\*SIN(0.01\*A1+0.0001\*A1\*A1)

=B1+10\*(1.00461579^A1)\*2\*(NORMINV(RAND(),0,1))

Now copy B1 and C1 into B2 and C2, then copy all columns to get the number of samples you want. Copy the three columns onto a new sheet, add any headers, e.g. X, Noiseless, Noisy, and save the sheet as a tab-separated text file.

When you Browse to this file and open it in ALNfitDeep, use the Previous connection and Remove buttons to exclude the actual values of the unknown function in column B from training. It will show up in the "E" output file for evaluating accuracy of the result.

The signal shrinks from 1000 to 0 and the noise grows exponentially from amplitude 10 to 1000.


