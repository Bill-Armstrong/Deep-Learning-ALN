# Coming Soon

ALNfitDeep is open-source machine learning software for Windows. This new version's main improvement is that it can handle a varying level of noise at different places in the domain of the unknown function. Please see the .png image in this directory for an example of how ALN learning with noise variance stopping can bring out details in the signal even under noisy conditions.  You can experiment with a pre-release version v2.0alpha of the software. There is a Help button on the interface with instructions. You can generate a test file using the  following description. Using the ALNfitDeep software you can change the line dblLimit = -1.0 to dblLimit = 10.0 in the approximation routine. The former does noise variance stopping and the latter overtrains.

## VeryNoisyGaussianSin2000.txt

* 2000 samples

* 1 input: a sequence from 0 to 999.5; start with A1=0, A2=A1+0.5

* output: (copy into B1 and C1 the expressions below)

=(1000-A1)\*SIN(0.01\*A1+0.0001\*A1\*A1)

=B1+10\*(1.00461579^A1)\*2\*(NORMINV(RAND(),0,1))

Now copy B1 and C1 into B2 and C2, then copy all columns to get 2000 samples.

When you Browse to this file and open it in ALNfitDeep, use the Previous connection and Remove buttons to exclude the correct unknown function in column B from training. It will show up in the "E" output file for evaluating accuracy of the result.

The signal shrinks from 1000 to 0 and the noise grows exponentially from amplitude 10 to 1000.


