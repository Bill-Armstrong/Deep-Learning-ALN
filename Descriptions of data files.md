

Flat.txt
has 4 inputs one output, labelled W X Y Z OUTPUT.
A is 3*RAND(), B is -5 * RAND(), C is 7*RAND(), D is 2*RAND()
The output is 6*A1-15*B1+4*C1-9*D1 with no noise added.
There are 4000 samples.
Experiments
Average Noise Variance before smoothing 3.718338 
 Average Noise Variance after smoothing 2.958296 
These two should have been zero! The deviation from zero This is obviously
because of the slope of the function. Slope is a significant
influence on the noise variance!!!

NoisyConstant.txt
has 3 inputs one output. The inputs are all 10 * RAND() and the output is 10*(RAND()-0.5).
The noise variance is constant 8.3333
RMSE of the best linear fit 2.8867
There are 10000 samples.
Average Noise Variance before smoothing 8.388240 
 Average Noise Variance after smoothing 8.308575 
These values are close to 8.33, showing that 0 slope makes things easier for
estimating noise variance.  Could we take a larger number of points, e.g. 2*nDim,
and estimate the local slope??? Could we use generalized cross-validation on 2 * nDim
points, then use only the one?

NoisyFlat.txt
has 8 inputs and one output labelled A B... H OUTPUT
The inputs are all rand() and the output is linear with
coefficients 1 ,2 ,-3, 4,-5, 7, 10, -1
The output has noise 10 * (RAND()-0.5)
The noise variance is constant 8.3333
RMSE of the best linear fit 2.8867
There are 5000 samples.
Experiment on NoisyFlat.txt
Noise variance after smoothing 8.83 on 4500 samples not used in testing.
Scatterplot has noise on the diagonal + or -5 as expected.
RMSE of linear regression 2.894 after iteration 0 -- very fast.
Weights after 29 epochs .82, 2.1, -2.9, 3.73, -5.09, 7.29, 9.82, -1.2
are close to the ones used in setting up the problem.
Smoothing reduced the noise variance from 9.23 to 8.83. The resulting samples
were concentrated mostly between 8.77 and 8.9.  This suggests that
the slope of the function has added about 0.5 to the estimated noise varinace.
There are no big outliers after smoothing.
The ALNs with the F-test did not change from the initialization, error 0.

NoisySinCos.txt
has two inputs one output
X and Y are both 62.8*RAND() while the output Z is
 =100* SIN(0.15*A1)*COS(0.15*B1) + 10.0 *(RAND() - 0.5)
There are 6000 samples.
Experiments
The linear regression RMSE went down to 51.38 in epoch 0 and then
varied around 50.
Average Noise Variance before smoothing 14.782750 
 Average Noise Variance after smoothing 12.319474 
This should have been around 8.3, so there is an influence of slope
and curvature.
The result is a good scatterplot with some noise.

SinNoise1000.txt
has one input and one output.
The input varies from 1 to 1000, and is given by
=(1000-A1)*SIN($E$1*A1+$E$2*A1*A1) with $E$1 = 0.01 $E$2 = 0.0001
The added noise is given by $E$3*2*A1*(RAND()-0.5)/1000 with $E$3 = 200
so the noise amplitude goes from 0.2 at A1=1 to 200 at A1 = 1000.
The formula for the amplitude is 0.2 A1.
The noise variance is 1/3 of the square of the amplitude.
 It's formula is 0.04 A1^2. This exact formula for noise variance
was used with no test set and still the fit was only good for the top seven
peaks. Could this be improved? Probably not.The training set has only 1000 samples,
and sometimes the noise deviates from the theoretical noise. If we had more training
samples, then the curve would get better, but the noise would only get closer to ideal.

NoisySin.txt
has one input and one output but
it has 10 times the density of sampling of SinNoise1000.txt
