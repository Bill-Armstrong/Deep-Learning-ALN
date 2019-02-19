
How to generate a test file using a spreadsheet

## NoisySin20000.txt

* 20000 samples

* 1 input: a sequence from 0 to 999.95; start with A1=0, A2=A1+0.05

* output: (copy into B1 and C1 the expressions below)

=(1000-A1)\*SIN(0.01\*A1+0.0001\*A1\*A1)

=B1+8\*(1.00461579^A1)\*2\*(RAND()-0.5)

Now copy B1 and C1 into B2 and C2, then copy all columns to get 20000 samples.

When you Browse to this file and open it in ALNfitDeep, use the Previous connection and Remove buttons to exclude the correct unknown function in column B from training. It will show up in the "E" output file for evaluating accuracy of the result.

The signal shrinks from 1000 to 0 and the noise grows exponentially from amplitude 8 to 800.
