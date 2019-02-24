The iris and abalone data are classification problems, but you can also do them as regressions.  You just look at the scatterplot output to see how well the ALN system performed. Classification is not a strength of this approach -- yet.

The Devil's Tower is a function looking like an odd-shaped mountain.  Of course you can see how well the system did on the scatterplot, but you can also use GnuPlot after choosing the right columns in the "E" output file with commands like this:

set view 60,30,0.85,1

set title "ALN output from Devil's Tower"

set dgrid3d 61,61

set hidden3d

splot "1140DevilsTowerALNoutput.txt" u 1:2:3 with lines