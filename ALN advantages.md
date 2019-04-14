25 ways in which ALNs (Adaptive Logic Networks) offer new possibilities for machine learning:

What is an ALN? An ALN is a type of neural net consisting of a first layer of affine functions whose outputs are combined in all other layers by two-input maximum and minimum functions to produce a single output.

1. ALN structure is very simple: Because all arithmetic is in the first layer, ALNs can be extremely deep without loss of speed.

2. ALNs are easy to understand and analyze: The output function is continuous. The graph is formed by pieces of the hyperplanes defined by the affine functions. Pieces are bounded by other pieces where they intersect.

3. ALNs can be smoothed: Inserting small quadratic or quartic "fillets" at max and min nodes makes the output once or twice continously differentiable.

4. Even smoothed ALN outputs are easy to analyze: The smoothed ALN can be made to differ from the piecewise linear surface by at most a prescribed amount.

5. ALNs elegantly solve the credit assignment problem: Given any input and the values computed in the ALN, it is easy to determine how the output is computed. At a max (or min) node, the input subtree with the greater (or lesser) input is responsible. Thus one affine function in the first layer is responsible for the output. It is the one whose weights are updated during training. In case of a (rare) tie, responsibility can be shared along two paths. This is also the case when smoothing is used. The ALN is not a "black box". There is no "vanishing gradient" problem.

6. ALN training is simple: When an affine function computing the output value is determined, an iterative form of least squares adjusts those weights. With tiny smoothing fillets, usually only one or a small number of affine functions have to be adjusted.

7. ALNs can capture convexity: If there are only max nodes in an ALN, the surface is convex down. Min nodes make it convex-up. The structure of an efficient ALN tree would follow the convexities of the ideal function. An inefficient ALN could fit any continuous function on a compact set with a maximum of minima (the tops of truncated pyramids). This is a simple proof of ALN's universal approximation capability.  

8. ALN evaluation is lazy: for example if you have a part of the ALN that looks like max(min(y1, y2), y3) with y1 = 5, y3 = 6, then the system doesn't have to evaluate y2. It can skip that whole subtree. This is like alpha-beta pruning in games.

9. ALN approximations can be quickly evaluated even when they are made up of a huge number of pieces:  The function's domain can be partitioned into blocks. Given an input vector, an algorithm quickly finds the right block by comparing components of the input vector to thresholds. It then computes only the very few affine functions affecting that block and the few max and min operators involved. What you have is a DTREE, a decision tree combined with a tiny ALN for each block. This is a huge advantage for speed of evaluation.

10. You can impose bounds on the partial derivatives of the learned function: During training, the weights of affine functions can be constrained not to go outside set bounds.  The constraint applies to all pieces and is thus global, even with smoothing. Hence even with sparse data, you can have control over ALN values between data points. If you don't have such control, your system may be unsafe!

11. ALN network architecture is built up during training, not chosen by a designer:  The ALN starts with one affine function and grows to fit the problem. This is done by splitting pieces in two when the training error is larger than the noise variance. The two new pieces are fed into a new max or min depending on the local data.

12. The learning rate for training an ALN has a clear meaning: Roughly, a learning rate of 0.2 reduces an error between what the ALN computes and the desired value by 20% at each training step. This makes it easy to automate training.

13. ALNs can produce results when the noise level is very high:  A first ALN can learn the noise variance function over the domain. This function can be used to stop splitting linear pieces of a second ALN when the mean squared training error is less than the noise variance. Overtraining is then not a problem, just as it isn't a problem with linear regression.

14. Training is stopped locally, not globally: Training is targeted at getting a good fit to data by each affine function. You don't need to end up with parts of an ALN which are overtrained or undertrained.

15. You don't need a validation set: If you use the noise variance function to stop pieces splitting, you can use all your data (except for a test set) for training. 

16. ALNs can be stored in memory: This is possible in either DTREE format or in a format that allows further training. If special parallel hardware is used, the same hardware can be used for training or evaluating many ALNs, reducing cost.

17. ALN training can be made automatic: For example, all you need is the ALNfitDeep program and properly formatted data. If you start by browsing to a tab-separated text file of numbers where the desired output is in the right-hand column, you can use the resulting trained ALN on new data to compute an output for each row. This saves on hiring expert help.

18. ALNs can distinguish helpful from unhelpful input variables: After training, an importance value is assigned to each ALN input. Unimportant input variables can be detected and removed in several steps.

19. ALNs are very appropriate for real-time operation: One successful project was control of a (simulated) vehicle active suspension system. Nothing forces training to be done batch-wise.

20. A trained ALN which is strongly monotonic in some input can be inverted: To see this, we have to look at the representation of an affine function as a hyperplane with the output weight normalized to minus 1.

 a0 (x0 - c0) + a1 (x1 - c1) + ... - (y - c(n-1)) = 0.

The ci are approximate centroids of the piece in the graph of the function. Dividing this equation by minus a1 puts the output weight, now -1, on variable x1. We do this for all the weights in the ALN. If all the a1's are negative, we are finished. If they are all positive, then we have to replace the max nodes by min nodes and vice-versa. Then x1 becomes the output of the new ALN which now has y as an input. No training is required. This has been used in control.

21. ALNs have simple properties when images are shifted, resized or have the intensity changed. To create an ALN that computes the same for a shifted image, shift the centroids except for the output centroid. To create a new ALN which gives the same output as an existing one for an image which is enlarged by factor C about the origin, divide all the weights by C and multiply all the centroids by C except for the output centroid. The output weight stays at -1.  To adjust for intensity, change the output centroid c(n-1). You can train on a standard image and it will be recognized at different sizes, light intensities, rotations, shifts... by ALNs created without further training.

22. ALN software is open-source: it can be used where cost is important, e.g. in rehabilitation. ALNs have been used to control a walking prosthesis using functional electrical stimulation.

23. Comparison of ALN output to a constant is easy: If the output of an ALN is to be tested against a threshold constant c, then one can replace the ALN by a collection of perceptrons connected to a network of boolean logic gates. e.g. max(y1, y2) > c iff (y1 > c) OR (y2 > c) and similarly for min, which converts to AND.

24. It is possible to analyze how changes of input will change the ALN response:  Say one affine function is responsible for the output. Consider the subtrees which join the path from the output to that function. Each of those has one function that determines the value of the subtree. We can determine how much change of the inputs can happen without changing the main path to the function computing the output.

25. Faster training by lazy evaluation: The program can store for each training sample, the responsible affine function. It evaluates that function first when the same training pattern is presented again. Less arithmetic evaluation may be needed because of lazy evaluation.

Conclusion: There is still plenty of low-hanging ALN fruit for researchers.  ALNs have been used to analyze tearing modes in plasmas in tokamaks -- can they be used to control tokamaks and help solve the world's energy problems? You can find the program ALNfitDeep, C++ software that works with the Community version of Visual Studio at github.com/Bill-Armstrong/Deep-Learning-ALN. 



   

  