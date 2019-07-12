# 30 advantages of using ALNs for machine learning.

What is an ALN? An ALN is a type of neural net consisting of a first layer of affine functions y = a0*x0 + a1*x1 ... + an * xn + c whose outputs are combined in all other layers by two-input maximum and minimum functions to produce a single output. If you spend a few seconds thinking about the following, you will understand ALNs.

1. ALN structure is very simple: Because all arithmetic is in the first layer, ALNs can be extremely deep without significant loss of speed.

2. ALNs are easy to understand and analyze: The output function is continuous. The graph is formed by linear pieces, parts of of the hyperplanes defined by the coefficients ("weights") of the affine functions. Pieces are bounded by other pieces where they intersect.The ALN is not a "black box". 

3. ALNs can be smoothed: Inserting small quadratic or quartic "fillets" at max and min nodes makes the output once or twice continously differentiable.

4. Even smoothed ALN outputs are easy to analyze: The smoothed ALN can be made to differ everywhere from the piecewise linear surface by at most a prescribed small amount.

5. ALNs elegantly solve the credit assignment problem: It is easy to determine how an ALN output value is computed, starting at the output. At a max (or min) node, the input subtree with the greater (or lesser) value is responsible. Ultimately, one affine function in the first layer is responsible for that output. It is the one whose weights are updated during training. (Ties at intersections of linear pieces are very rare.)

6. ALN supervised training is simple: When an affine function computing the output value is determined, an iterative form of least squares adjusts those weights to reduce the error.

7. ALNs can capture convexity: If there are only max nodes in an ALN, the surface is convex down. Min nodes make it convex-up. The structure of an efficient ALN tree would follow the convexities of the ideal function. 

8. Universal approximation: An ALN can fit any continuous function on a compact set to any degree of accuracy. For example, just take the maximum of well-chosen minima (e.g. the tops of truncated pyramids). 

9. ALN evaluation is lazy: for example if you have a part of the ALN that looks like max(min(5, y),6), then the system doesn't have to evaluate y to get the answer 6. It can skip the whole subtree that computes y, speeding the computation!

10. ALN approximations can be very rapidly evaluated even when they are made up of a huge number of pieces:  The function's domain can be partitioned into blocks. Given an input vector, an algorithm quickly finds the right block by comparing components of the input vector to thresholds. It then computes only the very few affine functions affecting that block and the few max and min operators involved. What you have is a DTREE, a decision tree combined with a tiny ALN for each block.

11. You can impose bounds on the weights (partial derivatives) of the learned function: During training, the weights of affine functions can be constrained to remain in set bounds.  The constraint applies to all pieces and is thus global, even with smoothing. Hence even with sparse data, you can have control over ALN values between data points. If you don't have such control, your system may be unsafe!

12. ALN network architecture is built up during training, not chosen by a designer:  The ALN starts with one affine function and grows to fit the problem. This is done by splitting linear pieces in two when, after several passes through the training data, the square error is still larger than the noise variance. The two new pieces are fed into a new max or min depending  on a measure of convexity of the local data.

13. The learning rate for training an ALN has a clear meaning: Roughly, a learning rate of 0.2 reduces an error between what the ALN computes and the desired value for the trained sample by 20%. The error for the other samples may increas. Still 20 passes through the data (an epoch) tends to reduce the error significanly so splitting can occur. 

14. ALNs can produce results when the noise level is very high: Overtraining is prevented by measuring the noise variance in various parts of the domain and stopping training on parts of the domain by doing an F-test on mean square training error and noise variance of the linear pieces. 

15. Training is stopped locally, not globally: Training is targeted at getting a good fit to data by each affine function. You don't need to end up with parts of an ALN which are overtrained or undertrained.

16. You don't need a validation set: If you use the noise variance function to stop pieces splitting, you can use all your data (except for a test set) for training. 

17. ALNs can be trained and evaluated efficiently on a computer without special hardware mainly due to lazy evaluation. 

18. ALNs can distinguish helpful from unhelpful input variables: After training, an importance value can be assigned to each ALN input. Unimportant input variables can be detected and removed.

19. ALNs are very appropriate for real-time operation: One successful project was control of a (simulated) vehicle active suspension system.  Even noise variance estimation on a stream of inputs can be done in time linear in training set size.

20. A trained ALN which is strongly monotonic in some input can be inverted without training: To see this, we have to look at the representation of an affine function as a hyperplane with the output weight normalized to minus 1.  a0 (x0 - c0) + a1 (x1 - c1) + ... - (y - c(n-1)) = 0.  The ci are approximate centroids of the piece in the graph of the function. Dividing this equation by minus a1 puts the output weight, now -1, on variable x1. We do this for all the weights in the ALN. If all the a1's are negative, we are finished. If they are all positive, then we have to replace the max nodes by min nodes and vice-versa. Then x1 becomes the output of the new ALN which now has y as an input.

21. ALNs have simple properties when images are shifted, resized or have the intensity changed. To create an ALN that computes the same for a shifted image, shift the centroids except for the output centroid. To create a new ALN which gives the same output as an existing one for an image which is enlarged by factor C about the origin, we change the weights systematically. You can train on a standard image and it will be recognized at different sizes, light intensities, rotations, shifts... by ALNs created without further training.

22. ALN software is open-source: it can be used where cost is important, e.g. in rehabilitation. ALNs have been used to control a walking prosthesis using functional electrical stimulation.

23. Comparison of the ALN output to a constant c is done by using a collection of perceptrons at the input layer connected to a network of boolean logic gates. e.g. max(y1, y2) > c iff (y1 > c) OR (y2 > c) and similarly for min, which converts to AND.

24. It is possible to analyze how changes of input will change the ALN response:  Say one affine function is responsible for the output. Consider the subtrees which join the path from the output to that function. Each of those subtrees has one function that determines the value of the subtree. We can determine how much change of the inputs can happen without changing the
value on the main path. (Try this with an ordinary neural net!)

25. Faster training by lazy evaluation: The program can store for each training sample, the responsible affine function. It evaluates that function first when the same training pattern is presented again. Less arithmetic evaluation may be needed because of lazy evaluation.

26. ALNs have been shown capable of doing reinforcement learning. When a function being learned by a neural network is trained using values of the same function, as solving Bellman's equation requires, the complexity of the function is not known in advance. ALN architecture grows to fit the problem.

27. It is possible to check all inputs that lead to an output value Y being greater than some constant C: Convert to a boolean logic tree (see 23 above) and then make that into a two-layer tree which is a max of mins. At the inputs there are now perceptrons of the form f>C. These are half-spaces which intersect to form convex polyhedra.  Their union is the set of all inputs satisfying Y>=C. You can do the same for Y<= D and find out where Y is in the interval [C,D]. It may not be easy, just possible.

28. You don't have to choose a squashing function: The non-linearity needed is in the max and min operators.  If you hold either input to a max or min fixed, the rest behaves like a rectified linear unit (ReLU). The whole ALN appears linear if you just add input vectors or multiply them by positive numbers. When you subtract them or multiply by negative numbers, linearity breaks down, as it must.

29. You don't have to normalize the data or regularize anything.

30. There is no problem of vanishing gradient as with networks depending on the chain rule. Max and min nodes do no multiplications.
  
Conclusion: There is still plenty of low-hanging ALN fruit for researchers.  ALNs have been used to analyze tearing modes in plasmas in tokamaks -- can they be used to control tokamaks and help solve the world's energy problems? You can find the program ALNfitDeep, C++ software that works with the Community version of Visual Studio at github.com/Bill-Armstrong/Deep-Learning-ALN. 



   

  