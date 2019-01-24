# Why the ideas in ALNfitDeep are important for reinforcement learning

An **ALN** is a neural network for deep learning that has a very simple functional form.  It is composed of **affine functions** defining planes, or in higher dimensions, hyperplanes, at the inputs and **maximum** and **minimum** operators in all other layers. So the graph of the function learned is composed of flat pieces joined to form a continuous surface.

Suppose now the function to be learned from samples indicates the value of a state resulting from a sequence of actions each of which is associated with a reinforcement (or punishment). We ask: "Can one change the inputs (from effectors) to improve the predicted value by changing the actions (effector outputs) involved?" Since any input vector to an ALN selects one flat part of the function surface, the one used to compute the output, a path towards greater predicted value is obtained by examining the coefficients of that flat piece. Following that path, perhaps along several joined flat pieces of the function surface, the inputs, e.g. proposed effector positions of a robot, can be changed to increase the value.

We are a long way from the final realization of this method, but the deep learning method used in ALNfitDeep is a start. A first step is developing techniques to make training even faster.  Faster training and execution will support very large ALNs that are necessary in complicated problems.

Ideas like **alpha-beta pruning** are used already in the program.

**Training could be localized** like evaluation already is: a DTREE idea partitions the input space so huge ALNs can be evaluated by quickly selecting one of a collection of much smaller ALNs to be evaluated. Locality and speed of computation is a distinct benefit of this kind of multilayer perceptron over the ones using standard back-propagation where a huge number of parameters can influence an output.

Further increases in speed can be obtained by **limiting the number of variables entering into each affine function**.  The "importance" values computed in ALNfitDeep select input variables for all affine functions in the whole ALN.  In general, different parts of a huge ALN could depend on different input variables. This is not the case for ALNfitDeep.

Several facts point to this plan for advancing reinforcement learning being realizable:

1. **ALNs can be grown to fit any continuous function** to any desired accuracy uniformly on any compact domain.

2. If an ALN is prescribed to be strongly monotonic in any input variable, then without further training, **an ALN can be constructed to compute the inverse function** by simple changes of coefficients and maximum and minimum operators.

3. **Value iteration and Q-learning have been shown to work** on problems like the inverted pendulum. This is probably due to the fact that the value function, for example, can be as "fine-grained" as desired, so Bellman's equation, where a function value is determined by those close to it, can be solved. Having an architecture which can grow without limit is important here.

wwarmstrong    AT    gmail.com 