# nca

An R implementation of the Neighborhood Components Analysis,
supporting both classification and regression problems.

Based on the [original paper](http://www.cs.nyu.edu/~roweis/papers/ncanips.pdf),
but modified in a few ways:

1. The gradient computation is vectorized, as in Python's implementation in
   [`sklearn`](https://github.com/scikit-learn/scikit-learn).
   The derivation is provided in one of the vignettes.

1. The problem is set up as a minimizing loss problem, instead of the original
   maximizing accuracy. This allows the algorithm to extend to regression tasks.
   
1. The implementation gives an optional penalty parameter, as in [NCFS](http://www.jcomputers.us/vol7/jcp0701-19.pdf).
