HMM
===

A C++ implementation of various order discrete Markov models for sequence analysis. 
HMM is coded in a generic programming style using C++ templates for efficiency
with a python interface. HMM implements the same algorithms on top of two distinct data mdoels:
one for dense hidden Markov models where transitions between most states are possible and another
for sparse HMMs where a limited number of transitions are possible. The library provides code to convert
between the models.
