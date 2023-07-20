# Periodic Orbit Theory with Examples

The method in [`periodicOrbitTheory.m`](periodicOrbitTheory.m) implements the Periodic Orbit Theory (POT) algorithm described in the Appendix of [\[arXiv:2307.09626\]](https://arxiv.org/abs/2307.09626). Given the metadata\* of a library of orbits, it computes the POT weights, which can be used to combine orbit averages of an observable into a long-time temporal average of the same observable.

As an example, [`lorenz.m`](lorenz.m) applies the method to the [orbit metadata provided by Viswanath for the Lorenz-63 attractor](https://dept.math.lsa.umich.edu/~divakar/lorenz), and uses it to compute the temporal average of the state space speed.

\* Orbit periods, symbolic lengths, and stability information
