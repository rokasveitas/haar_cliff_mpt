# haar_cliff_mpt

This is a repository containing all of my work from Summer 2021.

## Effect of Mølmer-Sørensen Angle on Critical Measurement Rate

My first question was extending the work in [this paper](https://arxiv.org/abs/2106.03769) by Stefanie Czischek et al.  The authors study a model of measurement-induced phase transitions in hybrid quantum circuits where the entangling two-qubit unitaries are Mølmer-Sørensen (MS) gates, defined in Eq. 1 therein.  They fix the parameter of the MS gate to be &Theta;=&pi;/4, generating maximal entanglement with every application, which allows the circuits to be simulated with fewer layers of depth before their entanglement entropy reaches its limit.  However, we hypothesize that for a different/lower value of &Theta;, the critical measurement rate p<sub>c</sub> may be lower.  This is valuable to know because high values of p are infeasible experimentally, so the experimental implementation of this circuit is often limited by the rate of measurement that must be reached before critical behavior can be observed.  By exploring the parameter space further, we can come closer to seeing a measurement-induced phase transition in experiment.

## Side-by-Side Simulation of Haar and Clifford Circuits

