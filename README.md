# haar_cliff_mpt

This is a repository containing all of my work in Roger Melko's group from Summer 2021.

## Effect of Mølmer-Sørensen Angle on Critical Measurement Rate

My first question was extending the work in [this paper](https://arxiv.org/abs/2106.03769) by Stefanie Czischek et al.  The authors study a model of measurement-induced phase transitions in hybrid quantum circuits where the entangling two-qubit unitaries are Mølmer-Sørensen (MS) gates, defined in Eq. 1 therein.  They fix the parameter of the MS gate to be &Theta;=&pi;/4, generating maximal entanglement with every application, which allows the circuits to be simulated with fewer layers of depth before their entanglement entropy reaches its limit.  However, we hypothesize that for a different/lower value of &Theta;, the critical measurement rate p<sub>c</sub> may be lower.  This is valuable to know because high values of p are infeasible experimentally, so the experimental implementation of this circuit is often limited by the rate of measurement that must be reached before critical behavior can be observed.  By exploring the parameter space further, we can come closer to seeing a measurement-induced phase transition in experiment.

## Side-by-Side Simulation of Haar and Clifford Circuits

Circuits of the type examined above can be studied for many different choices of the random unitary gates.  One option is to sample uniformly from the entire two-qubit unitary group *U(4)*: this is possible due to the existence of the Haar measure, which is the unique translation-invariant probability measure for any locally-compact topological group, like *U(4)*.  This requires the full capabilities of our matrix product state (MPS) formalism for accurate simulation, which we do in PastaQ.  Another option is to sample from the Clifford group, which is the finite group that permutes the states which are "stabilized" by tensor products of single-qubit Pauli operators.  This group is interesting because it is efficiently simulable on a classical computer according to the stabilizer formalism as developed by Gottesman (2004 or whatever).

