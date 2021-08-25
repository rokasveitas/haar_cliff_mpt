# Mutual information power law in Haar and Clifford systems

The main files of relevance here are `haar1.jl` and `clif1.jl`: they contain all of the functions used for Haar and Clifford computation, respectively.  The Haar calculations use the `PastaQ.jl` module, based off of `ITensors`, to perform matrix product state simulations, and the Clifford calculations are done in `QuantumClifford.jl`.

Currently, the Haar calculations either have a bug or haven't been run with enough iterations: only time will tell.

The Clifford results look like this:
![image](https://user-images.githubusercontent.com/5233686/130823038-98f6730f-8a82-4ee7-9e7f-19b0271a01af.png)
They are close to a -4 power law, but actually much more like a -3.
