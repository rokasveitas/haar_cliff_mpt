# Mutual information power law in Haar and Clifford systems

The main files of relevance here are `haar1.jl` and `clif1.jl`: they contain all of the functions used for Haar and Clifford computation, respectively.  The Haar calculations use the `PastaQ.jl` module, based off of `ITensors`, to perform matrix product state simulations, and the Clifford calculations are done in `QuantumClifford.jl`.

The Clifford results look like this:
![image](https://user-images.githubusercontent.com/5233686/130823038-98f6730f-8a82-4ee7-9e7f-19b0271a01af.png)
They are close to a -4 power law, but actually much more like a -3.  The power law was extracted with `find_pow`, which is a method in `clif1.jl` but would be equally applicable to Haar systems.  The data for this is in the file `clif_mutinf_tup`, which is a serialized tuple containing the output data and also all of the other parameters that defined the run (`N` values, number of iterations, etc.)

I have not had enough time to run the Haar calculations with enough iterations to get any good results—`haar1.jl` had a bad bug until late August, and then I started school and my Graham queue was full of Mølmer-Sørensen stuff.  But all of the machinery is there to get the data: it just has to be run.

Most of the files in this directory were created in order to make a dictionary of two-qubit Clifford operators and their correspoding matrix representation on the Hilbert space, the goal of which was to test that the mutual information calculations between `PastaQ.jl` and `QuantumClifford` would agree.  This took probably too long and didn't come to fruition because we decided that the test was not really necessary, and we were convinced by other tests that the mutual information for both cases was functional.
