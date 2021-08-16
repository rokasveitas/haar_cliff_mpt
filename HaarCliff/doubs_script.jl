using QuantumClifford
using Serialization
using Random

include("cms.jl")
using .CMs

#=

sings is the set of every single-qubit Clifford operator and its
size is 24.  doubs is the set of every two-qubit non-entangling
non-permuting Clifford operator, of which there are 24^2 / 4 = 144
(the factor of 4 is because only the product of the phases is preserved
by the tensor product.)  d is a dictionary between the stabilizer 
forms and matrix forms of these operators because gates_mats.jl needs them.
One could write something nice that would more directly compute the
required matrix to get a given two-qubit operator, but I found this easier
and it's faster.

=#


sings = Set()
#singsa = []

for i = 1:32
	a = CM(CliffordId, i2m)
	if i%2==1 (a *= z) end
	if i÷2%2==1  (a *= x) end
	if i÷4%2==1  (a *= ph) end
	if i÷8%2==1  (a *= h) end
	if i÷16%2==1  (a *= ph) end
	push!(sings, a)
	#push!(singsa, a)
end

# println(length(sings))
# 24


doubs = Set()

for i in sings
	for j in sings
		push!(doubs, i ⊗ j)
	end
end

# println(length(doubs))
# 144

d = Dict()
for a in doubs
	push!(d, a.c => a.m)
end

serialize("doubdict", d)