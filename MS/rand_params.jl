using Random
using Serialization

# p can be 28 values, from p=0.03 to p=0.30 on increments of 0.01.
# N can be 6 values, 8, 12, 16, 20, 24, 28
# θ can be 12 values, from π/48 to π/4 in increments of π/48

pinds = 1:28
Ninds = 1:6
θinds = 1:12

ordered = []

for pind in pinds
	for Nind in Ninds
		for θind in θinds
			push!(ordered, (pind, Nind, θind))
		end
	end
end

permed = shuffle(MersenneTwister(10), ordered)
println(permed)
serialize("shuf_params", permed)