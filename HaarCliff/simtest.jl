#=

hc_run doesn't work because my cmset isn't correct, but
I was using this to try to do a simulataneous test: I 
applied random Clifford gates to a stabilizer while applying
those same gates to an MPS and wanted to verify that their
resulting mutual informations were the same.  I ran som tests
manually in shell, applying known gates repeatedly, and I
was convinced that everything is working well except for cmset.

=#


include("cms.jl")   # QuantumClifford
include("clif1.jl") # LinearAlgebra, Random, QuantumClifford, Serialization


# ends(s::Stabilizer)
# canon_clip(s::Stabilizer)
# ent(s::Stabilizer, fir::Int, las::Int; canon=true)
# mutual_inf(s::Stabilizer, A::Union{Vector{Int}, Int}, B::Union{Vector{Int}, Int})
# run_brick(N::Int,
# 			depth::Int;
# 			p=0, 
# 			evol=false,
# 			f = x -> 0)
# avg_S(N::Int,
# 			depth::Int,
# 			it::Int;
# 			p=0, 
# 			evol=false,
# 			f = x -> 0)
# find_pow(x::Vector, y::Vector)

include("haar1.jl") # LinearAlgebra, Random, Statistics, ITensors, PastaQ 

# mutual_inf(psi0::MPS, a::Int, b::Int)
# stef_mutual_inf(psi0::MPS, aval::Int, bmax::Int, M::Int)
# measurement!(ψ0::MPS, site::Int, cutoff::Real, maxdim::Real)
# entanglement_entropy(ψ0::MPS, b::Int, N::Int)
# randomlayer(
#   gatename::AbstractString,
#   support::Union{Vector{<:Int},Vector{Tuple},AbstractRange};
#   rng=Random.GLOBAL_RNG,
#   kwargs...,
# )
# run_brick(
#   state::MPS,
#   N::Int,
#   depth::Int,
#   p=0;
#   twoqubitgates::Union{String,Vector{String}}="randU",
#   onequbitgates::Union{Nothing,String,Vector{String}}=nothing,
#   evol::Bool=false,
#   rng=Random.GLOBAL_RNG,
# )
# run_savg(N::Int, depth::Int, ite::Int, p::Real, evol::Bool)

using .CMs

cmset = deserialize("cmset")

gate(::GateName"a"; a::Matrix) = a  # PastaQ gate name allowing me to pass arbitrary matrices


function mutual_inf(a::Any)
	n = size(a)[1]

	Ss = zeros(n, n)

	for i=1:n, j=i+1:n
		if i != j
			Ss[i,j] = mutual_inf(a, i, j)
		end
	end

	return Ss
end

function gatepars(mat::Matrix, bond1::Int, bond2::Int)
	return ("a", (bond1, bond2), (a=mat,))
end


function hc_run(N::Int, depth::Int, p::Real)

	Ss = []

	ψ = productstate(N)
	s = one(Stabilizer, N)

	bonds = [[[j, j+1] for j=1:2:N-1-N%2],
			 [[j, j+1] for j=2:2:N-2+N%2]]

	zs = one(Stabilizer, N)

	for d in 1:depth

		layer = []

		# unitaries
		for bond in bonds[(d-1)%2+1]
			gate = rand(cmset)

			#println(gate)

			apply!(s, gate.c, bond) #clifford

			push!(layer, gatepars(gate.m, bond[1], bond[2])) #haar layer
		end

		ψ = runcircuit(ψ, layer, svd_alg="recursive") #haar run

		# measurements

		for b = 1:N
			r = rand()

			if r < p
				phase = measurement!(ψ, b, 1e-10, 10_000) # haar

				# clifford
				s, antiind, result = project!(s, zs[b])
				if isnothing(result)
					s.phases[antiind] = phase
				end
			end
		end



		push!(Ss, [s, ψ])
	end

	return [Ss[i][j] for i=1:depth, j=1:2]
end










