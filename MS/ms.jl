using PastaQ
import PastaQ: gate, randomparams, randomlayer

using ITensors
using LinearAlgebra
using Random
using Serialization

# Defining MS and rotation to match paper conventions

gate(::GateName"MS"; θ::Number) = gate("XX"; ϕ=θ)
gate(::GateName"R"; θ::Number, φ::Number) = gate("Rn"; θ=θ, ϕ=φ-π/2, λ=-φ+π/2)

gate(::GateName"Π0") = [ 1 0
			 0 0 ]

gate(::GateName"Π1") = [ 0 0
			 0 1 ]


randomparams(::GateName"MS", args...; θ=π/4, rng=Random.GLOBAL_RNG) = (θ=θ,)  #θ=2*π * rand(rng),)  # second one is for full-random MS
randomparams(::GateName"R", args...; rng=Random.GLOBAL_RNG) = (θ=π/2, φ=rand(rng, (0, π/2, π/4) ) )


function randomlayer(
  gatename::AbstractString,
  support::Union{Vector{<:Int},Vector{Tuple},AbstractRange},
  θ::Real;
  rng=Random.GLOBAL_RNG,
  kwargs...,
)
  layer = []
  for n in support
    pars = randomparams(gatename; θ=θ, rng=rng) 
    gatepars = (
      if isempty(pars)
        (isempty(kwargs) ? nothing : values(kwargs))
      else
        merge(pars, values(kwargs))
      end
    )
    g = (isnothing(gatepars) ? (gatename, n) : (gatename, n, gatepars))
    push!(layer, g)
  end
  return layer
end

# Function to perform a projective measurement in z-basis on two neighboring qubits
# (requires the definition of two additional gates)
#
# measurement! and entanglement_entropy were given by Stef

function measurement!(ψ0::MPS, site::Int, cutoff::Real, maxdim::Real)
	ψ = normalize!(copy(ψ0))

	orthogonalize!(ψ, site)
	ρ = prime(ψ[site], tags="Site") * dag(ψ[site])

	probs = real.(diag(array(ρ)))
	#print("sum(probs: ")
	#println(sum(probs))

	σ = Int(rand() > probs[1])

	ψ = runcircuit(ψ,("Π"*"$(σ)",site), svd_alg="recursive", cutoff=cutoff, maxdim=maxdim)


	normalize!(ψ)
	ψ0[:] = ψ

end

function entanglement_entropy(ψ0::MPS, b::Int, N::Int)
		ψ = normalize!(copy(ψ0))

		orthogonalize!(ψ, b)
		_, S = svd(ψ[b], (linkind(ψ, b-1), siteind(ψ, b)))

		SvN = 0.0
		for k = 1:dim(S, 1)
			p = S[k, k]^2
			SvN -= p * log(p + 1e-20)
		end

		return SvN
end


#=
if evol, the system will measure and output its entropy after each layer
=#
function run_brick(
  state::MPS,
  N::Int,
  depth::Int,
  p=0;
  θ=9,
  twoqubitgates::Union{String,Vector{String}}="MS",
  onequbitgates::Union{Nothing,String,Vector{String}}="R",
  evol::Bool,
  rng=Random.GLOBAL_RNG,
)
	ψ = state
	if evol
		S = []
	else
		S = 0
	end

	bonds = PastaQ.lineararray(N)

	for d in 1:depth
		
		layer = []

		# one-qubit gates
		append!(layer, randomlayer(onequbitgates, N; rng=rng))

		# two-qubit gates
		append!(layer, randomlayer(twoqubitgates, bonds[2 - d%2]; θ=θ, rng=rng))

		# run unitaries
		ψ = runcircuit(ψ, layer, svd_alg="recursive") # must be recursive, lest bad things happen

	
		# projective measurements

		for b = 1:N
			if rand() < p
				measurement!(ψ, b, 1e-10, 10_000)
			end
		end

		if evol
			push!(S, entanglement_entropy(ψ, N÷2, N))
		end
	end

	if !evol
		S = entanglement_entropy(ψ, N÷2, N)
	end

	return (ψ , S)
end




function run_savg(N::Int, depth::Int, ite::Int, p::Real, evol::Bool; θ=π/4)

	s_avg = 0

	for i = 1:ite
		ψ0 = productstate(N)
		_, S = run_brick(ψ0, N, depth, p, evol=evol, θ=θ)
		
		s_avg += S
	end

	return s_avg / ite
			
end

function main(args)

	 # The goal here is to find the critical point for several different values of θ.
	 # We will do this by running a job for each value of θ and N that we're
	 # interested in, of which there will be 6 for θ and 4 for N.  The values are
	 # θ = πn/24 and N = 4* (m+1) for n=1:6 and m=1:4

	 # The argument passed will be the index of the job in Graham, so it ranges on the
	 # natural numbers as high as we want to get enough iterations.
	 #
	 # 3 days is usually enough for most of these jobs, but 4 days is better for the
	 # high-N and high-θ ones.

	 iarg = try
			parse(Int, args[1])
		catch
			throw("sadness")
	 end;

	 # p can be 28 values, from p=0.03 to p=0.30 on increments of 0.01.
	 # N can be 6 values, 8, 12, 16, 20, 24, 28
	 # θ can be 12 values, from π/48 to π/4 in increments of π/48

	 Ss = zeros(28, 6, 3)

	 depths = Vector{Int}(ones(28)) * 30
	 depths[1:2] = [100, 75]
	 
	 for i = 1:504 # = 28*6*3
	   pind = (i-1) % 28 + 1
	   Nind = (i-1) ÷ 28 % 6 + 1
	   θind = (i-1) ÷ 28 ÷ 6 % 3 + 1

	   θarg = (iarg - 1) % 4

	   p = (pind + 2) * 0.01
	   N = (Nind + 1) * 4
	   θ = (θind + 3*θarg) * π/48

	   depth = depths[pind]

	   _, Ss[pind, Nind, θind] = run_brick(productstate(N), N, depth, p, evol=false, θ=θ)

	   serialize(string("msdata/ms_", args[1]), Ss)  # this puts the Ss entropy in a 
	   						 # binary file on the cluster that
	   						 # can be analyzed later
	 end

end

main(ARGS)

