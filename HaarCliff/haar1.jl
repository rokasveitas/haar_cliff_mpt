
using LinearAlgebra
using Random
using Statistics
using ITensors
using Serialization
using PastaQ
using Plots
import PastaQ: gate, randomlayer
using Dates


# two-qubit reduced density matrix portion translated from 
# http://www.itensor.org/docs.cgi?vers=cppv3&page=formulas/mps_two_rdm

function rho_two(psi0::MPS, a::Int, b::Int)
	psi = normalize!(copy(psi0))
	N = length(psi)

	# First: find the two-qubit density matrix for qubits a and b

	# shift orthogonality center to a
	orthogonalize!(psi, a)

	psidag = dag(psi)
	prime!(psidag, "Link")

	# index linkin a to a-1
	la = linkind(psi, a-1)

	rho = isnothing(la) ? psi[a] : prime(psi[a], la)

	rho *= prime(psidag[a], "Site")

	for k = a+1 : b-1
		rho *= psi[k]
		rho *= psidag[k]
	end

	lb = linkind(psi, b)

	rho *= isnothing(lb) ? psi[b] : prime(psi[b], lb)

	rho *= prime(psidag[b], "Site")

	return rho
end

function mutual_inf(psi0::MPS, a::Int, b::Int)
	psi = normalize!(copy(psi0))
	rho = rho_two(psi, a, b)

	# Now we have the density matrix rho.  rho has indices sa, sb, sa', sb'.
	# Next, compute the AB entanglement entropy.

	sa = siteind(psi, a)
	sb = siteind(psi, b)
	la = linkind(psi, a-1)
	lb = linkind(psi, b)

	U, S, V = svd(rho, sa, sb)

	AB = 0.0
	for i=1:dim(S, 1)
		p = S[i,i]^2
		AB -= p*log2(p + 1e-20)
	end
	# println("AB: ")
	# println(AB)

	# Now compute individual A and B entanglement entropies.

	orthogonalize!(psi, a)

	U, S, V = isnothing(la) ? svd(psi[a], sa) : svd(psi[a], (la, sa))
	A = 0.0
	for i=1:dim(S, 1)
		p = S[i,i]^2
		A -= p*log2(p + 1e-20)
	end

	# println("A: ")
	# println(A)

	orthogonalize!(psi, b)

	U, S, V = isnothing(lb) ? svd(psi[b], sb) : svd(psi[b], (lb, sb))
	B = 0.0
	for i=1:dim(S, 1)
		p = S[i,i]^2
		B -= p*log2(p + 1e-20)
	end
	# println("B: ")
	# println(B)

	return (A + B - AB)
end

gate(::GateName"Π0") = [ 1 0
					     0 0 ]

gate(::GateName"Π1") = [ 0 0
						 0 1 ]

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

	return σ == 0 ? 0x0 : 0x2

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

function randomlayer(
  gatename::AbstractString,
  support::Union{Vector{<:Int},Vector{Tuple},AbstractRange};
  rng=Random.GLOBAL_RNG,
  kwargs...,
)
  layer = []
  for n in support
    pars = randomparams(gatename, 4; rng=rng) # 4 is because it's a 4x4 matrix
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


function run_brick_haar(
  state::MPS,
  N::Int,
  depth::Int,
  p=0;
  twoqubitgates::Union{String,Vector{String}}="randU",
  onequbitgates::Union{Nothing,String,Vector{String}}=nothing,
  evol::Bool=false,
  rng=Random.GLOBAL_RNG,
)
	ψ = normalize!(copy(state))
	N = length(ψ)
	S = evol ? [] : 0

	bonds = PastaQ.lineararray(N)

	for d in 1:depth
		
		layer = []

		# one-qubit gates
		if !isnothing(onequbitgates)
			append!(layer, randomlayer(onequbitgates, N; rng=rng))
		end

		# two-qubit gates
		append!(layer, randomlayer(twoqubitgates, bonds[2 - d%2]; rng=rng))

		# run unitaries
		ψ = runcircuit(ψ, layer, svd_alg="recursive")

	
		# projective measurements

		for b = 1:N
			if rand() < p
				measurement!(ψ, b, 1e-10, 10_000)
			end
		end

		if evol && d%2 == 1
			push!(S, mutual_inf(ψ, N÷3, 2*N÷3))
		end
	end

	if !evol
		S = mutual_inf(ψ, N÷3, 2*N÷3)
	end
	#println(S)
	#println(typeof(S))
	return (ψ , S)
end

function haar_avg_S(N::Int, depth::Int, ite::Int, p::Real, evol::Bool)
	#global savg = zeros(depth)
	#println(savg)
	#ents = []

	s_avg = evol ? zeros(depth÷2) : 0

	for i = 1:ite
		ψ0 = productstate(N)
		_, S = run_brick_haar(ψ0, N, depth, p, evol=evol)
		
		if size(s_avg) != size(S)
			s_avg = zeros(size(S))
		end

		s_avg += S
	end

	if evol
		return s_avg ./ ite
	end

	return s_avg / ite
			
end


#if abspath(PROGRAM_FILE) == @__FILE__ # if run directly from shell

function main()
	beg = now()
	

	Ss = zeros(6, 4)
	times = zeros(Millisecond, (6, 4))

	las = now()

	for N = 8:4:20
		for p = 0.1:0.04:0.3
			iN = N÷4 - 1
			ip = Int(round(p/0.04 - 0.5)) - 1

			Ss[ip, iN] = haar_avg_S(N, 2*N, 200, p, false)

			this = now()

			times[ip, iN] = this - las

			println((iN, ip, times[ip, iN], Ss[ip, iN]))

			las = this
		end
	end

	fin = now()

	print("Total: ")
	println(fin - beg)
	serialize("haar_s_3", (Ss, times))
	# run(`say "all done"`)
	println(Ss)
	
	# show(plot(Ss))
end

main()












