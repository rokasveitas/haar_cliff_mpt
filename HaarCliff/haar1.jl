
using LinearAlgebra
using Random
using Statistics
using ITensors
using PastaQ
import PastaQ: gate, randomlayer


# two-qubit reduced density matrix portion translated from 
# http://www.itensor.org/docs.cgi?vers=cppv3&page=formulas/mps_two_rdm

function mutual_inf(psi0::MPS, a::Int, b::Int)
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

	# Now we have the density matrix rho.  rho has indices sa, sb, sa', sb'.
	# Next, compute the AB entanglement entropy.

	sa = siteind(psi, a)
	sb = siteind(psi, b)

	U, S, V = svd(rho, sa, sb)

	AB = 0.0
	for i=1:dim(S, 1)
		p = S[i,i]^2
		AB -= p*log(p + 1e-20)
	end
	# println("AB: ")
	# println(AB)

	# Now compute individual A and B entanglement entropies.

	U, S, V = isnothing(la) ? svd(psi[a], sa) : svd(psi[a], (la, sa))
	A = 0.0
	for i=1:dim(S, 1)
		p = S[i,i]^2
		A -= p*log(p + 1e-20)
	end

	# println("A: ")
	# println(A)

	orthogonalize!(psi, b)

	U, S, V = isnothing(lb) ? svd(psi[b], sb) : svd(psi[b], (lb, sb))
	B = 0.0
	for i=1:dim(S, 1)
		p = S[i,i]^2
		B -= p*log(p + 1e-20)
	end
	# println("B: ")
	# println(B)

	return (A + B - AB) * 2 / log(2)
end


function stef_mutual_inf(psi0::MPS, aval::Int, bmax::Int, M::Int)
	psi = normalize!(copy(psi0))
	N = length(psi)

    # array to store information
	I = zeros(4)

    # draw samples (four sets since we want I_A, I_B and I_{A,B})
	samps_1 = getsamples(psi, M, local_basis=["Z"])
	samps_2 = getsamples(psi, M, local_basis=["Z"])
	samp_data_1 = PastaQ.convertdatapoints(copy(samps_1); state=true)
	samp_data_2 = PastaQ.convertdatapoints(copy(samps_2); state=true)

    # evaluate entropies
	a = [aval]
	for (count, bval) in enumerate(aval+9:3:bmax)
		b = [bval]

		ab = zeros(2)
		ab[1] = aval
		ab[2] = bval

        # switch subsets in copies of system
		samp_data_a1_switch = copy(samp_data_1)
		samp_data_a2_switch = copy(samp_data_2)
		samp_data_b1_switch = copy(samp_data_1)
		samp_data_b2_switch = copy(samp_data_2)
		samp_data_ab1_switch = copy(samp_data_1)
		samp_data_ab2_switch = copy(samp_data_2)

		s = siteinds(psi)

		trace_rho_a = 0.
		trace_rho_b = 0.
		trace_rho_ab = 0.

        # get traces of density matrices from all samples
		for m in 1:M
			x1 = samp_data_1[m,:]
			x2 = samp_data_2[m,:]

			for n in a
				samp_data_a1_switch[m,Int(n)] = samp_data_2[m,Int(n)]
				samp_data_a2_switch[m,Int(n)] = samp_data_1[m,Int(n)]	
			end

			for n in b
				samp_data_b1_switch[m,Int(n)] = samp_data_2[m,Int(n)]
				samp_data_b2_switch[m,Int(n)] = samp_data_1[m,Int(n)]
			end

			for n in ab
				samp_data_ab1_switch[m,Int(n)] = samp_data_2[m,Int(n)]
				samp_data_ab2_switch[m,Int(n)] = samp_data_1[m,Int(n)]
			end

			xa1_switch = samp_data_a1_switch[m,:]
			xa2_switch = samp_data_a2_switch[m,:]
			xb1_switch = samp_data_b1_switch[m,:]
			xb2_switch = samp_data_b2_switch[m,:]
			xab1_switch = samp_data_ab1_switch[m,:]
			xab2_switch = samp_data_ab2_switch[m,:]

			psi_1 = dag(psi[1]) * state(x1[1], s[1])
			psi_2 = dag(psi[1]) * state(x2[1], s[1])
			psi_a1_switch = dag(psi[1]) * state(xa1_switch[1], s[1])
			psi_a2_switch = dag(psi[1]) * state(xa2_switch[1], s[1])
			psi_b1_switch = dag(psi[1]) * state(xb1_switch[1], s[1])
			psi_b2_switch = dag(psi[1]) * state(xb2_switch[1], s[1])
			psi_ab1_switch = dag(psi[1]) * state(xab1_switch[1], s[1])
			psi_ab2_switch = dag(psi[1]) * state(xab2_switch[1], s[1])

			for j in 2:N
				psi_1_r = dag(psi[j]) * state(x1[j], s[j])
				psi_2_r = dag(psi[j]) * state(x2[j], s[j])
				psi_a1_switch_r = dag(psi[j]) * state(xa1_switch[j], s[j])
				psi_a2_switch_r = dag(psi[j]) * state(xa2_switch[j], s[j])
				psi_b1_switch_r = dag(psi[j]) * state(xb1_switch[j], s[j])
				psi_b2_switch_r = dag(psi[j]) * state(xb2_switch[j], s[j])
				psi_ab1_switch_r = dag(psi[j]) * state(xab1_switch[j], s[j])
				psi_ab2_switch_r = dag(psi[j]) * state(xab2_switch[j], s[j])

				psi_1 = psi_1 * psi_1_r
				psi_2 = psi_2 * psi_2_r
				psi_a1_switch = psi_a1_switch * psi_a1_switch_r
				psi_a2_switch = psi_a2_switch * psi_a2_switch_r
				psi_b1_switch = psi_b1_switch * psi_b1_switch_r
				psi_b2_switch = psi_b2_switch * psi_b2_switch_r
				psi_ab1_switch = psi_ab1_switch * psi_ab1_switch_r
				psi_ab2_switch = psi_ab2_switch * psi_ab2_switch_r
			end

			trace_rho_a += (psi_a1_switch[] / psi_1[]) * (psi_a2_switch[] / psi_2[])
			trace_rho_b += (psi_b1_switch[] / psi_1[]) * (psi_b2_switch[] / psi_2[])
			trace_rho_ab += (psi_ab1_switch[] / psi_1[]) * (psi_ab2_switch[] / psi_2[])
		end

		trace_rho_a /= M
		trace_rho_b /= M
		trace_rho_ab /= M

		entropy_a = - real(log2(trace_rho_a))
		entropy_b = - real(log2(trace_rho_b))
		entropy_ab = - real(log2(trace_rho_ab))

		I[count] = entropy_a + entropy_b - entropy_ab

		println(ab)
	end

	return I
end

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



gate(::GateName"Π0") = [ 1 0
					     0 0 ]

gate(::GateName"Π1") = [ 0 0
						 0 1 ]

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
			push!(S, entanglement_entropy(ψ, N÷2, N))
		end
	end

	if !evol
		S = entanglement_entropy(ψ, N÷2, N)
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
		_, S = run_brick(ψ0, N, depth, p, evol=evol)
		
		s_avg += S
	end

	s_avg ./= ite
	return s_avg
			
end