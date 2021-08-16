
using LinearAlgebra
using Random
using QuantumClifford
using Serialization

clifset = deserialize("cliffgates")  # this is a 11520-element set containing all elements
									 # of the 2-qubit Clifford group, created by randomly 
									 # generating Clifford gates and adding them to
									 # a set until the number of elements was the order
									 # of the group.


# calculates the endpoints of all of the stabilizers in a tableau
#
# is only really meaningful when canon_clip has been applied first
function ends(s::Stabilizer)

	ns, nq = size(s)
	ends = zeros(Int, (ns, 2))

	tab = stab_to_gf2(s)

	dens = zeros(Int, (nq, 2))

	for stab = 1:ns
		for qub = 1:nq
			if tab[stab, qub] != 0 || tab[stab, qub+nq] != 0
				ends[stab, 1] = qub
				dens[qub, 1] += 1
				break
			end
		end

		for qub = nq:-1:1
			if tab[stab, qub] != 0 || tab[stab, qub+nq] != 0
				ends[stab, 2] = qub
				dens[qub, 2] += 1
				break
			end
		end
	end

	return (ends, dens)
end


# performs a row-reduction-like procedure to get a tableau into
# "clipped gauge," in which the lengths of the stabilizers are 
# meaningful
function canon_clip(s::Stabilizer)

	nq,ns = size(s)

	d = ends(s)[2]

	if sum([d[:,j] for j=1:2]) == 2 .* ones(Int, nq)
		return s
	end


	v = [paul for paul in s]


	# first, we kill everthing on the lower-left triangle

	startrow = 1

	for col = 1 : nq
		if startrow > ns - 1
			break
		end

		p1i = startrow  # the row where p1 is
		p1 = [0b0, 0b0]
		for _= startrow : ns
			p1 = [xbit(v[startrow])[col], zbit(v[startrow])[col]]
			if p1 != [0b0, 0b0]
				startrow = p1i+1
				break
			end
			insert!(v, startrow, pop!(v))
		end

		# Either we'll find a second pauli here, or the first pauli ends
		p2i = 0
		p2 = [0b0, 0b0]
		for _=startrow : ns
			p2 = [xbit(v[startrow])[col], zbit(v[startrow])[col]]
			if p2 != [0b0, 0b0] && p2 != p1
				p2i = startrow
				startrow = p2i+1
				break
			end
			insert!(v, startrow, pop!(v))
		end

		# at this point, we can tell if the first pauli ends because p2==[0b0,0b0]

		# Go down the line and fix the stabilizers whose first qubits are now non-identity
		if [p1, p2] != [[0b0, 0b0], [0b0, 0b0]]
			for row = startrow:ns

				prow = [xbit(v[row])[col], zbit(v[row])[col]]

				if prow == p1
					v[row] = v[row] * v[p1i]
				elseif prow == p2 && p2 != [0b0, 0b0]
					v[row] = v[row] * v[p2i]
				elseif prow != [0b0, 0b0] && p2i != 0
					v[row] = v[row] * v[p1i] * v[p2i]
				end
			end
		end
	end

	startrow = ns
	usedinds = Set{Int}()

	for col = nq:-1:1

		p1 = [0b0, 0b0]
		p1i = 0

		p2 = [0b0, 0b0]
		p2i = 0
		for row = startrow:-1:1

			prow = [xbit(v[row])[col], zbit(v[row])[col]]
			if prow != [0b0, 0b0]
				if p1i == 0 && !in(row, usedinds)
					p1i = row
					p1 = prow
					push!(usedinds, p1i)
				elseif p1 == prow
					v[row] = v[p1i] * v[row]
				elseif p2i == 0 && !in(row, usedinds)
					p2i = row
					p2 = prow
					push!(usedinds, p2i)
				elseif p2 == prow
					v[row] = v[p2i] * v[row]
				elseif p1i != 0 && p2i != 0
					v[row] = v[p1i] * v[p2i] * v[row]
				end
			end
		end
	end


	sout = Stabilizer(v)

	return sout
end

# computes the entanglement entropy of a contiguous region starting at fir and ending at las
function ent(s::Stabilizer, fir::Int, las::Int; canon=true)

	nq, ns = size(s)

	d = ends(s)[2]

	if canon && sum([d[:,j] for j=1:2]) != 2 .* ones(Int, nq) # checks if the tableau is already clipped
		s = canon_clip(s)
	end

	es = ends(s)[1]

	crossings = 0

	for stabi = 1:size(es)[1]
		e = es[stabi,:]
		if e[1] < fir && e[2] <= las && e[2] >= fir
			crossings += 1
		elseif e[2] > las && e[1] >= fir && e[1] <= las
			crossings += 1
		end
	end

	return crossings/2
end

# permutes the columns of a tableau, which was necessary because QuantumClifford only has colpermute!
function colpermute(s::Stabilizer, perm)
	r = copy(s)
	return colpermute!(r, perm)
end


# computes the mutual information of s between A and B, where A and B are either vectors describing
# not-necessarily-contiguous subsets of the qubits or just single qubits
function mutual_inf(s::Stabilizer, A::Union{Vector{Int}, Int}, B::Union{Vector{Int}, Int})

	if typeof(A) <: Int
		A = [A]
	end 

	if typeof(B) <: Int
		B = [B]
	end

	if A == B
		return 0
	end

	Al = size(A)[1]
	Bl = size(B)[1]

	perm = vcat(A, B, [i for i=1:s.nqubits if !in(i, A) && !in(i, B)])
	r = colpermute(s, perm)

	AB = ent(r, 1, Al + Bl)
	A = ent(r, 1, Al)
	B = ent(r, Al + 1, Al + Bl)

	return A + B - AB
end


# implements a measure that Tim suggested, which is not used or relevant but was interesting
function tim_meas(s::Stabilizer, a::Int, b::Int; print=false)

	sc = copy(s)
	nq = size(sc)[1]

	zs = one(Stabilizer, nq)

	

	for i=1:nq
		if i != a && i != b
			sc, antiind, result = project!(sc, zs[i])
			if isnothing(result)
				sc.phases[antiind] = rand([0x0, 0x2])
			end
		end
	end

	sc = canon_clip(sc)

	if print
		println(sc)
	end

	if ends(sc)[2][:,1] == ones(nq)
		return 0
	end
	return 1

end

# runs a brick-wall simulation with Clifford gates, returning 
# whatever quantity is calculated by f
function run_brick_clif(N::Int,
				   		depth::Int;
				   		p=0, 
				   		evol=false,
				 		f = x -> 0)

	s = one(Stabilizer, N)

	nq, ns = size(s)

	#ents = []

	bonds = [[[j, j+1] for j=1:2:N-1-N%2],
			 [[j, j+1] for j=2:2:N-2+N%2]]

	zs = one(Stabilizer, N)

	for d in 1:depth

		for bond in bonds[(d-1)%2+1]
			apply!(s, rand(clifset), bond)
		end

		for i=1:N
			if rand() < p
				s, antiind, result = project!(s, zs[i])
				if isnothing(result)
					s.phases[antiind] = rand([0x0, 0x2])
				end
			end
		end	

		if evol
			push!(ents, f(s))
		end

	end

	s = canon_clip(s)

	if !evol
		ents = f(s)
	end

	return (s, ents)
end

# computes the average value of f over "it" iterations
# 
# also works on vector-valued f.
function clif_avg_S(N::Int,
			   depth::Int,
			   it::Int;
			   p=0, 
			   evol=false,
			   f = x -> 0)
	
	siz = f(one(Stabilizer, N))

	if typeof(siz) <: Number
		if evol
			ents = zeros(depth)
		else
			ents = 0
		end
	else
		siz = size(siz)
		siz = evol ? Tuple(push!([i for i in siz], depth)) : siz
		ents = zeros(siz)
	end

	for i=1:it
		if typeof(ents) <: Number
			ents += run_brick(N, depth, p=p, evol=evol, f=f)[2]
		else
			ents .+= run_brick(N, depth, p=p, evol=evol, f=f)[2]
		end

		# if i%100 == 0
		# 	println(i)
		# end
	end

	if typeof(ents) <: Number
		ents /= it
	else
		ents ./= it
	end

	return ents
end


# finds the power law dependence of y on x by least-squares regression
# 
# this is used to find the mutual information fall-off power with the qubit separation
function find_pow(x::Vector, y::Vector)
	lx = copy(x)
	ly = copy(y)

	n = length(x)
	for i in 1:n
		if ly[i] == 0
			lx[i] = 0
		end

		if lx[i] == 0
			ly[i] = 0
		end
	end

	filter!(a -> a!=0, lx)
	filter!(a -> a!=0, ly)

	println(lx)
	println(ly)

	n = length(lx)

	if length(ly) != n
		throw("bad dims")
	end

	lx = log.(lx)
	ly = log.(ly)

	beta = (n * dot(lx, ly) - sum(lx) * sum(ly)) / (n * dot(lx, lx) - sum(lx)^2)
	alph = (sum(ly) - beta * sum(lx)) / n

	return (pow=beta, coef=exp(alph), n=n)
end
