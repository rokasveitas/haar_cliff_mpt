#= 
A CM is a composite type that has fields for a 2-qubit Clifford
    operator in the stabilizer representation and for a
    representation of the same operator as a 4x4 complex matrix
    acting on the Hilbert space of the 2-qubit system.  The is
    useful because we would like to replicate a given random
    Clifford circuit in the MPS formalism, so we need to be able
    to track the matrix representations of all 11520 2-qubit Clifford
    operators.

=#

module CMs

using QuantumClifford
using QuantumClifford: getallpaulis
import QuantumClifford: ⊗
export CM,
       id, cx, cxx, cz, swap,
       z, x, ph, h,
       i2m, QuantumClifford,
       i1, i2

struct CM
    c::CliffordOperator
    m::Matrix{ComplexF64}
end

Base.:(==)(a::CM, b::CM) = a.c == b.c

function ⊗(a::CM, b::CM)
    return CM(a.c ⊗ b.c, kron(a.m, b.m))
end

function ⊗(a::CliffordOperator, b::CliffordOperator)
    sizes = (size(a.tab)[2], size(b.tab)[2])

    if sizes == (1, 1)
        return CliffordOperator(permute!(a.tab ⊗ b.tab, [1, 3, 2, 4]))
    end

    if sizes == (2, 2)
        return CliffordOperator(permute!(a.tab ⊗ b.tab, [1, 2, 5, 6, 3, 4, 7, 8]))
    end

    if sizes == (1, 2)
        return CliffordOperator(permute!(a.tab ⊗ b.tab, [1, 3, 4, 2, 5, 6]))
    end

    if sizes == (2, 1)
        return CliffordOperator(permute!(a.tab ⊗ b.tab, [1, 2, 5, 3, 4, 6]))
    end

    if sizes == (3, 1)
        return CliffordOperator(permute!(a.tab ⊗ b.tab, [1, 2, 3, 7, 4, 5, 6, 8]))
    end

    throw("unsupported dimensions")

    # this is stupid, but QuantumClifford doesn't yet
    # implement ⊗ for CliffordOperators, but it does
    # for Stabilizers
end

function Base.hash(a::Stabilizer, h::UInt)
    return hash((a.phases, a.xzs), h)
end

function Base.hash(a::CliffordOperator, h::UInt)
    return hash(a.tab, h)
end

function Base.hash(a::CM, h::UInt)
    return hash(a.c)
end



Base.:(*)(a::CM, b::CM) = CM(a.c*b.c, a.m*b.m)

function print(a::CM)
    println(a.c)
    print("\n")
    print(a.m)
end



cxc = CNOT
cxm = Matrix([1. 0 0 0 
              0 1. 0 0 
              0 0 0 1.
              0 0 1. 0])

cx = CM(cxc, cxm)

cxxc = C"X_
         XX
         ZZ
         _Z"

cxxm = Matrix([1. 0 0 0
               0 0 0 1.
               0 0 1. 0
               0 1. 0 0])
cxx = CM(cxxc, cxxm)

czc = C"XZ
        ZX
        Z_
        _Z"

czm = Matrix([1. 0 0 0
              0 1. 0 0
              0 0 1. 0
              0 0 0 -1.])
cz = CM(czc, czm)

function kron(a, b)
   sa = size(a)
   sb = size(b)
   m = zeros(ComplexF64, sa[1]*sb[1], sa[2]*sb[2])
   for i=1:sa[1]
       for j=1:sa[2]
           for k=1:sb[1]
               for l=1:sb[2]
                   m[(i-1)*sb[1]+k, (j-1)*sb[2]+l] = a[i,j] * b[k,l]
               end
           end
       end
   end
   return m
end

h1c = Hadamard ⊗ CliffordId
h2c = CliffordId ⊗ Hadamard

i2m = Matrix([1. 0.
              0. 1.])

hm = 1/sqrt(2) * Matrix([1 1
                         1 -1])
h1m = kron(hm, i2m)
h2m = kron(i2m, hm)

h1 = CM(h1c, h1m)
h2 = CM(h2c, h2m)

id = CM(CliffordId ⊗ CliffordId, kron(i2m, i2m))

swap = cx * cxx * cx

ph = CM(Phase, Matrix([1. 0
                       0  1.0im]))

x = CM(C"+X
         -Z", Matrix([0 1.
                      1. 0]))
z = ph * ph

i1 = z * z
i2 = i1 ⊗ i1

h = CM(Hadamard, hm)


end