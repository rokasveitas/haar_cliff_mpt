using QuantumClifford
using Serialization
using Random

include("cms.jl")
using .CMs


ddict = deserialize("doubdict") # comes from gates_script.jl


# implementing section 2.2 of https://arxiv.org/pdf/2003.09412.pdf
# for 2-qubit situations
function mat_from_op(s::CliffordOperator)
    nq = 2
    bits = stab_to_gf2(s.tab)

    # First, apply lemma 5 if it needs to be disentangled
    b1 = id
    b2 = id
    b3 = id

    if [sum(bits[i,1:2] .| bits[i,3:4]) for i=1:4] != ones(Int64, 4) # true if s is entangling
        if sum(bits[3,1:2] .| bits[3,3:4]) != 1 # if Z_ needs to be single-qubitified
            alp = bits[3,1:2]
            bet = bits[3,3:4]

            if alp == [false, false] # alpha != 0^n
                if bet == [true, true]
                    b1 = cx
                end
            elseif alp == [false, true]
                if bet[1]
                    b1 = cz
                end
            elseif alp == [true, false]
                if bet[2]
                    b1 = cz
                end
            elseif alp == [true, true]
                if bet[2]
                    b1 = cx * cz
                else
                    b1 = cx
                end
            end
        end

        # now b1 * s * Z_ acts only on k

        u = b1.c * s # b * s is the new operator
        ubits = stab_to_gf2(u.tab)

        k = findall(ubits[3,1:2] .| ubits[3,3:4])[1] # the qubit that u * P"Z_" acts on

        b2 = ubits[4, k] | ubits[4, k+2] ? cx : id # if _Z acts on k qubit

        u = u * b2.c
        ubits = stab_to_gf2(u.tab)

        b3 = ubits[2,k] | ubits[2,k+2] ? cz : id

        # now b1 * s * b2 * b3 is non-entangling

    end

    v = b1.c * s * b2.c * b3.c
    vbits = stab_to_gf2(v.tab)

    # now use lemma 3 to find non-entangling representation

    perm = vbits[1,1] | vbits[1,3] ? id : swap

    v = perm.c * v

    squbs = get(ddict, v, nothing)

    # v = perm * b1 * s * b2 * b3
    # b1 * perm * v * b3 * b2

    return b1.m * perm.m * squbs * b3.m * b2.m

end

clifset = deserialize("cliffgates")

cmset = Set()

for gate in clifset
    push!(cmset, CM(gate, mat_from_op(gate)))
end

serialize("cmset", cmset)

# [3,1:2]=alpha: 00 go to beta
#              : 01 b1 = id, if bits[3,3] b2=cz
#              : 10 b1 = id, if bits[3,4] b2=cz
#              : 11 b1 = cx, if bits[3,4] b2=cz

