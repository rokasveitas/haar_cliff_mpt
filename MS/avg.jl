using Serialization
using LinearAlgebra

s = zeros(28, 6, 12)
denoms = zeros(28, 6, 12)  # this tracks how many actually were computed so
                           # that we can divide by the right number at the end

for i=1:960
    a = try
            deserialize(string("msdata/ms_", i))
        catch y
            continue
    end
    
    for pind=1:28
        for Nind=1:6
            for θind=1:12
                if a[pind, Nind, θind] != 0
                    s[pind, Nind, θind] += a[pind, Nind, θind]
                    denoms[pind, Nind, θind] += 1
                end
            end
        end
    end
end

denoms = map(x -> x==0 ? 1 : x, denoms)
s ./= denoms

# println("s is computed: ")
# println(s)

# Now we want to find out what the critical point is for each value of θ.  We will minimize
# the mean-squared error
#       R(p) = < ( S(N, p, θ) - f_p(N) )^2 >_N
#     f_p(N) = α_p pc ln(N) + β_p.

pcs = zeros(12)
logNs = log.(8:4:28)

pvec = 0.03:0.01:0.30

Rs = zeros(28, 12)

for θi = 1:12
    bestR = Inf
    bestpc = 0

    data = s[:,:,θi]

    αs = zeros(28)
    βs = zeros(28)
    # Rs = zeros(28)

    for pind = 1:28
        data = s[pind,:,θi]
        
        αs[pind] = (6 * dot(logNs, data) - sum(data) * sum(logNs)) / (6 * dot(logNs, logNs) - sum(logNs)^2)
        βs[pind] = (sum(data) - αs[pind] * sum(data)) / 6

        for Ni=1:6
            Rs[pind, θi] += (data[Ni] - αs[pind] * logNs[Ni] + βs[pind])^2
        end
    end
    # println("αs: ")
    # println(αs)
    # println("βs: ")
    # println(βs)
    # println("Rs: ")
    # println(Rs)
end

serialize("Rs", Rs)

serialize("denoms", denoms)



