using Serialization

s = zeros(28, 6, 12)
denoms = zeros(28, 6, 12)

for i=1:960
    a = try
            deserialize(string("ms_", i))
        catch y
            continue
    end
    
    θarg = (i - 1) % 4
    
    for pind=1:28
        for Nind=1:6
            for θind=1:3
                if a[pind, Nind, θind] != 0
                    s[pind, Nind, θind + 3*θarg] += a[pind, Nind, θind]
                    denoms[pind, Nind, θind + 3*θarg] += 1
                end
            end
        end
    end
end

denoms = map(x -> x==0 ? 1 : x, denoms)
s ./= denoms

serialize("S_ave", s)
