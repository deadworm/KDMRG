#Compute s⋅s correlation
using Plots
plotlyjs()

@time begin
    mlt = lt
    mrt = rt

    szc = zeros(l, l)
    rnMPS!(mps, l, 1)
    szSpace = Rep[U₁×U₁]((0, 0) => 1, (1, 0) => 1, (-1, 0) => 1)
    szsz = MPOss(2, 3, l, phySpace, szSpace)

    ltrt!(mlt, mrt, mps, szsz, l)
    for ie = 3:l-1
        mrt[ie] = cmpsright(mrt[ie+1], mps[ie], szsze(phySpace, szSpace), permute(mps[ie]', (3,), (1, 2)))
    end
    for ii = 2:l-1
        ltb = mlt[ii]
        for ie = ii+1:l-1
            @tensor szb = mlt[ie][a, b, c] * mrt[ie][a, b, c]
            szc[ii, ie] = szb
            println("s⋅s between $(ii) and $(ie) is $(szb)")
            mlt[ie+1] = cmpsleft(mlt[ie], mps[ie], szszo(phySpace, szSpace), permute(mps[ie]', (3,), (1, 2)))
        end
        mlt[ii+1] = cmpsleft(ltb, mps[ii], szszo(phySpace, szSpace), permute(mps[ii]', (3,), (1, 2)))
        mlt[ii+2] = cmpsleft(mlt[ii+1], mps[ii+1], szszi(phySpace, szSpace), permute(mps[ii+1]', (3,), (1, 2)))
    end
end

szc = szc + szc'
for il = 2:l-1
    szc[il, il] = 1
end
szc = szc[2:l-1, 2:l-1]

szf = zeros(ComplexF64, l - 2)
for iq = 1:l-2
    for il1 = 1:l-2
        for il2 = 1:l-2
            szf[iq] = szf[iq] + exp(im * 2 * pi * (iq - 1) / (l - 2) * (il1 - il2)) * szc[il1, il2]
        end
    end
end
szf

p1 = plot(szc[l÷2, l÷2:end])
p2 = plot(1:l-2, real(szf))
p3 = plot(real(eei))

plot(p1, p2, p3, layout=(3, 1), legend=false)
