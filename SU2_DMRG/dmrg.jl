#dmrg test for Heisenberg model
using TensorKit
using KrylovKit
using Random
using Combinatorics

@time begin
    include("mpsinitial.jl")
    include("mpoinitial.jl")
    include("contract.jl")
    include("mpsupdate.jl")

    # lx = 2
    # ly = 2
    # l = lx * ly
    # lx=2
    # ly=2
    l = 18

    Rm = Rep[SU₂×U₁]
    phySpace = Rm((0, 0) => 1, (1 / 2, 1) => 1, (0, 2) => 1)
    opzo = Irrep[SU₂×U₁](0, 0)

    #define the state string with length l+1
    # config_nfs = [isodd(n) ? Rm((0, l - n + 1) => 1) : Rm((1 / 2, l - n + 1) => 1) for n = 1:l+1]
    # [Rm((0, 4) => 1); Rm((1 / 2, 3) => 1); Rm((0, 2) => 1); Rm((1 / 2, 1) => 1); Rm((0, 0) => 1)]

    bondd = 800
    mpod = 100
    nsweep = 10

    #generate qn according to config_sz for MPSbg or MPSfg

    p1 = 1:l
    p2 = []
    opfn = Irrep[SU₂×U₁]((l - 2) / 2, l)

    # p1 = 2:l
    # p2 = [1]
    # opfn = Irrep[SU₂×U₁]((l - 1) / 2, l + 1)

    # p1 = 1:l
    # p2 = []
    # opfn = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}]((l ) / 2, l , 0, 0)

    # p1=[]
    # p2=[(ik - 1) * 2 + 1]
    # opfn = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 2, 2*lklist[p2[1], 1], 2*lklist[p2[1], 2])

    # p1 = [2]
    # p2 = []
    # opfn = Irrep[SU₂×U₁](1 / 2, 1)

    statelist = Vector{Irrep[SU₂×U₁]}(undef, l)
    for il = 1:l
        if il in p1
            statelist[il] = Irrep[SU₂×U₁](1 / 2, 1)
        elseif il in p2
            statelist[il] = Irrep[SU₂×U₁](0, 2)
        else
            statelist[il] = opzo
        end
    end
    mps = MPSfg(l, statelist, phySpace, opzo, opfn)
    rdMPS!(l, mps, 100, 100)
    rnMPS!(mps, l, 1)

    #generate hubbard mpo
    # os = []
    # for ix = 1:lx-1
    #     for iy = 1:ly-1
    #         push!(os, [1 "cd" (ix - 1) * ly + iy "c" (ix - 1) * ly + iy + 1])
    #         push!(os, [1 "cd" (ix - 1) * ly + iy + 1 "c" (ix - 1) * ly + iy])

    #         push!(os, [1 "cd" (ix - 1) * ly + iy "c" (ix) * ly + iy])
    #         push!(os, [1 "cd" (ix) * ly + iy "c" (ix - 1) * ly + iy])
    #     end
    #     push!(os, [1 "cd" (ix - 1) * ly + ly "c" (ix - 1) * ly + 1])
    #     push!(os, [1 "cd" (ix - 1) * ly + 1 "c" (ix - 1) * ly + ly])
    # end
    # push!(os, [1 "cd" (lx - 1) * ly + ly "c" (lx - 1) * ly + 1])
    # push!(os, [1 "cd" (lx - 1) * ly + 1 "c" (lx - 1) * ly + ly])
    # for ix = 1:lx-1
    #     push!(os, [1 "cd" (ix - 1) * ly + ly "c" (ix) * ly + ly])
    #     push!(os, [1 "cd" (ix) * ly + ly "c" (ix - 1) * ly + ly])
    # end
    # for iy = 1:ly-1
    #     push!(os, [1 "cd" (lx - 1) * ly + iy "c" (lx - 1) * ly + iy + 1])
    #     push!(os, [1 "cd" (lx - 1) * ly + iy + 1 "c" (lx - 1) * ly + iy])
    # end
    # for iy = 1:ly-2
    #     push!(os, [1 "cd" (lx - 1) * ly + iy "c" (lx - 1) * ly + iy + 2])
    #     push!(os, [1 "cd" (lx - 1) * ly + iy + 2 "c" (lx - 1) * ly + iy])
    # end

    # for il = 1:l
    #     push!(os, [8 "nud" il])
    # end
    #define dictionary for this model
    dicop = op2ten(phySpace)

    mpolist = []
    for inh in eachindex(hlist)
        # inh=1
        os = []
        tsla = hlist[inh][1]
        for ik1 = 1:Nk
            for ik2 = 1:Nk
                if tsla[(ik1-1)*2+1, (ik2-1)*2+1] != 0
                    if ik1 != ik2
                        push!(os, [tsla[(ik1-1)*2+1, (ik2-1)*2+1] "cd" (ik1 - 1) * 2 + 1 "c" (ik2 - 1) * 2 + 1])
                    else
                        push!(os, [tsla[(ik1-1)*2+1, (ik2-1)*2+1] "nall" (ik1 - 1) * 2 + 1])
                    end
                end
                if tsla[(ik1-1)*2+1, (ik2-1)*2+2] != 0
                    push!(os, [tsla[(ik1-1)*2+1, (ik2-1)*2+2] "cd" (ik1 - 1) * 2 + 1 "c" (ik2 - 1) * 2 + 2])
                end
                if tsla[(ik1-1)*2+2, (ik2-1)*2+1] != 0
                    push!(os, [tsla[(ik1-1)*2+2, (ik2-1)*2+1] "cd" (ik1 - 1) * 2 + 2 "c" (ik2 - 1) * 2 + 1])
                end
                if tsla[(ik1-1)*2+2, (ik2-1)*2+2] != 0
                    if ik1 != ik2
                        push!(os, [tsla[(ik1-1)*2+2, (ik2-1)*2+2] "cd" (ik1 - 1) * 2 + 2 "c" (ik2 - 1) * 2 + 2])
                    else
                        push!(os, [tsla[(ik1-1)*2+2, (ik2-1)*2+2] "nall" (ik1 - 1) * 2 + 2])
                    end
                end
            end
        end
        opfn = [Irrep[SU₂×U₁](0, 0)]
        mpo = autompo(l, os, phySpace, Rm)

        cm = cmpo(l, hlist[inh][2], phySpace, Rm)
        mpo = mpo₊mpo(l, cm, mpo)
        mpo = mpo_mpo⁺(l, mpo)
        push!(mpolist, mpo)
    end
    mpo2list = Array{TensorMap}(undef, length(mpolist), l - 1)
    for im in eachindex(mpolist)
        for il = 1:l-1
            @tensor mpo2list[im, il][o1 p11 p21; p12 p22 o2] := mpolist[im][il][o1 p11; p12 o0] * mpolist[im][il+1][o0 p21; p22 o2]
        end
    end

    ltlist = Array{TensorMap}(undef, length(mpolist), l + 1)
    rtlist = Array{TensorMap}(undef, length(mpolist), l + 1)
    ltrt!(ltlist, rtlist, mps, mpolist, l)

    eei = zeros(ComplexF64, l - 1)
    e1 = 0
    for is = 1:nsweep
        for il = 1:l-1    #sweep from 1 to l-1
            @tensor mps2[b1 p1 p2 b2] := mps[il][b1; p1 b0] * mps[il+1][b0; p2 b2]
            e1, v1 = eigsolve([ltlist[:, il], mpo2list[:, il], rtlist[:, il+2]], mps2, 1, :SR; ishermitian=true, tol=1e-14, krylovdim=10, maxiter=1)
            mps[il], S, V, ϵ = tsvd(v1[1], (1, 2), (3, 4); trunc=truncdim(bondd), alg=TensorKit.SVD())
            mps[il+1] = S * V
            for im in eachindex(mpolist)
                ltlist[im, il+1] = cmpsleft(ltlist[im, il], mps[il], mpolist[im][il], permute(mps[il]', (2,), (3, 1)))
            end
            println("Sweep $(is) Sites $(il) $(il+1): Energy $(e1[1]), trunc $(ϵ) D $(dim(domain(mps[il])))")

            if is == nsweep #measure the entanglement entropy for the last sweep
                ps = convert(Array, S)
                for id = axes(ps, 1)
                    p = ps[id, id]^2
                    eei[il] -= p * log(p)
                end
            end
        end

        for il = l:-1:2     #sweep from l to 2
            @tensor mps2[b1 p1 p2 b2] := mps[il-1][b1; p1 b0] * mps[il][b0; p2 b2]
            e1, v1 = eigsolve([ltlist[:, il-1], mpo2list[:, il-1], rtlist[:, il+1]], mps2, 1, :SR; ishermitian=true, tol=1e-14, krylovdim=10, maxiter=1)
            U, S, mps[il], ϵ = tsvd(v1[1], (1, 2), (3, 4); trunc=truncdim(bondd), alg=TensorKit.SVD())
            mps[il-1] = U * S
            for im in eachindex(mpolist)
                rtlist[im, il] = cmpsright(rtlist[im, il+1], mps[il], mpolist[im][il], permute(mps[il]', (3,), (1, 2)))
            end
            println("Sweep $(is) Sites $(il) $(il-1): Energy $(e1[1]), trunc $(ϵ) D $(dim(codomain(mps[il])))")
        end
    end
end

include("measure.jl")