#dmrg test for Heisenberg model
using TensorKit
using KrylovKit
using Random
using Combinatorics
using Plots

@time begin
    include("mpsinitial.jl")
    include("mpoinitial.jl")
    include("contract.jl")
    include("mpsupdate.jl")
    include("readdata.jl")

    l = 2 * Nk
    Rm = Rep[SU₂×U₁×ℤ{nx}×ℤ{ny}]
    #define the state string with length l+1
    # p1 = 1
    # config_nfs = [Rm((0, 2) => 1) for n = 1:p1]
    # config_nfs = [config_nfs; [Rm((0, 0) => 1) for n = p1+1:l+1]]

    phySpace = Vector{Any}(undef, l)
    lklist = zeros(Int64, l, 2)
    for il = l:-1:1
        k1b = (il - 1) ÷ 2 ÷ ny
        k2b = mod((il - 1) ÷ 2, ny)
        phySpace[il] = Rm((0, 0, 0, 0) => 1, (1 / 2, 1, k1b, k2b) => 1, (0, 2, k1b * 2, k2b * 2) => 1)
        lklist[il, :] = [k1b, k2b]
    end

    opzo = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 0, 0, 0)

    bondd = 80
    mpod = 10
    nsweep = 10

    #define dictionary for this model
    dicop = op2ten(phySpace)
    #generate mpo
    mpolist = []
    for inh in eachindex(hlist)
        # inh = 1
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
        opfn = [Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 0, hlist[inh][3], hlist[inh][4])]
        st = setsiteterms(l, os, opzo, opfn)

        mpo1 = autompo(l, st, phySpace, Rm, opzo, opfn)
        if opzo == opfn[1]
            cm = cmpo(l, hlist[inh][2], phySpace, Rm)
            mpo1 = mpo₊mpo(l, cm, mpo1)
        end

        mpo1 = mpo_mpo⁺(l, mpo1)
        push!(mpolist, mpo1)
    end

    #test mpo bond contraction and dimension
    for im in eachindex(mpolist)
        mpo = mpolist[im]
        for il = 1:l-1
            vs = (sectors(domain(mpo[il], 2))...,)
            alvs = 0
            for ivs in vs
                alvs += blockdim(domain(mpo[il], 2), ivs)
            end
            if domain(mpo[il], 2) == codomain(mpo[il+1], 1)
            else
                println("dimension wrong!", alvs)
            end
        end
    end

    mpo2list = Array{TensorMap}(undef, length(mpolist), l - 1)
    for im in eachindex(mpolist)
        for il = 1:l-1
            @tensor mpo2list[im, il][o1 p11 p21; p12 p22 o2] := mpolist[im][il][o1 p11; p12 o0] * mpolist[im][il+1][o0 p21; p22 o2]
        end
    end

    #set state
    Ek = zeros(2, Nk)
    for ik = 1:Nk
        # ik=1
        # if ik == 1
        #     p1 = 1:l
        #     p2 = []
        #     opfn = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}]((l - 2) / 2, l, 0, 0)
        # else
        #     p1 = [2:((ik-1)*2); ((ik-1)*2+2):l]
        #     p2 = [(ik - 1) * 2 + 1]
        #     opfn = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}]((l - 2) / 2, l, lklist[p2[1], 1], lklist[p2[1], 2])
        # end

        p1 = [1:((ik-1)*2); ((ik-1)*2+2):l]
        p2 = [(ik - 1) * 2 + 1]
        opfn = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}]((l - 1) / 2, l + 1, lklist[p2[1], 1], lklist[p2[1], 2])

        # p1 = 1:l
        # p2 = []
        # opfn = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}]((l ) / 2, l , 0, 0)

        # p1=[]
        # p2=[(ik - 1) * 2 + 1]
        # opfn = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 2, 2*lklist[p2[1], 1], 2*lklist[p2[1], 2])

        # p1 = [(ik - 1) * 2 + 2]
        # p2 = []
        # opfn = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](1 / 2, 1, lklist[p1[1], 1], lklist[p1[1], 2])

        statelist = Vector{Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}]}(undef, l)
        for il = 1:l
            if il in p1
                statelist[il] = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](1 / 2, 1, lklist[il, 1], lklist[il, 2])
            elseif il in p2
                statelist[il] = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 2, lklist[il, 1] * 2, lklist[il, 2] * 2)
            else
                statelist[il] = opzo
            end
        end

        #generate qn according to config_sz for MPSbg or MPSfg
        mps = MPSfg(l, statelist, phySpace, opzo, opfn)
        rdMPS!(l, mps, 100, 100)
        rnMPS!(mps, l, 1)

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
        Ek[1, ik] = real(e1[1])
    end
    @show Ek
    plot(Ek[1, 1:end], legend=false)
end