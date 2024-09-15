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

    lx = 2
    ly = 2
    l = lx * ly
    config_nfsz = [isodd(n) ? Rep[U₁×U₁]((1 / 2, 1) => 1) : Rep[U₁×U₁]((-1 / 2, 1) => 1) for n = 1:l]#; Rep[U₁×U₁]((0, 2) => 1)]

    # l = 2 * Nk
    # p1 = 1
    # config_nfsz = [Rep[U₁×U₁]((0, 0) => 1) for n = 1:p1-1]
    # config_nfsz = [config_nfsz; Rep[U₁×U₁]((1 / 2, 1) => 1)]
    # config_nfsz = [config_nfsz; [Rep[U₁×U₁]((0, 0) => 1) for n = p1+1:l]]

    bondd = 800
    mpod = 100
    nsweep = 20

    #generate qn according to config_sz for MPSbg or MPSfg
    phySpace = Rep[U₁×U₁]((0, 0) => 1, (1 / 2, 1) => 1, (-1 / 2, 1) => 1, (0, 2) => 1)
    # phySpace = Rep[U₁](0 => 1, 1 => 2, 2 => 1)
    mps = MPSfg(l, config_nfsz, phySpace)
    rnMPS!(mps, l, 1)

    os = []
    for ix = 1:lx-1
        for iy = 1:ly-1
            push!(os, [1 "cdu" (ix - 1) * ly + iy "cu" (ix - 1) * ly + iy + 1])
            push!(os, [1 "cdd" (ix - 1) * ly + iy "cd" (ix - 1) * ly + iy + 1])
            push!(os, [1 "cdu" (ix - 1) * ly + iy + 1 "cu" (ix - 1) * ly + iy])
            push!(os, [1 "cdd" (ix - 1) * ly + iy + 1 "cd" (ix - 1) * ly + iy])

            push!(os, [1 "cdu" (ix - 1) * ly + iy "cu" (ix) * ly + iy])
            push!(os, [1 "cdd" (ix - 1) * ly + iy "cd" (ix) * ly + iy])
            push!(os, [1 "cdu" (ix) * ly + iy "cu" (ix - 1) * ly + iy])
            push!(os, [1 "cdd" (ix) * ly + iy "cd" (ix - 1) * ly + iy])
        end
        push!(os, [1 "cdu" (ix - 1) * ly + ly "cu" (ix - 1) * ly + 1])
        push!(os, [1 "cdd" (ix - 1) * ly + ly "cd" (ix - 1) * ly + 1])
        push!(os, [1 "cdu" (ix - 1) * ly + 1 "cu" (ix - 1) * ly + ly])
        push!(os, [1 "cdd" (ix - 1) * ly + 1 "cd" (ix - 1) * ly + ly])
    end
    push!(os, [1 "cdu" (lx - 1) * ly + ly "cu" (lx - 1) * ly + 1])
    push!(os, [1 "cdd" (lx - 1) * ly + ly "cd" (lx - 1) * ly + 1])
    push!(os, [1 "cdu" (lx - 1) * ly + 1 "cu" (lx - 1) * ly + ly])
    push!(os, [1 "cdd" (lx - 1) * ly + 1 "cd" (lx - 1) * ly + ly])
    for ix = 1:lx-1
        push!(os, [1 "cdu" (ix - 1) * ly + ly "cu" (ix) * ly + ly])
        push!(os, [1 "cdd" (ix - 1) * ly + ly "cd" (ix) * ly + ly])
        push!(os, [1 "cdu" (ix) * ly + ly "cu" (ix - 1) * ly + ly])
        push!(os, [1 "cdd" (ix) * ly + ly "cd" (ix - 1) * ly + ly])
    end
    for iy = 1:ly-1
        push!(os, [1 "cdu" (lx - 1) * ly + iy "cu" (lx - 1) * ly + iy + 1])
        push!(os, [1 "cdd" (lx - 1) * ly + iy "cd" (lx - 1) * ly + iy + 1])
        push!(os, [1 "cdu" (lx - 1) * ly + iy + 1 "cu" (lx - 1) * ly + iy])
        push!(os, [1 "cdd" (lx - 1) * ly + iy + 1 "cd" (lx - 1) * ly + iy])
    end
    # for iy = 1:ly-2
    #     push!(os, [1 "cdu" (lx - 1) * ly + iy "cu" (lx - 1) * ly + iy + 2])
    #     push!(os, [1 "cdd" (lx - 1) * ly + iy "cd" (lx - 1) * ly + iy + 2])
    #     push!(os, [1 "cdu" (lx - 1) * ly + iy + 2 "cu" (lx - 1) * ly + iy])
    #     push!(os, [1 "cdd" (lx - 1) * ly + iy + 2 "cd" (lx - 1) * ly + iy])
    # end
    for il = 1:l
        push!(os, [8 "nud" il])
    end

    # os = []
    # tsla1 = Lamda[1, 2]
    # tsla2 = Lamda[1, 33]
    # for ik1 = 1:Nk
    #     for ik2 = 1:Nk
    #         for ia1 = 1:2
    #             for ia2 = 1:2
    #                 for iq1 = 1:Nk
    #                     for iq2 = 1:Nk
    #                         for ib1 = 1:2
    #                             for ib2 = 1:2
    #                                 if tsla1[(ik1-1)*2+ia1, (ik2-1)*2+ia2] != 0 && tsla2[(iq1-1)*2+ib1, (iq2-1)*2+ib2] != 0
    #                                     push!(os, [tsla1[(ik1-1)*2+ia1, (ik2-1)*2+ia2] * tsla2[(iq1-1)*2+ib1, (iq2-1)*2+ib2] "cdu" (ik1 - 1) * 2 + ia1 "cu" (ik2 - 1) * 2 + ia2 "cdu" (iq1 - 1) * 2 + ib1 "cu" (iq2 - 1) * 2 + ib2])
    #                                 end
    #                             end
    #                         end
    #                     end
    #                 end
    #             end
    #         end
    #     end
    # end

    #define dictionary for this model
    dicop = Dict(
        "cdu" => [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 1 0],
        "cdd" => [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 -1 0 0],
        "cu" => [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0],
        "cd" => [0 0 1 0; 0 0 0 -1; 0 0 0 0; 0 0 0 0],
        "fi" => [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1],
        "ui" => [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],
        "nu" => [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1],
        "nd" => [0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1],
        "nud" => [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1],
    )
    push!(dicop, "sz" => (dicop["nu"] - dicop["nd"]))
    push!(dicop, "sp" => (dicop["cdu"] * dicop["cd"]))
    push!(dicop, "sm" => (dicop["cdd"] * dicop["cu"]))
    push!(dicop, "cducu" => (dicop["cdu"] * dicop["cu"]))
    push!(dicop, "cucdu" => (dicop["cu"] * dicop["cdu"]))
    push!(dicop, "cducd" => (dicop["cdu"] * dicop["cd"]))
    push!(dicop, "cdcdu" => (dicop["cd"] * dicop["cdu"]))
    push!(dicop, "cddcu" => (dicop["cdd"] * dicop["cu"]))
    push!(dicop, "cucdd" => (dicop["cu"] * dicop["cdd"]))
    push!(dicop, "cddcd" => (dicop["cdd"] * dicop["cd"]))
    push!(dicop, "cdcdd" => (dicop["cd"] * dicop["cdd"]))
    dicqn = Dict(
        "cdu" => Irrep[U₁×U₁](1 / 2, 1),
        "cdd" => Irrep[U₁×U₁](-1 / 2, 1),
        "cu" => Irrep[U₁×U₁](-1 / 2, -1),
        "cd" => Irrep[U₁×U₁](1 / 2, -1),
        "fi" => Irrep[U₁×U₁](0, 0),
        "ui" => Irrep[U₁×U₁](0, 0),
        "nu" => Irrep[U₁×U₁](0, 0),
        "nd" => Irrep[U₁×U₁](0, 0),
        "nud" => Irrep[U₁×U₁](0, 0),
        "sz" => Irrep[U₁×U₁](0, 0),
        "sp" => Irrep[U₁×U₁](1, 0),
        "sm" => Irrep[U₁×U₁](-1, 0),
    )
    push!(dicqn, "cducu" => first(dicqn["cdu"] ⊗ dicqn["cu"]))
    push!(dicqn, "cucdu" => first(dicqn["cu"] ⊗ dicqn["cdu"]))
    push!(dicqn, "cducd" => first(dicqn["cdu"] ⊗ dicqn["cd"]))
    push!(dicqn, "cdcdu" => first(dicqn["cd"] ⊗ dicqn["cdu"]))
    push!(dicqn, "cddcu" => first(dicqn["cdd"] ⊗ dicqn["cu"]))
    push!(dicqn, "cucdd" => first(dicqn["cu"] ⊗ dicqn["cdd"]))
    push!(dicqn, "cddcd" => first(dicqn["cdd"] ⊗ dicqn["cd"]))
    push!(dicqn, "cdcdd" => first(dicqn["cd"] ⊗ dicqn["cdd"]))
    dicbf = Dict(
        "cdu" => "f",
        "cdd" => "f",
        "cu" => "f",
        "cd" => "f",
        "fi" => "b",
        "ui" => "b",
        "nu" => "b",
        "nd" => "b",
        "nud" => "b",
        "sz" => "b",
        "sp" => "b",
        "sm" => "b",
        "cducu" => "b",
        "cucdu" => "b",
        "cducd" => "b",
        "cdcdu" => "b",
        "cddcu" => "b",
        "cucdd" => "b",
        "cddcd" => "b",
        "cdcdd" => "b"
    )
    dicphys = Dict(
        Irrep[U₁×U₁](0, 0) => 1,
        Irrep[U₁×U₁](1 / 2, 1) => 2,
        Irrep[U₁×U₁](-1 / 2, 1) => 3,
        Irrep[U₁×U₁](0, 2) => 4
    )

    mpo = autompo(l, os, phySpace)

    # mpo_mpo⁺!(l, mpo)

    lt = Vector{TensorMap}(undef, l + 1)
    rt = Vector{TensorMap}(undef, l + 1)
    ltrt!(lt, rt, mps, mpo, l)

    eei = zeros(ComplexF64, l - 1)
    for is = 1:nsweep
        for il = 1:l-1    #sweep from 1 to l-1
            @tensor mpo2[o1 p11 p21; p12 p22 o2] := mpo[il][o1 p11; p12 o0] * mpo[il+1][o0 p21; p22 o2]
            @tensor mps2[b1 p1 p2 b2] := mps[il][b1; p1 b0] * mps[il+1][b0; p2 b2]
            e1, v1 = eigsolve([lt[il], mpo2, rt[il+2]], mps2, 1, :SR; ishermitian=true, tol=1e-14, krylovdim=3, maxiter=1)
            mps[il], S, V, ϵ = tsvd(v1[1], (1, 2), (3, 4); trunc=truncdim(bondd), alg=TensorKit.SVD())
            mps[il+1] = S * V
            lt[il+1] = cmpsleft(lt[il], mps[il], mpo[il], permute(mps[il]', (2,), (3, 1)))
            println("Sweep $(is) Sites $(il) $(il+1): Energy $(e1[1]), trunc $(ϵ) D $(dim(codomain(mps[il])))")

            if is == nsweep #measure the entanglement entropy for the last sweep
                ps = convert(Array, S)
                for id = axes(ps, 1)
                    p = ps[id, id]^2
                    eei[il] -= p * log(p)
                end
            end
        end

        for il = l:-1:2     #sweep from l to 2
            @tensor mpo2[o1 p11 p21; p12 p22 o2] := mpo[il-1][o1 p11; p12 o0] * mpo[il][o0 p21; p22 o2]
            @tensor mps2[b1 p1 p2 b2] := mps[il-1][b1; p1 b0] * mps[il][b0; p2 b2]
            e1, v1 = eigsolve([lt[il-1], mpo2, rt[il+1]], mps2, 1, :SR; ishermitian=true, tol=1e-14, krylovdim=3, maxiter=1)
            U, S, mps[il], ϵ = tsvd(v1[1], (1, 2), (3, 4); trunc=truncdim(bondd), alg=TensorKit.SVD())
            mps[il-1] = U * S
            rt[il] = cmpsright(rt[il+1], mps[il], mpo[il], permute(mps[il]', (3,), (1, 2)))
            println("Sweep $(is) Sites $(il) $(il-1): Energy $(e1[1]), trunc $(ϵ) D $(dim(codomain(mps[il])))")
        end
    end
end

include("measure.jl")

