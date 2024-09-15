#For preparing a MPS
function MPSfg(l, statelist, phySpace, opzo, opfn)
    MPS = Vector{TensorMap}(undef, l)
    inbond = first(collect(fusiontrees((reverse(statelist)...,), opfn))).innerlines
    Vr = statelist[l]
    MPS[l] = TensorMap(ones, Float64, Rm(Vr => 1) ← phySpace[l] ⊗ Rm(opzo => 1))
    for i = l-1:-1:2
        Vl = inbond[l-i]
        MPS[i] = TensorMap(ones, Float64, Rm(Vl => 1) ← phySpace[i] ⊗ Rm(Vr => 1))
        Vr = Vl
    end
    MPS[1] = TensorMap(ones, Float64, Rm(opfn => 1) ← phySpace[1] ⊗ Rm(Vr => 1))
    #   MPS diagram:
    #  ← o1 ← o2 ← ... ← ol ←
    #    ↑    ↑    ...   ↑
    return MPS
end

#For randomize a MPS
function rdMPS!(l, mps, mdim, nsweep)
    for ins = 1:nsweep
        for sl = 1:l-1
            tt = TensorMap(rand, ComplexF64, phySpace[sl]' ⊗ phySpace[sl+1]' ← phySpace[sl]' ⊗ phySpace[sl+1]')
            q, r = leftorth(tt)
            qb = isometry(domain(q) ← domain(tt))
            @tensor mps2[p1 p2; b1 b2] := mps[sl][b1 p1; b0] * mps[sl+1][b0 p2; b2]
            mps2n = permute(q * qb * mps2, (3,1),(2,4))
            U, S, V, ϵ = tsvd(mps2n, (1, 2), (3, 4); trunc=truncdim(mdim), alg=TensorKit.SVD())
            mps[sl] = permute(U, (1,), (2, 3))
            mps[sl+1] = S * V
            # println("Sweep $(ins) Sites $(sl) $(sl+1): trunc $(ϵ) D $(dim(domain(mps[sl])))")
        end
        for sl = l-1:-1:1
            tt = TensorMap(rand, ComplexF64, phySpace[sl]' ⊗ phySpace[sl+1]' ← phySpace[sl]' ⊗ phySpace[sl+1]')
            q, r = leftorth(tt)
            qb = isometry(domain(q) ← domain(tt))
            @tensor mps2[p1 p2; b1 b2] := mps[sl][b1 p1; b0] * mps[sl+1][b0 p2; b2]
            mps2n = q * qb * mps2
            U, S, V, ϵ = tsvd(mps2n, (3, 1), (2, 4); trunc=truncdim(mdim), alg=TensorKit.SVD())
            mps[sl] = permute(U * S, (1,), (2, 3))
            mps[sl+1] = V
            # println("Sweep $(ins) Sites $(sl) $(sl+1): trunc $(ϵ) D $(dim(domain(mps[sl])))")
        end
    end
end

#For Leftnorm MPS up to i
#   MPS diagram:
#  ← o1 ← o2 ← ... ← oi ← ...
#    ↓    ↓    ...   ↓    ...
function lnMPS!(mps, L, i)
    q, r = leftorth(mps[1], (1, 2), (3,))
    mps[1] = q
    for il = 2:i
        q, r = leftorth(r * permute(mps[il], (1,), (2, 3)), (1, 2), (3,))
        mps[il] = q
    end
end

#For Rightnorm MPS up to i
#   MPS diagram:
#  ...  ← oi ← oi ← ... ← ol ←
#  ...    ↑    ↑    ...   ↑
function rnMPS!(mps, L, i)
    l, q = rightorth(mps[L], (1,), (2, 3))
    mps[L] = q
    for il = (L-1):-1:i
        l, q = rightorth(permute(mps[il], (1, 2), (3,)) * l, (1,), (2, 3))
        mps[il] = q
    end
end