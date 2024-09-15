#For preparing a MPS
function MPSfg(l, config_nfsz, phySpace)
    MPS = Vector{TensorMap}(undef, l)
    flow_qn = Rep[U₁×U₁]((0, 0) => 1)
    for i = l:-1:1
        Vr = flow_qn
        flow_qn = fuse(flow_qn, config_nfsz[i])
        Vl = flow_qn
        MPS[i] = TensorMap(ones, Float64, Vl ← phySpace ⊗ Vr)
    end
    #   MPS diagram:
    #  ← o1 ← o2 ← ... ← ol ←
    #    ↑    ↑    ...   ↑
    return MPS
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