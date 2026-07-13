#For SRC compression of boundary MPS
function SRCcontract(pepsi1, pepsi2, li, bond_cut)
    # bond_cut = 10
    # li = lx
    # pepsi1 = Vector{TensorMap}(undef, li)
    # pepsi2 = Vector{TensorMap}(undef, li)
    # for ix = 1:li
    #     pepsi1[ix] = peps[(ix-1)*ly+1]
    #     pepsi2[ix] = peps[(ix-1)*ly+2]
    # end
    svddim = bond_cut
    bond_cut = svddim * 2
    Ωlist = Vector{TensorMap}(undef, li - 1)
    for il = 1:li-1
        Ωlist[il] = TensorMap(randn, ComplexF64, domain(pepsi2[il], 1) ← Rep[U₁](blocksectors(domain(pepsi2[il], 1))[1] => bond_cut))
    end
    Blist = Vector{TensorMap}(undef, li)
    for il = 1:li
        piso = isometry(fuse(codomain(pepsi2[il], 1) ⊗ codomain(pepsi1[il], 1)) ← codomain(pepsi2[il], 1) ⊗ codomain(pepsi1[il], 1))
        @tensor Blist[il][p l2 l1 d1; u2 r2 r1] := pepsi2[il][p2 l2 o1; u2 r2] * pepsi1[il][p1 l1 d1; o1 r1] * piso[p; p2 p1]
    end
    Clist = Vector{TensorMap}(undef, li - 1)
    liso = isometry(fuse(codomain(Blist[1], 2) ⊗ codomain(Blist[1], 3)) ← codomain(Blist[1], 2) ⊗ codomain(Blist[1], 3))
    @tensor Clist[1][p l d; u r2 r1] := Ωlist[1][o1; u] * Blist[1][p l2 l1 d; o1 r2 r1] * liso[l; l2 l1]
    for il = 2:li-1
        data = zeros(ComplexF64, bond_cut, bond_cut, bond_cut)
        for i = 1:bond_cut
            data[i, i, i] = 1.0
        end
        dbash = TensorMap(data, domain(Clist[il-1], 1) ⊗ domain(Ωlist[il], 1) ← Rep[U₁](blocksectors(fuse(domain(Clist[il-1], 1), domain(Ωlist[il], 1)))[1] => bond_cut))
        @tensor Cbash[p l d2 d1; u r2 r1] := Clist[il-1][p l d1; o1 r2 r1] * Ωlist[il][d2; o2] * dbash[o1 o2; u]
        piso = isometry(fuse(codomain(Blist[il], 1) ⊗ codomain(Cbash, 1)) ← codomain(Blist[il], 1) ⊗ codomain(Cbash, 1))
        diso = isometry(fuse(codomain(Blist[il], 4) ⊗ codomain(Cbash, 4)) ← codomain(Blist[il], 4) ⊗ codomain(Cbash, 4))
        @tensor Clist[il][p l d; u r2 r1] := Cbash[p1 l o3 d1; u o2 o1] * Blist[il][p2 o2 o1 d2; o3 r2 r1] * piso[p; p2 p1] * diso[d; d2 d1]
    end
    #   Clist[p l d; u r2 r1] diagram:
    #     ↓
    #  ←  o(i) ←
    #   ↙ ↓    ←
    #  p
    #   Blist[p l2 l1 d; u r2 r1] diagram:
    #      ↓
    #  ←   o(i) ←
    #  ← ↙ ↓    ←
    #  p

    #QB decomposition and truncation
    mps = Vector{TensorMap}(undef, li)
    @tensor Ynow[p1 l d1 u1; p2 d2 u2 r2 r1] := Clist[li-1][p1 l d1; u1 o2 o1] * Blist[li][p2 o2 o1 d2; u2 r2 r1]
    l, q = rightorth(Ynow, (1, 2, 3, 4), (5, 6, 7, 8, 9))
    riso = isometry(domain(q, 4) ⊗ domain(q, 5) ← fuse(domain(q, 4) ⊗ domain(q, 5)))
    @tensor mps[li][p l d; u r] := q[l; p d u r2 r1] * riso[r2 r1; r]
    qp = permute(q, (2, 1, 3), (4, 5, 6))'
    @tensor Rnow[l2 l1; r] := qp[o0 r2 r1; p r d] * Blist[li][p l2 l1 d; o0 r2 r1]
    for il = (li-1):-1:2
        @tensor Rnow[p l2 l1 d; u1 r] := Rnow[o2 o1; r] * Blist[il][p l2 l1 d; u1 o2 o1]
        @tensor Ynow[p1 l d1 u; p2 d2 u1 r] := Rnow[p2 o2 o1 d2; u1 r] * Clist[il-1][p1 l d1; u o2 o1]
        l, q = rightorth(Ynow, (1, 2, 3, 4), (5, 6, 7, 8))
        mps[il] = permute(q, (2, 1, 3), (4, 5))
        qp = mps[il]'
        @tensor Rnow[l2 l1; r] := qp[o1 o2; p r d] * Rnow[p l2 l1 d; o1 o2]
    end
    liso = isometry(fuse(codomain(Blist[1], 2) ⊗ codomain(Blist[1], 3)) ← codomain(Blist[1], 2) ⊗ codomain(Blist[1], 3))
    @tensor mps[1][p l d; u r] := Rnow[o2 o1; r] * Blist[1][p l2 l1 d; u o2 o1] * liso[l; l2 l1]

    #SVD compression
    for il = 1:li-1
        @tensor cmps[p1 l1 d1 u1; p2 d2 u2 r2] := mps[il][p1 l1 d1; u1 o0] * mps[il+1][p2 o0 d2; u2 r2]
        mps[il], S, V, ϵ = tsvd(cmps, (1, 2, 3, 4), (5, 6, 7, 8); trunc=truncdim(svddim), alg=TensorKit.SVD())
        println("Truncate bond between $(il) and $(il+1): Error=$(ϵ) D=$(dim(domain(mps[il])))")
        mps[il] = permute(mps[il], (1, 2, 3), (4, 5))
        mps[il+1] = permute(S * V, (2, 1, 3), (4, 5))
    end
    return mps
end

