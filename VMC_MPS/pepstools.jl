#For preparing a PEPS with physical index fixed by config_sz
function PEPSinit(lx, ly, config_sz, phySpace)
    PEPS = Vector{TensorMap}(undef, lx * ly)
    zero_qn = Rep[U₁]((0) => 1)
    qn_bash = [zero_qn for i = 1:ly]
    qn_u = zero_qn
    for ix = 1:lx
        for iy = 1:ly
            i = (ix - 1) * ly + iy
            if ix != lx
                qn_r = Rep[U₁](blocksectors(fuse(config_sz[i], qn_bash[iy]))[1] => 1)
                if iy == 1
                    qn_d = zero_qn
                else
                    qn_d = Rep[U₁]((0) => 1)
                end
                if iy == ly
                    qn_u = zero_qn
                else
                    qn_u = Rep[U₁]((0) => 1)
                end
                PEPS[i] = TensorMap(ones, ComplexF64, phySpace ⊗ qn_bash[iy] ⊗ qn_d ← qn_u ⊗ qn_r) + TensorMap(ones, ComplexF64, phySpace ⊗ qn_bash[iy] ⊗ qn_d ← qn_u ⊗ qn_r)
                qn_bash[iy] = qn_r
            else
                if iy == 1
                    qn_d = zero_qn
                else
                    qn_d = qn_u
                end
                if iy == ly
                    qn_u = zero_qn
                else
                    qn_u = Rep[U₁](blocksectors(fuse(config_sz[i], qn_bash[iy], qn_d))[1] => 1)
                end
                PEPS[i] = TensorMap(ones, ComplexF64, phySpace ⊗ qn_bash[iy] ⊗ qn_d ← qn_u ⊗ zero_qn) + TensorMap(ones, ComplexF64, phySpace ⊗ qn_bash[iy] ⊗ qn_d ← qn_u ⊗ zero_qn)
            end
            PEPS[i] = PEPS[i] / 2
            # PEPS[i] = PEPS[i] / sqrt(size(PEPS[i].data, 1))  #normalize the tensor
        end
    end
    #   PEPS diagram:
    #    ↓        ↓        ...   ↓
    #  ← o(1,ly)← o(2,ly)← ... ← o(lx,ly)←
    #   ...      ...       ...  ...
    #    ↓        ↓        ...   ↓
    #  ← o(1,2) ← o(2,2) ← ... ← o(lx,2) ←
    #    ↓        ↓        ...   ↓
    #  ← o(1,1) ← o(2,1) ← ... ← o(lx,1) ←
    #    ↓        ↓        ...   ↓
    # PEPS[i] = TensorMap(ones, Float64, phySpace ⊗ qn_l ⊗ qn_d ← qn_u ⊗ qn_r)
    return PEPS
end

#For randomize a PEPS
function rdPEPS!(lx, ly, peps, mdim, dt, nsweep)
    for ins = 1:nsweep
        for ix = 1:lx
            for iy = 1:2:ly-1
                i1 = (ix - 1) * ly + iy
                i2 = (ix - 1) * ly + iy + 1
                tt = TensorMap(rand, ComplexF64, phySpace ⊗ phySpace ← phySpace ⊗ phySpace)
                qb = exp(im * (tt' + tt) * dt)
                # q, r = leftorth(tt)
                # qb = q * isometry(domain(q) ← domain(tt))
                @tensor pepsb[p1 l1 d1 r1; p2 l2 u2 r2] := peps[i1][o1 l1 d1; o0 r1] * peps[i2][o2 l2 o0; u2 r2] * qb[p1 p2; o1 o2]
                U, S, V, ϵ = tsvd(pepsb, (1, 2, 3, 4), (5, 6, 7, 8); trunc=truncdim(mdim), alg=TensorKit.SVD())
                peps[i1] = permute(U * sqrt(S), (1, 2, 3), (5, 4))
                peps[i2] = permute(sqrt(S) * V, (2, 3, 1), (4, 5))
                println("Sweep $(ins) Sites $(i1) $(i2): trunc $(ϵ) D $(dim(domain(U)))")
            end
            for iy = ly-2:-2:2
                i1 = (ix - 1) * ly + iy
                i2 = (ix - 1) * ly + iy + 1
                tt = TensorMap(rand, ComplexF64, phySpace ⊗ phySpace ← phySpace ⊗ phySpace)
                qb = exp(im * (tt' + tt) * dt)
                # q, r = leftorth(tt)
                # qb = q * isometry(domain(q) ← domain(tt))
                @tensor pepsb[p1 l1 d1 r1; p2 l2 u2 r2] := peps[i1][o1 l1 d1; o0 r1] * peps[i2][o2 l2 o0; u2 r2] * qb[p1 p2; o1 o2]
                U, S, V, ϵ = tsvd(pepsb, (1, 2, 3, 4), (5, 6, 7, 8); trunc=truncdim(mdim), alg=TensorKit.SVD())
                peps[i1] = permute(U * sqrt(S), (1, 2, 3), (5, 4))
                peps[i2] = permute(sqrt(S) * V, (2, 3, 1), (4, 5))
                println("Sweep $(ins) Sites $(i1) $(i2): trunc $(ϵ) D $(dim(domain(U)))")
            end
        end
        for iy = 1:ly
            for ix = 1:2:lx-1
                i1 = (ix - 1) * ly + iy
                i2 = (ix) * ly + iy
                tt = TensorMap(rand, ComplexF64, phySpace ⊗ phySpace ← phySpace ⊗ phySpace)
                qb = exp(im * (tt' + tt) * dt)
                # q, r = leftorth(tt)
                # qb = q * isometry(domain(q) ← domain(tt))
                @tensor pepsb[p1 l1 d1 u1; p2 d2 u2 r2] := peps[i1][o1 l1 d1; u1 o0] * peps[i2][o2 o0 d2; u2 r2] * qb[p1 p2; o1 o2]
                U, S, V, ϵ = tsvd(pepsb, (1, 2, 3, 4), (5, 6, 7, 8); trunc=truncdim(mdim), alg=TensorKit.SVD())
                peps[i1] = permute(U * sqrt(S), (1, 2, 3), (4, 5))
                peps[i2] = permute(sqrt(S) * V, (2, 1, 3), (4, 5))
                println("Sweep $(ins) Sites $(i1) $(i2): trunc $(ϵ) D $(dim(domain(U)))")
            end
            for ix = lx-2:-2:2
                i1 = (ix - 1) * ly + iy
                i2 = (ix) * ly + iy
                tt = TensorMap(rand, ComplexF64, phySpace ⊗ phySpace ← phySpace ⊗ phySpace)
                qb = exp(im * (tt' + tt) * dt)
                # q, r = leftorth(tt)
                # qb = q * isometry(domain(q) ← domain(tt))
                @tensor pepsb[p1 l1 d1 u1; p2 d2 u2 r2] := peps[i1][o1 l1 d1; u1 o0] * peps[i2][o2 o0 d2; u2 r2] * qb[p1 p2; o1 o2]
                U, S, V, ϵ = tsvd(pepsb, (1, 2, 3, 4), (5, 6, 7, 8); trunc=truncdim(mdim), alg=TensorKit.SVD())
                peps[i1] = permute(U * sqrt(S), (1, 2, 3), (4, 5))
                peps[i2] = permute(sqrt(S) * V, (2, 1, 3), (4, 5))
                println("Sweep $(ins) Sites $(i1) $(i2): trunc $(ϵ) D $(dim(domain(U)))")
            end
        end
    end
end

#For preparing the contraction of a partial PEPS at ith line
function PEPSlc!(lt, rt, mps_u, mps_d, peps, iy, lx)
    if mps_d == 0
        @tensor ltbash[p1 p2 l1 l2 d2; u1 r1 r2] := mps_u[1][p1 l1 o0; u1 r1] * peps[iy][p2 l2 d2; o0 r2]
        piso = isometry(fuse(codomain(ltbash, 1) ⊗ codomain(ltbash, 2)), codomain(ltbash, 1) ⊗ codomain(ltbash, 2))
        liso = isometry(fuse(codomain(ltbash, 3) ⊗ codomain(ltbash, 4)), codomain(ltbash, 3) ⊗ codomain(ltbash, 4))
        @tensor lt[1][p l d; u r1 r2] := ltbash[p1 p2 l1 l2 d; u r1 r2] * piso[p; p1 p2] * liso[l; l1 l2]
        #   lt[p l d; u r1 r2] diagram:
        #     ↓    ←
        #  ←  o(i) ←
        #   ↙ ↓    
        #  p

        @tensor rtbash[p1 p2 l1 l2 d2; u1 r1 r2] := mps_u[lx][p1 l1 o0; u1 r1] * peps[(lx-1)*lx+iy][p2 l2 d2; o0 r2]
        piso = isometry(fuse(codomain(rtbash, 1) ⊗ codomain(rtbash, 2)), codomain(rtbash, 1) ⊗ codomain(rtbash, 2))
        riso = isometry(domain(rtbash, 2) ⊗ domain(rtbash, 3), fuse(domain(rtbash, 2) ⊗ domain(rtbash, 3)))
        @tensor rt[lx][p l1 l2 d; u r] := rtbash[p1 p2 l1 l2 d; u r1 r2] * piso[p; p1 p2] * riso[r1 r2; r]
        #   rt[p l1 l2 d; u r] diagram:
        #  ←   ↓    
        #  ←   o(i) ←
        #    ↙ ↓    
        #  p

        for ix = 2:lx
            i = (ix - 1) * lx + iy
            ir = (lx - ix) * lx + iy
            @tensor ltbash[p1 p2 l d1 d2; u1 u2 r1 r2] := lt[ix-1][p1 l d1; u1 o1 r2] * mps_u[ix][p2 o1 d2; u2 r1]
            @tensor ltbash[p1 p2 p3 l d1 d2; u1 u2 r1 r2] := ltbash[p1 p2 l d1 o0; u1 u2 r1 o2] * peps[i][p3 o2 d2; o0 r2]
            piso = isometry(fuse(codomain(ltbash, 1) ⊗ codomain(ltbash, 2) ⊗ codomain(ltbash, 3)), codomain(ltbash, 1) ⊗ codomain(ltbash, 2) ⊗ codomain(ltbash, 3))
            diso = isometry(fuse(codomain(ltbash, 5) ⊗ codomain(ltbash, 6)), codomain(ltbash, 5) ⊗ codomain(ltbash, 6))
            uiso = isometry(domain(ltbash, 1) ⊗ domain(ltbash, 2), fuse(domain(ltbash, 1) ⊗ domain(ltbash, 2)))
            @tensor lt[ix][p l d; u r1 r2] := ltbash[p1 p2 p3 l d1 d2; u1 u2 r1 r2] * piso[p; p1 p2 p3] * diso[d; d1 d2] * uiso[u1 u2; u]

            @tensor rtbash[p1 p2 l1 l2 d1 d2; u1 u2 r] := rt[lx-ix+2][p1 o1 l2 d1; u1 r] * mps_u[lx-ix+1][p2 l1 d2; u2 o1]
            @tensor rtbash[p1 p2 p3 l1 l2 d1 d2; u1 u2 r] := rtbash[p1 p2 l1 o2 d1 o0; u1 u2 r] * peps[ir][p3 l2 d2; o0 o2]
            piso = isometry(fuse(codomain(rtbash, 1) ⊗ codomain(rtbash, 2) ⊗ codomain(rtbash, 3)), codomain(rtbash, 1) ⊗ codomain(rtbash, 2) ⊗ codomain(rtbash, 3))
            diso = isometry(fuse(codomain(rtbash, 6) ⊗ codomain(rtbash, 7)), codomain(rtbash, 6) ⊗ codomain(rtbash, 7))
            uiso = isometry(domain(rtbash, 1) ⊗ domain(rtbash, 2), fuse(domain(rtbash, 1) ⊗ domain(rtbash, 2)))
            @tensor rt[lx-ix+1][p l1 l2 d; u r] := rtbash[p1 p2 p3 l1 l2 d1 d2; u1 u2 r] * piso[p; p1 p2 p3] * diso[d; d1 d2] * uiso[u1 u2; u]
        end
    elseif mps_u == 0
        @tensor ltbash[p2 p3 l2 l3 d3; u2 r2 r3] := mps_d[1][p3 l3 d3; o0 r3] * peps[iy][p2 l2 o0; u2 r2]
        piso = isometry(fuse(codomain(ltbash, 1) ⊗ codomain(ltbash, 2)), codomain(ltbash, 1) ⊗ codomain(ltbash, 2))
        liso = isometry(fuse(codomain(ltbash, 3) ⊗ codomain(ltbash, 4)), codomain(ltbash, 3) ⊗ codomain(ltbash, 4))
        @tensor lt[1][p l d; u r2 r3] := ltbash[p2 p3 l2 l3 d; u r2 r3] * piso[p; p2 p3] * liso[l; l2 l3]
        #   lt[p l d; u r2 r3] diagram:
        #     ↓    
        #  ←  o(i) ←
        #   ↙ ↓    ←
        #  p

        @tensor rtbash[p2 p3 l2 l3 d3; u2 r2 r3] := mps_d[lx][p3 l3 d3; o0 r3] * peps[(lx-1)*lx+iy][p2 l2 o0; u2 r2]
        piso = isometry(fuse(codomain(rtbash, 1) ⊗ codomain(rtbash, 2)), codomain(rtbash, 1) ⊗ codomain(rtbash, 2))
        riso = isometry(domain(rtbash, 2) ⊗ domain(rtbash, 3), fuse(domain(rtbash, 2) ⊗ domain(rtbash, 3)))
        @tensor rt[lx][p l2 l3 d; u r] := rtbash[p2 p3 l2 l3 d; u r2 r3] * piso[p; p2 p3] * riso[r2 r3; r]
        #   rt[p l2 l3 d; u r] diagram:
        #      ↓    
        #  ←   o(i) ←
        #  ← ↙ ↓    
        #  p

        for ix = 2:lx
            i = (ix - 1) * lx + iy
            ir = (lx - ix) * lx + iy
            @tensor ltbash[p1 p2 l d1 d2; u1 u2 r2 r3] := lt[ix-1][p1 l d1; u1 r2 o3] * mps_d[ix][p2 o3 d2; u2 r3]
            @tensor ltbash[p1 p2 p3 l d1 d2; u1 u2 r2 r3] := ltbash[p1 p2 l d1 d2; u1 o0 o2 r3] * peps[i][p3 o2 o0; u2 r2]
            piso = isometry(fuse(codomain(ltbash, 1) ⊗ codomain(ltbash, 2) ⊗ codomain(ltbash, 3)), codomain(ltbash, 1) ⊗ codomain(ltbash, 2) ⊗ codomain(ltbash, 3))
            diso = isometry(fuse(codomain(ltbash, 5) ⊗ codomain(ltbash, 6)), codomain(ltbash, 5) ⊗ codomain(ltbash, 6))
            uiso = isometry(domain(ltbash, 1) ⊗ domain(ltbash, 2), fuse(domain(ltbash, 1) ⊗ domain(ltbash, 2)))
            @tensor lt[ix][p l d; u r2 r3] := ltbash[p1 p2 p3 l d1 d2; u1 u2 r2 r3] * piso[p; p1 p2 p3] * diso[d; d1 d2] * uiso[u1 u2; u]

            @tensor rtbash[p1 p2 l2 l3 d1 d2; u1 u2 r] := rt[lx-ix+2][p1 l2 o3 d1; u1 r] * mps_d[lx-ix+1][p2 l3 d2; u2 o3]
            @tensor rtbash[p1 p2 p3 l2 l3 d1 d2; u1 u2 r] := rtbash[p1 p2 o2 l3 d1 d2; u1 o0 r] * peps[ir][p3 l2 o0; u2 o2]
            piso = isometry(fuse(codomain(rtbash, 1) ⊗ codomain(rtbash, 2) ⊗ codomain(rtbash, 3)), codomain(rtbash, 1) ⊗ codomain(rtbash, 2) ⊗ codomain(rtbash, 3))
            diso = isometry(fuse(codomain(rtbash, 6) ⊗ codomain(rtbash, 7)), codomain(rtbash, 6) ⊗ codomain(rtbash, 7))
            uiso = isometry(domain(rtbash, 1) ⊗ domain(rtbash, 2), fuse(domain(rtbash, 1) ⊗ domain(rtbash, 2)))
            @tensor rt[lx-ix+1][p l2 l3 d; u r] := rtbash[p1 p2 p3 l2 l3 d1 d2; u1 u2 r] * piso[p; p1 p2 p3] * diso[d; d1 d2] * uiso[u1 u2; u]
        end
    else
        @tensor ltbash[p1 p2 l1 l2 d2; u1 r1 r2] := mps_u[1][p1 l1 o0; u1 r1] * peps[iy][p2 l2 d2; o0 r2]
        @tensor ltbash[p1 p2 p3 l1 l2 l3 d3; u1 r1 r2 r3] := mps_d[1][p3 l3 d3; o0 r3] * ltbash[p1 p2 l1 l2 o0; u1 r1 r2]
        piso = isometry(fuse(codomain(ltbash, 1) ⊗ codomain(ltbash, 2) ⊗ codomain(ltbash, 3)), codomain(ltbash, 1) ⊗ codomain(ltbash, 2) ⊗ codomain(ltbash, 3))
        liso = isometry(fuse(codomain(ltbash, 4) ⊗ codomain(ltbash, 5) ⊗ codomain(ltbash, 6)), codomain(ltbash, 4) ⊗ codomain(ltbash, 5) ⊗ codomain(ltbash, 6))
        @tensor lt[1][p l d; u r1 r2 r3] := ltbash[p1 p2 p3 l1 l2 l3 d; u r1 r2 r3] * piso[p; p1 p2 p3] * liso[l; l1 l2 l3]
        #   lt[p l d; u r1 r2 r3] diagram:
        #     ↓    ←
        #  ←  o(i) ←
        #   ↙ ↓    ←
        #  p

        @tensor rtbash[p1 p2 l1 l2 d2; u1 r1 r2] := mps_u[lx][p1 l1 o0; u1 r1] * peps[(lx-1)*lx+iy][p2 l2 d2; o0 r2]
        @tensor rtbash[p1 p2 p3 l1 l2 l3 d3; u1 r1 r2 r3] := mps_d[lx][p3 l3 d3; o0 r3] * rtbash[p1 p2 l1 l2 o0; u1 r1 r2]
        piso = isometry(fuse(codomain(rtbash, 1) ⊗ codomain(rtbash, 2) ⊗ codomain(rtbash, 3)), codomain(rtbash, 1) ⊗ codomain(rtbash, 2) ⊗ codomain(rtbash, 3))
        riso = isometry(domain(rtbash, 2) ⊗ domain(rtbash, 3) ⊗ domain(rtbash, 4), fuse(domain(rtbash, 2) ⊗ domain(rtbash, 3) ⊗ domain(rtbash, 4)))
        @tensor rt[lx][p l1 l2 l3 d; u r] := rtbash[p1 p2 p3 l1 l2 l3 d; u r1 r2 r3] * piso[p; p1 p2 p3] * riso[r1 r2 r3; r]
        #   rt[p l1 l2 l3 d; u r] diagram:
        #  ←   ↓    
        #  ←   o(i) ←
        #  ← ↙ ↓    
        #  p

        for ix = 2:lx
            i = (ix - 1) * lx + iy
            ir = (lx - ix) * lx + iy
            @tensor ltbash[p1 p2 l d1 d2; u1 u2 r1 r2 r3] := lt[ix-1][p1 l d1; u1 o1 r2 r3] * mps_u[ix][p2 o1 d2; u2 r1]
            @tensor ltbash[p1 p2 p3 l d1 d2; u1 u2 r1 r2 r3] := ltbash[p1 p2 l d1 o0; u1 u2 r1 o2 r3] * peps[i][p3 o2 d2; o0 r2]
            @tensor ltbash[p1 p2 p3 p4 l d1 d2; u1 u2 r1 r2 r3] := ltbash[p1 p2 p3 l d1 o0; u1 u2 r1 r2 o3] * mps_d[ix][p4 o3 d2; o0 r3]
            piso = isometry(fuse(codomain(ltbash, 1) ⊗ codomain(ltbash, 2) ⊗ codomain(ltbash, 3) ⊗ codomain(ltbash, 4)), codomain(ltbash, 1) ⊗ codomain(ltbash, 2) ⊗ codomain(ltbash, 3) ⊗ codomain(ltbash, 4))
            diso = isometry(fuse(codomain(ltbash, 6) ⊗ codomain(ltbash, 7)), codomain(ltbash, 6) ⊗ codomain(ltbash, 7))
            uiso = isometry(domain(ltbash, 1) ⊗ domain(ltbash, 2), fuse(domain(ltbash, 1) ⊗ domain(ltbash, 2)))
            @tensor lt[ix][p l d; u r1 r2 r3] := ltbash[p1 p2 p3 p4 l d1 d2; u1 u2 r1 r2 r3] * piso[p; p1 p2 p3 p4] * diso[d; d1 d2] * uiso[u1 u2; u]

            @tensor rtbash[p1 p2 l1 l2 l3 d1 d2; u1 u2 r] := rt[lx-ix+2][p1 o1 l2 l3 d1; u1 r] * mps_u[lx-ix+1][p2 l1 d2; u2 o1]
            @tensor rtbash[p1 p2 p3 l1 l2 l3 d1 d2; u1 u2 r] := rtbash[p1 p2 l1 o2 l3 d1 o0; u1 u2 r] * peps[ir][p3 l2 d2; o0 o2]
            @tensor rtbash[p1 p2 p3 p4 l1 l2 l3 d1 d2; u1 u2 r] := rtbash[p1 p2 p3 l1 l2 o3 d1 o0; u1 u2 r] * mps_d[lx-ix+1][p4 l3 d2; o0 o3]
            piso = isometry(fuse(codomain(rtbash, 1) ⊗ codomain(rtbash, 2) ⊗ codomain(rtbash, 3) ⊗ codomain(rtbash, 4)), codomain(rtbash, 1) ⊗ codomain(rtbash, 2) ⊗ codomain(rtbash, 3) ⊗ codomain(rtbash, 4))
            diso = isometry(fuse(codomain(rtbash, 8) ⊗ codomain(rtbash, 9)), codomain(rtbash, 8) ⊗ codomain(rtbash, 9))
            uiso = isometry(domain(rtbash, 1) ⊗ domain(rtbash, 2), fuse(domain(rtbash, 1) ⊗ domain(rtbash, 2)))
            @tensor rt[lx-ix+1][p l1 l2 l3 d; u r] := rtbash[p1 p2 p3 p4 l1 l2 l3 d1 d2; u1 u2 r] * piso[p; p1 p2 p3 p4] * diso[d; d1 d2] * uiso[u1 u2; u]
        end
    end
end

#For contracting a partial PEPS at ith line
function pPEPSlc!(ppeps, lt, rt, mps_u, mps_d, ix, iy, lx)
    i = (ix - 1) * lx + iy
    if mps_d == 0
        if lt == 0
            @tensor pbash[p1 p2 l d n2 n3; u1 u2 r] := mps_u[ix][p1 l n2; u1 o4] * rt[ix+1][p2 o4 n3 d; u2 r]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2)), codomain(pbash, 1) ⊗ codomain(pbash, 2))
            uiso = isometry(domain(pbash, 1) ⊗ domain(pbash, 2), fuse(domain(pbash, 1) ⊗ domain(pbash, 2)))
            @tensor ppeps[i][p l d n2 n3; u r] := pbash[p1 p2 l d n2 n3; u1 u2 r] * piso[p; p1 p2] * uiso[u1 u2; u]
        elseif rt == 0
            @tensor pbash[p1 p2 l d n2; u1 u2 r n1] := lt[ix-1][p1 l d; u1 o1 n1] * mps_u[ix][p2 o1 n2; u2 r]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2)), codomain(pbash, 1) ⊗ codomain(pbash, 2))
            uiso = isometry(domain(pbash, 1) ⊗ domain(pbash, 2), fuse(domain(pbash, 1) ⊗ domain(pbash, 2)))
            @tensor ppeps[i][p l d n2; u r n1] := pbash[p1 p2 l d n2; u1 u2 r n1] * piso[p; p1 p2] * uiso[u1 u2; u]
        else
            @tensor pbash[p1 p2 p3 l d1 d2 n2 n3; u1 u2 u3 r n1] := lt[ix-1][p1 l d1; u1 o1 n1] * mps_u[ix][p2 o1 n2; u2 o4] * rt[ix+1][p3 o4 n3 d2; u3 r]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3)), codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3))
            diso = isometry(fuse(codomain(pbash, 5) ⊗ codomain(pbash, 6)), codomain(pbash, 5) ⊗ codomain(pbash, 6))
            uiso = isometry(domain(pbash, 1) ⊗ domain(pbash, 2) ⊗ domain(pbash, 3), fuse(domain(pbash, 1) ⊗ domain(pbash, 2) ⊗ domain(pbash, 3)))
            @tensor ppeps[i][p l d n2 n3; u r n1] := pbash[p1 p2 p3 l d1 d2 n2 n3; u1 u2 u3 r n1] * piso[p; p1 p2 p3] * diso[d; d1 d2] * uiso[u1 u2 u3; u]
        end
    elseif mps_u == 0
        if lt == 0
            @tensor pbash[p1 p2 l d1 d2 n3; u r n4] := rt[ix+1][p1 n3 o6 d1; u r] * mps_d[ix][p2 l d2; n4 o6]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2)), codomain(pbash, 1) ⊗ codomain(pbash, 2))
            diso = isometry(fuse(codomain(pbash, 4) ⊗ codomain(pbash, 5)), codomain(pbash, 4) ⊗ codomain(pbash, 5))
            @tensor ppeps[i][p l d n3; u r n4] := pbash[p1 p2 l d1 d2 n3; u r n4] * piso[p; p1 p2] * diso[d; d1 d2]
        elseif rt == 0
            @tensor pbash[p1 p2 l d1 d2; u r n1 n4] := lt[ix-1][p1 l d1; u n1 o3] * mps_d[ix][p2 o3 d2; n4 r]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2)), codomain(pbash, 1) ⊗ codomain(pbash, 2))
            diso = isometry(fuse(codomain(pbash, 4) ⊗ codomain(pbash, 5)), codomain(pbash, 4) ⊗ codomain(pbash, 5))
            @tensor ppeps[i][p l d; u r n1 n4] := pbash[p1 p2 l d1 d2; u r n1 n4] * piso[p; p1 p2] * diso[d; d1 d2]
        else
            @tensor pbash[p1 p2 p3 l d1 d2 d3 n3; u1 u2 r n1 n4] := lt[ix-1][p1 l d1; u1 n1 o3] * rt[ix+1][p2 n3 o6 d2; u2 r] * mps_d[ix][p3 o3 d3; n4 o6]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3)), codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3))
            diso = isometry(fuse(codomain(pbash, 5) ⊗ codomain(pbash, 6) ⊗ codomain(pbash, 7)), codomain(pbash, 5) ⊗ codomain(pbash, 6) ⊗ codomain(pbash, 7))
            uiso = isometry(domain(pbash, 1) ⊗ domain(pbash, 2), fuse(domain(pbash, 1) ⊗ domain(pbash, 2)))
            @tensor ppeps[i][p l d n3; u r n1 n4] := pbash[p1 p2 p3 l d1 d2 d3 n3; u1 u2 r n1 n4] * piso[p; p1 p2 p3] * diso[d; d1 d2 d3] * uiso[u1 u2; u]
        end
    else
        if lt == 0
            @tensor pbash[p1 p2 p3 l1 l2 d1 d2 n2 n3; u1 u2 r n4] := mps_u[ix][p1 l1 n2; u1 o4] * rt[ix+1][p2 o4 n3 o6 d1; u2 r] * mps_d[ix][p3 l2 d2; n4 o6]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3)), codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3))
            diso = isometry(fuse(codomain(pbash, 6) ⊗ codomain(pbash, 7)), codomain(pbash, 6) ⊗ codomain(pbash, 7))
            uiso = isometry(domain(pbash, 1) ⊗ domain(pbash, 2), fuse(domain(pbash, 1) ⊗ domain(pbash, 2)))
            liso = isometry(fuse(codomain(pbash, 4) ⊗ codomain(pbash, 5)), codomain(pbash, 4) ⊗ codomain(pbash, 5))
            @tensor ppeps[i][p l d n2 n3; u r n4] := pbash[p1 p2 p3 l1 l2 d1 d2 n2 n3; u1 u2 r n4] * piso[p; p1 p2 p3] * diso[d; d1 d2] * uiso[u1 u2; u] * liso[l; l1 l2]
        elseif rt == 0
            @tensor pbash[p1 p2 p3 l d1 d2 n2; u1 u2 r1 r2 n1 n4] := lt[ix-1][p1 l d1; u1 o1 n1 o3] * mps_u[ix][p2 o1 n2; u2 r1] * mps_d[ix][p3 o3 d2; n4 r2]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3)), codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3))
            diso = isometry(fuse(codomain(pbash, 5) ⊗ codomain(pbash, 6)), codomain(pbash, 5) ⊗ codomain(pbash, 6))
            uiso = isometry(domain(pbash, 1) ⊗ domain(pbash, 2), fuse(domain(pbash, 1) ⊗ domain(pbash, 2)))
            riso = isometry(domain(pbash, 3) ⊗ domain(pbash, 4), fuse(domain(pbash, 3) ⊗ domain(pbash, 4)))
            @tensor ppeps[i][p l d n2; u r n1 n4] := pbash[p1 p2 p3 l d1 d2 n2; u1 u2 r1 r2 n1 n4] * piso[p; p1 p2 p3] * diso[d; d1 d2] * uiso[u1 u2; u] * riso[r1 r2; r]
        else
            @tensor pbash[p1 p2 p3 p4 l d1 d2 d3 n2 n3; u1 u2 u3 r n1 n4] := lt[ix-1][p1 l d1; u1 o1 n1 o3] * mps_u[ix][p2 o1 n2; u2 o4] * rt[ix+1][p3 o4 n3 o6 d2; u3 r] * mps_d[ix][p4 o3 d3; n4 o6]
            piso = isometry(fuse(codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3) ⊗ codomain(pbash, 4)), codomain(pbash, 1) ⊗ codomain(pbash, 2) ⊗ codomain(pbash, 3) ⊗ codomain(pbash, 4))
            diso = isometry(fuse(codomain(pbash, 6) ⊗ codomain(pbash, 7) ⊗ codomain(pbash, 8)), codomain(pbash, 6) ⊗ codomain(pbash, 7) ⊗ codomain(pbash, 8))
            uiso = isometry(domain(pbash, 1) ⊗ domain(pbash, 2) ⊗ domain(pbash, 3), fuse(domain(pbash, 1) ⊗ domain(pbash, 2) ⊗ domain(pbash, 3)))
            @tensor ppeps[i][p l d n2 n3; u r n1 n4] := pbash[p1 p2 p3 p4 l d1 d2 d3 n2 n3; u1 u2 u3 r n1 n4] * piso[p; p1 p2 p3 p4] * diso[d; d1 d2 d3] * uiso[u1 u2 u3; u]
        end
    end
end