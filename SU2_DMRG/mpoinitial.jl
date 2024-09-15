#For two-site measurment s⋅s
# using StridedViews
# using Strided
function MPOss(inl, enl, l, phySpace, bondSpace)
    u0Space = Rep[U₁×U₁]((0, 0) => 1)
    dp = dim(phySpace)
    db = dim(bondSpace)

    szsz = Vector{TensorMap}(undef, l)
    data1 = zeros(1, dp, dp, db)
    data1[1, :, :, 1] = dicop["ui"]
    szsz[1] = TensorMap(data1, u0Space ⊗ phySpace, phySpace ⊗ bondSpace)
    datae = zeros(db, dp, dp, 1)
    datae[1, :, :, 1] = dicop["ui"]
    szsz[l] = TensorMap(datae, bondSpace ⊗ phySpace, phySpace ⊗ u0Space)
    for il = 2:l-1
        szsz[il] = szszo(phySpace, bondSpace)
    end

    szsz[inl] = szszi(phySpace, bondSpace)
    szsz[enl] = szsze(phySpace, bondSpace)
    return szsz
end

function szszi(phySpace, bondSpace)
    dp = dim(phySpace)
    db = dim(bondSpace)

    data1 = zeros(db, dp, dp, db)
    data1[1, :, :, 1] = dicop["sz"]
    data1[1, :, :, 2] = 0.5 * dicop["sp"]
    data1[1, :, :, 3] = 0.5 * dicop["sm"]
    szsz = TensorMap(data1, bondSpace ⊗ phySpace, phySpace ⊗ bondSpace)

    return szsz
end

function szsze(phySpace, bondSpace)
    dp = dim(phySpace)
    db = dim(bondSpace)

    datae = zeros(db, dp, dp, db)
    datae[1, :, :, 1] = dicop["sz"]
    datae[2, :, :, 1] = dicop["sm"]
    datae[3, :, :, 1] = dicop["sp"]
    szsz = TensorMap(datae, bondSpace ⊗ phySpace, phySpace ⊗ bondSpace)

    return szsz
end

function szszo(phySpace, bondSpace)
    dp = dim(phySpace)
    db = dim(bondSpace)

    datai = zeros(db, dp, dp, db)
    datai[1, :, :, 1] = dicop["ui"]
    datai[2, :, :, 2] = dicop["ui"]
    datai[3, :, :, 3] = dicop["ui"]
    szsz = TensorMap(datai, bondSpace ⊗ phySpace, phySpace ⊗ bondSpace)

    return szsz
end

function op2ten(phySpace)
    pspace = phySpace
    dicop = Dict()

    zospace = Irrep[SU₂×U₁](0, 0)
    isspace = Irrep[SU₂×U₁](1 / 2, 1)
    sspace = Irrep[SU₂×U₁](1 / 2, -1)
    #"cdl"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(isspace => 1))
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)] .= -sqrt(2)
    blocks(top)[SU2Irrep(1 / 2)⊠U1Irrep(1)] .= 1
    # at=zeros(4,4,2)
    # at[:,:,1]=[0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 1 0]
    # at[:,:,2]=[0 0 0 0; 0 0 0 0; 1 0 0 0; 0 -1 0 0]
    # top = TensorMap(at, pspace ← pspace ⊗ Rm(isspace=>1))
    push!(dicop, "cdl" => [top, "f", isspace, zospace])
    #"cl"
    top = permute(-top', ((3, 1), (4, 2)))
    push!(dicop, "cl" => [top, "f", sspace, zospace])
    #"cdfr"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(isspace => 1))
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)] .= sqrt(2)
    blocks(top)[SU2Irrep(1/2)⊠U1Irrep(1)] .= 1
    top=permute(top, ((4, 2), (3, 1)))
    push!(dicop, "cdfr" => [top, "f", isspace, sspace])
    #"cfr"
    top = permute(top', ((3, 1), (4, 2)))
    push!(dicop, "cfr" => [top, "f", sspace, isspace])

    isspace = Irrep[SU₂×U₁](0, 0)
    #"ui"
    top = TensorMap(ones, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(isspace => 1))
    push!(dicop, "ui" => [top, "b", isspace, zospace])
    #"fi"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(isspace => 1))
    blocks(top)[SU2Irrep(0)⊠U1Irrep(0)] .= 1
    blocks(top)[SU2Irrep(1 / 2)⊠U1Irrep(1)] .= -1
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)] .= 1
    push!(dicop, "fi" => [top, "b", isspace, zospace])
    #"nall"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(isspace => 1))
    blocks(top)[SU2Irrep(1 / 2)⊠U1Irrep(1)] .= 1
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)] .= 2
    push!(dicop, "nall" => [top, "b", isspace, zospace])
    #"nud"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(isspace => 1))
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)] .= 1
    push!(dicop, "nud" => [top, "b", isspace, zospace])

    #"cd"
    push!(dicop, "cd" => ["", "f", "", ""])
    #"c"
    push!(dicop, "c" => ["", "f", "", ""])
    #"cdf"
    push!(dicop, "cdf" => ["", "f", "", ""])
    #"cf"
    push!(dicop, "cf" => ["", "f", "", ""])
    return dicop
end

function autompo(l, mpochain, phySpace, Rm)
    #define the storage of operator strings with the form "coefficient, term, position, term, position"
    os = mpochain

    #define the storage for each site with the form "quantumnumberin, qunatumnumberout, coefficient, operators
    #outputposition, inputposition, J_W_matrix"
    siteterms = Vector(undef, l)
    #define the amount of nodes in siteterms
    opzo = dicop["ui"][3]
    dimsiteterms = Vector(undef, l + 1)
    siteterms[1] = [opzo opzo 1 "ui" 1 1]
    dimsiteterms[1] = 1
    for il = 2:l-1
        siteterms[il] = [opzo opzo 1 "ui" 1 1]
        siteterms[il] = [siteterms[il]; opzo opzo 1 "ui" 2 2]
        dimsiteterms[il] = 2
    end
    siteterms[l] = [opzo opzo 1 "ui" 2 1]
    dimsiteterms[l] = 2
    dimsiteterms[l+1] = 1
    # o -I-> o -I-> o -I-> o
    #   \                    \
    #    0                     0
    #     -> o -I-> o -I-> o -I-> o 
    #initialize the siteterms

    for (inum, it) in enumerate(os)
    # it=os[18]
        coe = it[1]
        nop = (length(it) - 1) ÷ 2
        opbash = []
        opindex = []
        for inop = 1:nop
            push!(opbash, it[inop*2])
            push!(opindex, it[inop*2+1])
        end
        isort = sortperm(opindex)
        o_parity = levicivita(isort)
        opindex = opindex[isort]
        opbash = opbash[isort]
        iopt = 0
        iopf = 0
        for iop = opindex[1]:opindex[end]
            if !(iop in opindex)
                if (-1)^iopf == 1
                    push!(opbash, "ui")
                    push!(opindex, iop)
                else
                    push!(opbash, "fi")
                    push!(opindex, iop)
                end
            else
                iopt += 1
                if dicop[opbash[iopt]][2] == "f"
                    iopf += 1
                end
                if (-1)^iopf == 1 && dicop[opbash[iopt]][2] == "f"
                    opbash[iopt] = opbash[iopt] * "f"
                end
            end
        end
        isort = sortperm(opindex)
        opindex = opindex[isort]
        opbash = opbash[isort]
        #finish operator string and use operator string to construct siteterms

        opac = dicop["ui"][3]
        for iin = eachindex(opindex)
            iop = opindex[iin]
            if iop == opindex[1] #use the first row as starting and the second row as the ending
                if iop == opindex[end]
                    siteterms[iop] = [siteterms[iop]; [opzo opzo coe opbash[iin] 1 min(2, dimsiteterms[iop+1])]]
                else
                    dimsiteterms[iop+1] += 1
                    if dicop[opbash[iin]][2] == "b"
                        opnow = opbash[iin]
                        siteterms[iop] = [siteterms[iop]; [opzo dicop[opnow][3] coe opnow 1 dimsiteterms[iop+1]]]
                    end
                    if dicop[opbash[iin]][2] == "f"
                        opnow = opbash[iin] * "l"
                        siteterms[iop] = [siteterms[iop]; [opzo dicop[opnow][3] o_parity * coe opnow 1 dimsiteterms[iop+1]]]
                    end
                    opac = first(opac ⊗ dicop[opnow][3])
                end
            elseif iop != opindex[end]
                dimsiteterms[iop+1] += 1
                opnow = opbash[iin]
                siteterms[iop] = [siteterms[iop]; [opac first(opac ⊗ dicop[opnow][3]) 1 opnow dimsiteterms[iop] dimsiteterms[iop+1]]]
                opac = first(opac ⊗ dicop[opnow][3])
            else
                siteterms[iop] = [siteterms[iop]; [opac opzo 1 opbash[iin] * "r" dimsiteterms[iop] min(2, dimsiteterms[iop+1])]]
            end
        end
    end
    #finish siteterms

    #use siteterms to construct mpo
    MPO = Vector{TensorMap}(undef, l)

    for il = 1:l
        # il=2
        ml = maximum(siteterms[il][:, 5])
        mr = maximum(siteterms[il][:, 6])
        mllist = []
        mrlist = []
        for iml = 1:ml
            sortml = siteterms[il][:, 5] .== iml
            if sum(sortml) == 0
                push!(mllist, dicop["ui"][3])
            else
                push!(mllist, siteterms[il][sortml, 1][1])
            end
        end
        for imr = 1:mr
            sortmr = siteterms[il][:, 6] .== imr
            if sum(sortmr) == 0
                push!(mrlist, dicop["ui"][3])
            else
                push!(mrlist, siteterms[il][sortmr, 2][1])
            end
        end
        #push 0 for il=1 to match the dimension
        if il == 1 && mr < 2
            push!(mrlist, 0)
            mr += 1
        end

        #sort to keep the dimension correctly
        bld = Dict()
        brd = Dict()
        dllist = Dict()
        drlist = Dict()
        for iml in unique(mllist)
            push!(bld, iml => count(x -> x == iml, mllist))
            push!(dllist, iml => findall(x -> x == iml, mllist))
        end
        for imr in unique(mrlist)
            push!(brd, imr => count(x -> x == imr, mrlist))
            push!(drlist, imr => findall(x -> x == imr, mrlist))
        end

        blspace = Rm(bld)
        brspace = Rm(brd)

        MPO[il] = TensorMap(zeros, ComplexF64, blspace ⊗ phySpace, phySpace ⊗ brspace)
        for (f1, f2) in fusiontrees(MPO[il])
            dim1 = bld[f1.uncoupled[1]]
            dim4 = brd[f2.uncoupled[2]]
            datab = zeros(Complex{Float64}, dim1, 1, 1, dim4)
            for (i1, id1) in enumerate(dllist[f1.uncoupled[1]])
                for (i4, id4) in enumerate(drlist[f2.uncoupled[2]])
                    if Bool((siteterms[il][:, 5] .== id1)' * (siteterms[il][:, 6] .== id4))
                        po = siteterms[il][siteterms[il][:, 5].==id1.&&siteterms[il][:, 6].==id4, :]
                        phyr = first(dicop[po[4]][4] ⊗ f1.uncoupled[2])
                        if po[4] == "fi" && (f1.coupled == Irrep[SU₂×U₁](0, 2) || f1.coupled == Irrep[SU₂×U₁](0, 0)) #special case for ccd and cdc (one can check "fi" in this space)
                            poc = -po[3]
                        else
                            poc = po[3]
                        end
                        phym = first(block(dicop[po[4]][1], phyr))
                        datab[i1, 1, 1, i4] = poc * phym
                    else
                        datab[i1, 1, 1, i4] = 0
                    end
                end
            end
            MPO[il][f1, f2] = datab
        end
    end
    #finish MPO and use svd to compress
    for il = 1:l-1
        MPO[il], S, V, ϵ = tsvd(MPO[il], (1, 2, 3), (4,); trunc=truncerr(1e-10), alg=TensorKit.SVD())
        # println("MPO $(il) $(il+1): trunc $(ϵ) D $(dim(domain(MPO[il])))")
        MPO[il] = permute(MPO[il], (1, 2), (3, 4))
        mbash = S * V
        @tensor MPO[il+1][o1 p11; p21 o2] := mbash[o1 o0] * MPO[il+1][o0 p11 p21 o2]
    end
    for il = l:-1:2
        U, S, MPO[il], ϵ = tsvd(MPO[il], (1,), (2, 3, 4); trunc=truncerr(1e-10), alg=TensorKit.SVD())
        # println("MPO $(il) $(il-1): trunc $(ϵ) D $(dim(codomain(MPO[il])))")
        MPO[il] = permute(MPO[il], (1, 2), (3, 4))
        mbash = U * S
        @tensor MPO[il-1][o1 p11; p21 o2] := MPO[il-1][o1 p11 p21 o0] * mbash[o0 o2]
    end
    return MPO
end

#return mpo*mpo'  
function mpo_mpo⁺(l, mpo)
    newmpo = Vector{TensorMap}(undef, l)
    for il = 1:l
        mpolp = permute(mpo[il]',(3,1),(4,2))
        @tensor cpmpo[o1 o1p p1; p1p o2 o2p] := mpo[il][o1 p1; p0 o2] * mpolp[o1p p0; p1p o2p]
        isol = isometry(fuse(codomain(cpmpo, 1) ⊗ codomain(cpmpo, 2)), codomain(cpmpo, 1) ⊗ codomain(cpmpo, 2))
        isor = isometry(domain(cpmpo, 2) ⊗ domain(cpmpo, 3), fuse(domain(cpmpo, 2) ⊗ domain(cpmpo, 3)))
        @tensor newmpo[il][o1n p1; p1p o2n] := isol[o1n; o1 o1p] * cpmpo[o1 o1p p1; p1p o2 o2p] * isor[o2 o2p; o2n]
    end
    return newmpo
end

#return mpo*mpo
function mpo_mpo(l, mpo1, mpo2)
    newmpo = Vector{TensorMap}(undef, l)
    for il = 1:l
        mpolp = mpo2[il]
        @tensor cmpo[o1 o1p p1; p1p o2 o2p] := mpo1[il][o1 p1; p0 o2] * mpolp[o1p p0; p1p o2p]
        isol = isometry(fuse(codomain(cmpo, 1) ⊗ codomain(cmpo, 2)), codomain(cmpo, 1) ⊗ codomain(cmpo, 2))
        isor = isometry(domain(cmpo, 2) ⊗ domain(cmpo, 3), fuse(domain(cmpo, 2) ⊗ domain(cmpo, 3)))
        @tensor newmpo[il][o1n p1; p1p o2n] := isol[o1n; o1 o1p] * cmpo[o1 o1p p1; p1p o2 o2p] * isor[o2 o2p; o2n]
    end
    return newmpo
end

#return constant number mpo
function cmpo(l, n, phySpace, Rm)
    MPO = Vector{TensorMap}(undef, l)
    MPO[1] = n * isometry(Rm(opzo => 1) ⊗ phySpace, phySpace ⊗ Rm(opzo => 1))
    for il = 2:l
        MPO[il] = isometry(Rm(opzo => 1) ⊗ phySpace, phySpace ⊗ Rm(opzo => 1))
    end
    return MPO
end