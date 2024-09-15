#tensor⊕tensor
function tn₊tn(tn1, tn2)
    phy = domain(tn1, 1)
    cd1 = codomain(tn1, 1)
    cd2 = codomain(tn2, 1)
    d1 = domain(tn1, 2)
    d2 = domain(tn2, 2)

    da = d1 ⊕ d2
    cda = cd1 ⊕ cd2

    tnall = TensorMap(zeros, ComplexF64, cda ⊗ phy, phy ⊗ da)
    for (f1, f2) in fusiontrees(tnall)
        flag1 = 0
        for (f11, f12) in fusiontrees(tn1)
            if f1 == f11 && f2 == f12
                flag1 = 1
                s1 = size(tn1[f11, f12])[1]
                s2 = size(tn1[f11, f12])[4]
                tnall[f1, f2][1:s1, 1, 1, 1:s2] = tn1[f11, f12]
            end
        end
        for (f21, f22) in fusiontrees(tn2)
            if f1 == f21 && f2 == f22 && flag1 == 0
                s1 = size(tn2[f21, f22])[1]
                s2 = size(tn2[f21, f22])[4]
                tnall[f1, f2][end-s1+1:end, 1, 1, end-s2+1:end] = tn2[f21, f22]
            end
            if f1 == f21 && f2 == f22 && flag1 == 1
                s1 = size(tn2[f21, f22])[1]
                s2 = size(tn2[f21, f22])[4]
                tnall[f1, f2][end-s1+1:end, 1, 1, end-s2+1:end] = tn2[f21, f22]
            end
        end
    end
    return tnall
end

function op2ten(phySpace)
    dicop = Dict()
    pspace = phySpace[1] #using zeros moment phySpace
    zospace = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 0, 0, 0)

    #"ui"
    top = TensorMap(ones, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(zospace => 1))
    push!(dicop, "ui" => [top, "b", zospace, zospace])
    #"fi"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(zospace => 1))
    blocks(top)[SU2Irrep(0)⊠U1Irrep(0)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= 1
    blocks(top)[SU2Irrep(1 / 2)⊠U1Irrep(1)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= -1
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= 1
    push!(dicop, "fi" => [top, "b", zospace, zospace])
    #"nall"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(zospace => 1))
    blocks(top)[SU2Irrep(1 / 2)⊠U1Irrep(1)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= 1
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= 2
    push!(dicop, "nall" => [top, "b", zospace, zospace])
    #"nud"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(zospace => 1))
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= 1
    push!(dicop, "nud" => [top, "b", zospace, zospace])

    isspace = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](1 / 2, 1, 0, 0)
    ivspace = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](1 / 2, -1, 0, 0)
    #"cdl"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(isspace => 1))
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= -sqrt(2)
    blocks(top)[SU2Irrep(1 / 2)⊠U1Irrep(1)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= 1
    push!(dicop, "cdl" => [top, "f", isspace, zospace])
    #"cl"
    topct = permute(-top', ((3, 1), (4, 2)))
    push!(dicop, "cl" => [topct, "f", ivspace, zospace])

    #"cdfl"
    top = TensorMap(zeros, Rm(zospace => 1) ⊗ pspace ← pspace ⊗ Rm(isspace => 1))
    blocks(top)[SU2Irrep(0)⊠U1Irrep(2)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= sqrt(2)
    blocks(top)[SU2Irrep(1 / 2)⊠U1Irrep(1)⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)] .= 1
    #"cdfr"
    topt = permute(top, ((4, 2), (3, 1)))
    push!(dicop, "cdfr" => [topt, "f", isspace, ivspace])
    #"cfr"
    topct=permute(topt', ((3, 1), (4, 2)))
    push!(dicop, "cfr" => [topct, "f", ivspace, isspace])

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

function setsiteterms(l, mpochain, opzo, opfn)
    #define the storage of operator strings with the form "coefficient, term, position, term, position"
    os = mpochain

    #define the storage for each site with the form "quantumnumberin, qunatumnumberout, coefficient, operators
    #outputposition, inputposition"
    siteterms = Vector(undef, l)

    #define the amount of nodes in siteterms
    dimsiteterms = Vector(undef, l + 1)
    siteterms[1] = [opzo opzo 1 "ui" 1 1]
    dimsiteterms[1] = 1
    for il = 2:l-1
        siteterms[il] = [opzo opzo 1 "ui" 1 1]
        siteterms[il] = [siteterms[il]; [opfn[1] opfn[1] 1 "ui" 2 2]]
        dimsiteterms[il] = 2
    end
    siteterms[l] = [opfn[1] opfn[1] 1 "ui" 2 1]
    dimsiteterms[l] = 2
    dimsiteterms[l+1] = 1
    # o -I-> o -I-> o -I-> o
    #   \                    \
    #    0                     0
    #     -> o -I-> o -I-> o -I-> o
    #initialize the siteterms

    for (inum, it) in enumerate(os)
        # it=os[2]
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

        opac = []
        push!(opac, [opzo, 1])
        for iin = eachindex(opindex)
            opnow = opbash[iin]
            ik1 = lklist[opindex[iin], 1]
            ik2 = lklist[opindex[iin], 2]
            iop = opindex[iin]
            #fill kspace
            ikspace = 0
            if dicop[opnow][3] == Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 0, 0, 0)
                ikspace = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 0, 0, 0)
            elseif opnow == "c" || opnow == "cf"
                ikspace = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 0, -ik1, -ik2)
            else
                ikspace = Irrep[SU₂×U₁×ℤ{nx}×ℤ{ny}](0, 0, ik1, ik2)
            end
            if iop == opindex[1] #use the first row as starting and the second row as the ending
                if iop == opindex[end]
                    siteterms[iop] = [siteterms[iop]; [opac[1][1] opzo coe opnow 1 min(2, dimsiteterms[iop+1])]]
                else
                    dimsiteterms[iop+1] += 1
                    if dicop[opnow][2] == "b"
                        opspace = first(dicop[opnow][3] ⊗ ikspace)
                        siteterms[iop] = [siteterms[iop]; [opac[1][1] opspace coe opnow 1 dimsiteterms[iop+1]]]
                    end
                    if dicop[opnow][2] == "f"
                        opnow = opnow * "l"
                        opspace = first(dicop[opnow][3] ⊗ ikspace)
                        siteterms[iop] = [siteterms[iop]; [opac[1][1] opspace o_parity * coe opnow 1 dimsiteterms[iop+1]]]
                    end
                    opac[1][1] = first(opac[1][1] ⊗ opspace)
                    opac[1][2] = dimsiteterms[iop+1]
                end
            elseif iop != opindex[end]
                if dicop[opnow][2] == "f"
                    if opnow[end] == 'f'
                        opnow = opnow * "l"
                    else
                        opnow = opnow * "r"
                    end
                end
                opspace = first(dicop[opnow][3] ⊗ ikspace)
                acbash = []
                for iac in eachindex(opac)
                    fc = [opac[iac][1] ⊗ opspace...,]
                    fcd = zeros(length(fc))
                    for ifc in eachindex(fc)
                        if fc[ifc][1] != Irrep[SU₂](1) #physical cutoff for final space opzo
                            dimsiteterms[iop+1] += 1
                            siteterms[iop] = [siteterms[iop]; [opac[iac][1] fc[ifc] 1 opnow opac[iac][2] dimsiteterms[iop+1]]]
                            fcd[ifc] = dimsiteterms[iop+1]
                            push!(acbash, [fc[ifc], Int(fcd[ifc])])
                        end
                    end
                end
                opac = acbash
            else
                opnow = opnow * "r"
                opspace = first(dicop[opnow][3] ⊗ ikspace)
                for iac in eachindex(opac)
                    fc = first(opac[iac][1] ⊗ opspace)
                    siteterms[iop] = [siteterms[iop]; [opac[iac][1] fc 1 opnow opac[iac][2] min(2, dimsiteterms[iop+1])]]
                end
            end
        end
    end
    #finish siteterms
    return siteterms
end

function autompo(l, siteterms, phySpace, Rm, opzo, opfn)
    #use siteterms to construct mpo
    # siteterms = st
    MPO = Vector{TensorMap}(undef, l)

    for il = 1:l
        # @show il
        # il=2
        ml = maximum(siteterms[il][:, 5])
        mr = maximum(siteterms[il][:, 6])
        mllist = []
        mrlist = []
        for iml = 1:ml
            sortml = siteterms[il][:, 5] .== iml
            if sum(sortml) == 0
                push!(mllist, opzo)
            else
                push!(mllist, siteterms[il][sortml, 1][1])
            end
        end
        for imr = 1:mr
            sortmr = siteterms[il][:, 6] .== imr
            if sum(sortmr) == 0
                if il == 1
                    push!(mrlist, opfn[imr-1])
                else
                    push!(mrlist, opzo)
                end
            else
                push!(mrlist, siteterms[il][sortmr, 2][1])
            end
        end
        #push mrlist[1] for il=1 to match the dimension
        while il == 1 && mr < length(opfn) + 1
            push!(mrlist, opfn[mr])
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

        MPO[il] = TensorMap(zeros, ComplexF64, blspace ⊗ phySpace[il], phySpace[il] ⊗ brspace)
        for (f1, f2) in fusiontrees(MPO[il])
            dim1 = bld[f1.uncoupled[1]]
            dim4 = brd[f2.uncoupled[2]]
            datab = zeros(Complex{Float64}, dim1, 1, 1, dim4)
            for (i1, id1) in enumerate(dllist[f1.uncoupled[1]])
                for (i4, id4) in enumerate(drlist[f2.uncoupled[2]])
                    if Bool((siteterms[il][:, 5] .== id1)' * (siteterms[il][:, 6] .== id4))
                        po = siteterms[il][siteterms[il][:, 5].==id1.&&siteterms[il][:, 6].==id4, :]
                        phyr = first(dicop[po[4]][4] ⊗ f1.uncoupled[2])
                        if po[4] == "fi" && (f1.coupled[1]⊠f1.coupled[2] == Irrep[SU₂×U₁](0, 2) || f1.coupled[1]⊠f1.coupled[2] == Irrep[SU₂×U₁](0, 0)) #special case for ccd and cdc (one can check "fi" in this space)
                            poc = -po[3]
                        else
                            poc = po[3]
                        end
                        phym = first(block(dicop[po[4]][1], phyr[1]⊠phyr[2]⊠ZNIrrep{nx}(0)⊠ZNIrrep{ny}(0)))
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
        MPO[il], S, V, ϵ = tsvd(MPO[il], (1, 2, 3), (4,); trunc=truncerr(1e-16), alg=TensorKit.SVD())
        # println("MPO $(il) $(il+1): trunc $(ϵ) D $(dim(domain(MPO[il])))")
        MPO[il] = permute(MPO[il], (1, 2), (3, 4))
        mbash = S * V
        @tensor MPO[il+1][o1 p11; p21 o2] := mbash[o1 o0] * MPO[il+1][o0 p11 p21 o2]
    end
    for il = l:-1:2
        U, S, MPO[il], ϵ = tsvd(MPO[il], (1,), (2, 3, 4); trunc=truncerr(1e-16), alg=TensorKit.SVD())
        # println("MPO $(il) $(il-1): trunc $(ϵ) D $(dim(codomain(MPO[il])))")
        MPO[il] = permute(MPO[il], (1, 2), (3, 4))
        mbash = U * S
        @tensor MPO[il-1][o1 p11; p21 o2] := MPO[il-1][o1 p11 p21 o0] * mbash[o0 o2]
    end
    return MPO
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

#return mpo*mpo'  
function mpo_mpo⁺(l, mpo)
    newmpo = Vector{TensorMap}(undef, l)
    for il = 1:l
        mpolp = permute(mpo[il]', (3, 1), (4, 2))
        @tensor cnmpo[o1 o1p p1; p1p o2 o2p] := mpo[il][o1 p1; p0 o2] * mpolp[o1p p0; p1p o2p]
        isol = isometry(fuse(codomain(cnmpo, 1) ⊗ codomain(cnmpo, 2)), codomain(cnmpo, 1) ⊗ codomain(cnmpo, 2))
        isor = isometry(domain(cnmpo, 2) ⊗ domain(cnmpo, 3), fuse(domain(cnmpo, 2) ⊗ domain(cnmpo, 3)))
        @tensor newmpo[il][o1n p1; p1p o2n] := isol[o1n; o1 o1p] * cnmpo[o1 o1p p1; p1p o2 o2p] * isor[o2 o2p; o2n]
    end
    MPO = newmpo
    #finish MPO and use svd to compress
    for il = l:-1:2
        U, S, MPO[il], ϵ = tsvd(MPO[il], (1,), (2, 3, 4); trunc=truncerr(1e-10), alg=TensorKit.SVD())
        # println("MPO $(il) $(il-1): trunc $(ϵ) D $(dim(codomain(MPO[il])))")
        MPO[il] = permute(MPO[il], (1, 2), (3, 4))
        mbash = U * S
        @tensor MPO[il-1][o1 p11; p21 o2] := MPO[il-1][o1 p11 p21 o0] * mbash[o0 o2]
    end
    for il = 1:l-1
        MPO[il], S, V, ϵ = tsvd(MPO[il], (1, 2, 3), (4,); trunc=truncerr(1e-10), alg=TensorKit.SVD())
        # println("MPO $(il) $(il+1): trunc $(ϵ) D $(dim(domain(MPO[il])))")
        MPO[il] = permute(MPO[il], (1, 2), (3, 4))
        mbash = S * V
        @tensor MPO[il+1][o1 p11; p21 o2] := mbash[o1 o0] * MPO[il+1][o0 p11 p21 o2]
    end
    return newmpo
end

#return mpo+mpo
function mpo₊mpo(l, mpo1, mpo2)
    # mpo1=mpo
    MPO = Vector{TensorMap}(undef, l)
    for il = 1:l
        MPO[il] = tn₊tn(mpo1[il], mpo2[il])
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

#return constant number mpo
function cmpo(l, n, phySpace, Rm)
    MPO = Vector{TensorMap}(undef, l)
    MPO[1] = n * isometry(Rm(opzo => 1) ⊗ phySpace[1], phySpace[1] ⊗ Rm(opzo => 1))
    for il = 2:l
        MPO[il] = isometry(Rm(opzo => 1) ⊗ phySpace[il], phySpace[il] ⊗ Rm(opzo => 1))
    end
    return MPO
end