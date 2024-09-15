#For two-site measurment s⋅s
using StridedViews
using Strided
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

function autompo(l, mpochain, phySpace)

    #define the storage of operator strings with the form "coefficient, term, position, term, position"
    os = mpochain

    #define the storage for each site with the form "quantumnumberin, qunatumnumberout, coefficient, operators
    #outputposition, inputposition, J_W_matrix"
    siteterms = Vector(undef, l)
    #define the amount of nodes in siteterms
    dimsiteterms = Vector(undef, l + 1)
    siteterms[1] = [dicqn["ui"] dicqn["ui"] 1 "ui" 1 1 "ui"]
    dimsiteterms[1] = 1
    for il = 2:l-1
        siteterms[il] = [dicqn["ui"] dicqn["ui"] 1 "ui" 1 1 "ui"]
        siteterms[il] = [siteterms[il]; [dicqn["ui"] dicqn["ui"] 1 "ui" 2 2 "ui"]]
        dimsiteterms[il] = 2
    end
    siteterms[l] = [dicqn["ui"] dicqn["ui"] 1 "ui" 2 1 "ui"]
    dimsiteterms[l] = 2
    dimsiteterms[l+1] = 1
    # o -I-> o -I-> o -I-> o
    #   \                    \
    #    0                     0
    #     -> o -I-> o -I-> o -I-> o 
    #initialize the siteterms

    for (inum, it) in enumerate(os)
        # it=os[8]
        coe = it[1]
        nop = (length(it) - 1) ÷ 2
        opbash = []
        opindex = []
        opqn = []
        for inop = 1:nop
            push!(opbash, it[inop*2])
            push!(opindex, it[inop*2+1])
            push!(opqn, dicqn[it[inop*2]])
        end
        isort = sortperm(opindex)
        o_parity = levicivita(isort)
        opindex = opindex[isort]
        opbash = opbash[isort]
        opqn = opqn[isort]

        opi = []
        opb = []
        opq = []
        for iu in unique(opindex)
            if sum(opindex .== iu) > 1
                st = join(opbash[opindex.==iu])
                if st == "cucu" || st == "cdcd" || st == "cducdu" || st == "cddcdd"
                    coe=0
                else
                    push!(opi, iu)
                    push!(opb, st)
                    push!(opq, dicqn[st])
                end
            else
                push!(opi, iu)
                push!(opb, opbash[opindex.==iu][1])
                push!(opq, opqn[opindex.==iu][1])
            end
        end
        opindex = opi
        opbash = opb
        opqn = opq

        if coe!=0
            iopt = 0
            iopf = 0
            for iop = opindex[1]:opindex[end]
                if !(iop in opindex)
                    if (-1)^iopf == 1
                        push!(opbash, "ui")
                        push!(opindex, iop)
                        push!(opqn, dicqn["ui"])
                    else
                        push!(opbash, "fi")
                        push!(opindex, iop)
                        push!(opqn, dicqn["fi"])
                    end
                else
                    iopt += 1
                    if dicbf[opbash[iopt]] == "f"
                        iopf += 1
                    end
                end
            end
            isort = sortperm(opindex)
            opindex = opindex[isort]
            opbash = opbash[isort]
            opqn = opqn[isort]
        else
            continue
        end
        #finish operator string and use operator string to construct siteterms

        opac = dicqn["ui"]
        for iin = eachindex(opindex)
            iop = opindex[iin]
            if iop == opindex[1] #use the first row as starting and the second row as the ending
                if iop == opindex[end]
                    siteterms[iop] = [siteterms[iop]; [dicqn["ui"] dicqn["ui"] coe opbash[iin] 1 min(2, dimsiteterms[iop+1]) "ui"]]
                else
                    dimsiteterms[iop+1] += 1
                    if dicbf[opbash[iin]] == "b"
                        siteterms[iop] = [siteterms[iop]; [dicqn["ui"] opqn[iin] coe opbash[iin] 1 dimsiteterms[iop+1] "ui"]]
                    end
                    if dicbf[opbash[iin]] == "f"
                        siteterms[iop] = [siteterms[iop]; [dicqn["ui"] opqn[iin] o_parity * coe opbash[iin] 1 dimsiteterms[iop+1] "fi"]]
                    end
                    opac = first(opac ⊗ opqn[iin])
                end
            elseif iop != opindex[end]
                dimsiteterms[iop+1] += 1
                if dicbf[opbash[iin]] == "f" && opbash[iin+1] == "fi"
                    siteterms[iop] = [siteterms[iop]; [opac first(opac ⊗ opqn[iin]) 1 opbash[iin] dimsiteterms[iop] dimsiteterms[iop+1] "fi"]]
                else
                    siteterms[iop] = [siteterms[iop]; [opac first(opac ⊗ opqn[iin]) 1 opbash[iin] dimsiteterms[iop] dimsiteterms[iop+1] "ui"]]
                end
                opac = first(opac ⊗ opqn[iin])
            else
                siteterms[iop] = [siteterms[iop]; [opac dicqn["ui"] 1 opbash[iin] dimsiteterms[iop] min(2, dimsiteterms[iop+1]) "ui"]]
            end
        end
    end
    #finish siteterms

    #use siteterms to construct mpo
    MPO = Vector{TensorMap}(undef, l)

    for il = 1:l
        # il=1
        ml = maximum(siteterms[il][:, 5])
        mr = maximum(siteterms[il][:, 6])
        mllist = Vector(undef, ml)
        mrlist = Vector(undef, mr)
        for iml = 1:ml
            sortml = siteterms[il][:, 5] .== iml
            if sum(sortml) == 0
                mllist[iml] = dicqn["ui"]
            else
                mllist[iml] = siteterms[il][sortml, 1][1]
            end
        end
        for imr = 1:mr
            sortmr = siteterms[il][:, 6] .== imr
            if sum(sortmr) == 0
                mrlist[imr] = dicqn["ui"]
            else
                mrlist[imr] = siteterms[il][sortmr, 2][1]
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

        blspace = Rep[U₁×U₁](bld)
        brspace = Rep[U₁×U₁](brd)

        MPO[il] = TensorMap(zeros, ComplexF64, blspace ⊗ phySpace, phySpace ⊗ brspace)
        for (f1, f2) in fusiontrees(MPO[il])
            phyr = dicphys[f1.uncoupled[2]]
            phyc = dicphys[f2.uncoupled[1]]
            dim1 = bld[f1.uncoupled[1]]
            dim4 = brd[f2.uncoupled[2]]
            datab = zeros(Complex{Float64}, dim1, 1, 1, dim4)
            for (i1, id1) in enumerate(dllist[f1.uncoupled[1]])
                for (i4, id4) in enumerate(drlist[f2.uncoupled[2]])
                    if ((siteterms[il][:, 5] .== id1)' * (siteterms[il][:, 6] .== id4)) != 0
                        po = siteterms[il][siteterms[il][:, 5].==id1.&&siteterms[il][:, 6].==id4, :]
                        @show
                        phym = dicop[po[4]] * dicop[po[7]]
                        datab[i1, 1, 1, i4] = po[3] * phym[phyr, phyc]
                    else
                        datab[i1, 1, 1, i4] = 0
                    end
                end
            end
            MPO[il][f1, f2] = datab
        end
    end

    #finish MPO and use svd to compress
    # for il = 1:l-1
    #     MPO[il], S, V, ϵ = tsvd(MPO[il], (1, 2, 3), (4,); trunc=truncdim(mpod), alg=TensorKit.SVD())
    #     println("MPO $(il) $(il+1): trunc $(ϵ) D $(dim(domain(MPO[il])))")
    #     MPO[il] = permute(MPO[il], (1, 2), (3, 4))
    #     mbash = S * V
    #     @tensor MPO[il+1][o1 p11; p21 o2] := mbash[o1 o0] * MPO[il+1][o0 p11 p21 o2]
    # end

    return MPO
end