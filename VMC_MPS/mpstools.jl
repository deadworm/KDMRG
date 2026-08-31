#For initial random mps
function rnd_mps(l, bond_dim)
    B = Vector{typeof(zb)}(undef, l + 1)
    B[1] = zb
    # 构造中间虚键：B[i + 1] 对应第 i 个站点之后的键
    for i in 1:(l-1)
        qmin = max(0, n - 2 * (l - i))
        qmax = min(2 * i, n)
        charges = collect(qmin:qmax)
        nsectors = length(charges)
        # 将最大维数尽量平均分配给各个允许的电荷扇区
        base, remainder = divrem(bond_dim, nsectors)
        sector_dims = [
            (q % 2, q) => base + (j <= remainder ? 1 : 0)
            for (j, q) in enumerate(charges)
        ]
        B[i+1] = Vect[FermionParity⊠U1Irrep](sector_dims...)
    end
    # 固定右边界，从而保证整个 MPS 的总粒子数为 n
    B[end] = nb
    tensors = Vector{TensorMap}(undef, l)
    for i in 1:l
        tensors[i] = randn(ComplexF64, B[i] ⊗ phySpace ← B[i+1])
    end
    m = FiniteMPS([tensors[i] for i in 1:l])
    m = changebonds(m, RandExpand(trscheme=truncdim(bond_dim)))
    m = changebonds(m, SvdCut(trscheme=truncdim(bond_dim)))
    @show [dim(space(m[i], 3)) for i in 1:l]
    return m
end

#For initial direct-product mps
function prod_mps(config)
    B = Vector{typeof(zb)}(undef, l + 1)
    B[1] = zb
    tensors = Vector{TensorMap}(undef, l)
    for i = 1:l
        pindex=s2p[(config[2*i-1],config[2*i])]
        B[i+1] = fuse(B[i], phylist[pindex])
        tenow = zeros(Float64, B[i] ⊗ phySpace ← B[i+1])
        if pindex == 3
            tenow.data[2] = 1.0
        else
            tenow.data[1] = 1.0
        end
        tensors[i]=tenow
    end
    mf = FiniteMPS([tensors[i] for i in 1:l])
    return mf
end

#For computing observables [coefficient, site1, site2, ...] for each term
function observable_config(ostring, m, cm, config, Dall)
    oterms=size(ostring, 1)
    olength=size(ostring, 2)
    oc=ostring[:, 1]
    Dc = Dict{Vector{Int},ComplexF64}()
    for it=1:oterms
        nconfig=copy(config)
        for il=0:(olength-2)
            oind=Int(ostring[it, end-il])
            if oind==0
                continue
            elseif oind>0&&nconfig[oind]==0
                oc[it]=oc[it]*(-1)^sum(@view nconfig[1:(oind-1)])
                nconfig[oind]=1
            elseif oind<0&&nconfig[-oind]==1
                oc[it]=oc[it]*(-1)^sum(@view nconfig[1:(-oind-1)])
                nconfig[-oind]=0
            else
                oc[it]=0
                break
            end
        end
        if oc[it]!=0
            if !haskey(Dc, nconfig)
                Dc[nconfig]=oc[it]
            else
                Dc[nconfig]+=oc[it]
            end
        end
    end
    osum=0
    for cf in keys(Dc)
        if haskey(Dall, cf)
            osum+=Dc[cf] * Dall[cf]
        else
            mf=prod_mps(cf)
            Dall[cf] = MPSKit.dot(mf, m)
            osum+=Dc[cf] * Dall[cf]
        end
    end
    olocal=osum / cm
    return olocal
end

#For mps gradient-normalize update
function gn_update(m, wm, em, gm, Dm, dm)
    nm=Vector{TensorMap}(undef, l)
    configs=collect(keys(Dm))
    for il=1:l
        nm[il] = copy(m[il])
    end
    for cf in configs
        (cn, _, enm, gradm) = Dm[cf]
        for il=1:l
            δm = cn * (gradm[il]-gm[il]) * (enm-em) / wm
            nm[il] = nm[il] - dm * δm
        end
    end
    δm_norm=0
    for il=1:l
        δm_norm+=norm(nm[il]-m[il])
    end
    @show δm_norm
    return FiniteMPS([nm[il] for il in 1:l])
end

#For mps stochastic-reconfiguration update
function sr_update(m, wm, em, gm, Dm, dt)
    lk=length(Dm)
    configs=collect(keys(Dm))
    cflist=zeros(ComplexF64, lk, lk)
    enlist = [sqrt(Dm[cf][1]/wm) * (Dm[cf][3]-em) for cf in configs]
    glist = [sqrt(Dm[cf][1]/wm) * (Dm[cf][4]-gm) for cf in configs]
    for i1=1:lk
        for i2=1:i1
            cft=0
            if i2 == i1
                for il=1:l
                    @planar cfb = glist[i1][il]'[r; l u] * glist[i2][il][l u; r]
                    cft+=cfb
                end
                cflist[i1, i2] = real(cft)
            else
                for il=1:l
                    @planar cfb = glist[i1][il]'[r; l u] * glist[i2][il][l u; r]
                    cft+=cfb
                end
                cflist[i1, i2] = cft
                cflist[i2, i1] = conj(cft)
            end
        end
    end
    ylist, cglog = cg(Hermitian(cflist+1e-4*sum(diag(cflist))/lk*I), enlist; reltol=1e-5, maxiter=5lk, log=true)
    @show cglog
    nm=Vector{TensorMap}(undef, l)
    for il=1:l
        nm[il] = copy(m[il])
    end
    ∂m_norm=0
    for il=1:l
        for ik=1:lk
            ∂m = glist[ik][il] * ylist[ik]
            nm[il] = nm[il] - dt * ∂m
        end
        ∂m_norm+=norm(nm[il]-m[il])
    end
    @show ∂m_norm
    return FiniteMPS([nm[il] for il in 1:l])
end

#For computing the entanglement entropy
function mps_entropy(m, l)
    entropy_list = Vector{Float64}(undef, l-1)
    for b in 1:(l-1)
        λ = entanglement_spectrum(m, b)
        # collect singular values from all symmetry sectors
        λ_all = vcat(values(λ)...)
        # probabilities
        p = λ_all .^ 2
        p ./= sum(p)   # safe normalization
        # von Neumann entropy
        S = -sum(p .* log.(p .+ 1e-15))
        entropy_list[b] = S
    end
    return entropy_list
end