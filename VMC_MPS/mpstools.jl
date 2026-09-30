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

#For expanding the active bond of mps
function expand_active(m, active; add=4, noise=1e-3)
    # Frozen exterior bonds have rank one, with fixed cumulative U(1) charges.
    sub=FiniteMPS([ComplexF64(1)*copy(m[j]) for j in active])
    sub=changebonds(sub, RandExpand(trscheme=truncdim(add)))
    raw=[ComplexF64(1)*copy(m[j]) for j in 1:length(m)]
    for (k, j) in enumerate(active)
        raw[j]=copy(sub[k])
        raw[j].data .+= noise*randn(ComplexF64, length(raw[j].data))/sqrt(length(raw[j].data))
    end
    return normalize!(FiniteMPS(raw))
end

#For initial direct-product mps
function prod_mps(config0)
    B = Vector{typeof(zb)}(undef, l + 1)
    B[1] = zb
    tensors = Vector{TensorMap}(undef, l)
    for i = 1:l
        pindex=s2p[(config0[2*i-1], config0[2*i])]
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

#For computing the slice of mps
function mps_slice(md, config0)
    v=ComplexF64[1]
    for il=1:l
        v=v*md[il][:, s2p[config0[(il-1)*2+1], config0[(il-1)*2+2]], :]
    end
    return v[1]
end

#For computing observables [coefficient, site1, site2, ...] for each term
function observable_config(ostring, md, cm, config0, Dall)
    oterms=size(ostring, 1)
    olength=size(ostring, 2)
    oc=ostring[:, 1]
    Dc = Dict{Vector{Int},ComplexF64}()
    for it=1:oterms
        nconfig=copy(config0)
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
            Dall[cf] = mps_slice(md, cf)
            osum+=Dc[cf] * Dall[cf]
        end
    end
    olocal=osum / cm
    return olocal
end

#稠密数组中对称性允许的元素的线性指标（其余位置恒为 0，不是独立参数）
function sr_param_index(t)
    u=zeros(ComplexF64, space(t))
    u.data .= 1
    idx=findall(!iszero, vec(convert(Array, u)))
    length(idx) == dim(space(t)) || error("sr_param_index: found $(length(idx)) entries, expected $(dim(space(t)))")
    return idx
end

#求解 (A + shift*I) x = b，A 为 Hermitian 半正定；优先 Cholesky，失败时退回 Bunch-Kaufman
function sr_solve(A, shift, b)
    H=Hermitian(A + shift*I)
    F=cholesky(H; check=false)
    return issuccess(F) ? F \ b : H \ b
end

#For mps stochastic-reconfiguration update
#regularization: 对角偏移 λ 相对于 S 对角元平均值 tr(S)/Np 的比例
function sr_update(mt, wm, em, gm, Dm, dt; regularization=1e-2)
    lk=length(Dm)
    lk > 0 || error("sr_update: no samples")
    configs=collect(keys(Dm))
    idx=[sr_param_index(mt[il]) for il=1:l]
    offs=cumsum([0; length.(idx)])
    np=offs[end]

    # Y[:, k] = sqrt(w_k/wm) * (O_k - <O>)，只取对称性允许的参数；e[k] = sqrt(w_k/wm) * (E_k - <E>)
    Y=Matrix{ComplexF64}(undef, np, lk)
    e=Vector{ComplexF64}(undef, lk)
    for (k, cf) in enumerate(configs)
        (wk, _, enm, gradm) = Dm[cf]
        s=sqrt(wk/wm)
        e[k]=s*(enm-em)
        for il=1:l
            g=gradm[il]; g0=gm[il]; ix=idx[il]; o=offs[il]
            @inbounds for j in eachindex(ix)
                Y[o+j, k]=s*(g[ix[j]]-g0[ix[j]])
            end
        end
    end

    # S = Y Y'（np×np），F = Y e；δ = (S+λ)^-1 F = Y (Y'Y+λ)^-1 e。在较小的空间里求解
    if np <= lk
        S=Y*Y'
        shift=regularization*max(real(tr(S))/np, eps())
        δ=sr_solve(S, shift, Y*e)
        solved="parameter space"
    else
        T=Y'*Y
        shift=regularization*max(real(tr(T))/np, eps())
        δ=Y*sr_solve(T, shift, e)
        solved="sample space"
    end
    all(isfinite, δ) || error("sr_update: non-finite update (shift=$shift)")
    println("SR: parameter space Np=", np, " (dense ", sum(length, gm), "), sample space lk=", lk,
        ", solved in ", solved, ", shift=", shift)

    nm=map(1:l) do il
        ∂m=zeros(ComplexF64, size(gm[il]))
        ∂m[idx[il]] .= @view δ[(offs[il]+1):offs[il+1]]
        mt[il] - dt * TensorMap(∂m, space(mt[il]))
    end
    @show [norm(nm[il]-mt[il]) for il=1:l]
    return FiniteMPS(nm)
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
