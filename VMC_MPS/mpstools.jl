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
function observable_config(ostring, md, co0, config0, Dall)
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
    return osum/co0, Dc
end

# ---------------------------------------------------------------------------
# MPS梯度下降：用采样构型计算MPS每个参数的自然梯度
# ---------------------------------------------------------------------------

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
        (wk, _, _, enm, gradm) = Dm[cf]
        s=sqrt(wk/wm)
        e[k]=s*(enm-em)
        for il=1:l
            g=gradm[il]
            g0=gm[il]
            ix=idx[il]
            o=offs[il]
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

# ---------------------------------------------------------------------------
# 残差扩维：用采样构型经 H 散射得到的新构型扩充 MPS 的键维
# ---------------------------------------------------------------------------

#稠密指标对应的对称性扇区：GradedSpace 的稠密排布按 sectors(V) 的顺序，每个扇区连续占 dim(V, c) 个位置
dense_sectors(V) = [c for c in sectors(V) for _ in 1:dim(V, c)]

#一步虚时演化 |T⟩ = |ψ⟩ - Δτ (H-E)|ψ⟩ = cψ|ψ⟩ - Δτ H|ψ⟩，cψ = 1 + Δτ E
#H|ψ⟩ 由采样给出：h(b) = ‖ψ‖² ĥ(b)，只在 hb 的构型上非零；-E|ψ⟩ 在全空间上精确保留在 cψ|ψ⟩ 中
#从左到右逐键做 TT-SVD：第 i 步已有新的左正交基 B[1..i-1]，把 |T⟩ 投影到 span{B}⊗(第 i 个格点) 上，
#对约化密度矩阵 ρ = T T† 按扇区对角化，保留 min(D_i + Δb, maxdim) 个最大本征态作为新的 B[i]。
#ψ 取中心在 i 的混合正则形式 ψ = AL[1..i-1] AC[i] AR[i+1..l]，右侧基 |R_β⟩ 正交归一，于是
#  ρ = Y Y† - Δτ (Y G† + G Y†) + Δτ² V V†
#  Y[(α,s),β] = cψ Σ_α' E[α,α'] AC[α',s,β]，               E[α,α'] = ⟨B_α|AL_α'⟩
#  G[(α,s),β] = Σ_k h_k conj(λ_k[α]) δ_{s,s_k} conj(R_k[β])，λ_k = 构型 k 前缀在 B 下的行向量，R_k = 后缀在 AR 下的列向量
#  V[(α,s),σR] = Σ_{k: 后缀为 σR} h_k conj(λ_k[α]) δ_{s,s_k}
function residual_expand(mt, hb, em, Δτ, Δb, maxdim; tol=1e-14)
    ψ=FiniteMPS([ComplexF64(1)*copy(mt[il]) for il=1:l])
    ALd=[convert(Array, ψ.AL[il]) for il=1:l]
    ARd=[convert(Array, ψ.AR[il]) for il=1:l]
    ACd=[convert(Array, ψ.AC[il]) for il=1:l]
    d=dim(phySpace)
    psec=dense_sectors(phySpace)
    cψ=1 + Δτ*real(em)
    nrm2=norm(ψ)^2

    cbs=collect(keys(hb))
    nk=length(cbs)
    rv=[nrm2*hb[cb] for cb in cbs]
    ps=[[s2p[cb[2*il-1], cb[2*il]] for il=1:l] for cb in cbs]

    #R[k][il]：构型 k 在格点 il 右侧的后缀在 ψ 右正交基下的列向量（长度 = 第 il 个键的维数）
    R=[Vector{Vector{ComplexF64}}(undef, l) for _ in 1:nk]
    for k=1:nk
        R[k][l]=ComplexF64[1]
        for il=l:-1:2
            R[k][il-1]=ARd[il][:, ps[k][il], :] * R[k][il]
        end
    end

    λ=[ComplexF64[1] for _ in 1:nk]   #构型前缀在新基 B 下的行向量 (存为列向量)
    E=ones(ComplexF64, 1, 1)          #E[α,α'] = ⟨B_α|AL_α'⟩
    Vl=codomain(ψ.AL[1])[1]
    Bs=Vector{TensorMap}(undef, l)
    for il=1:l
        Dl=size(E, 1)
        Dr=size(ACd[il], 3)
        Y=cψ * reshape(E * reshape(ACd[il], size(ACd[il], 1), :), Dl*d, Dr)

        #U[:, k] = h_k conj(λ_k) ⊗ e_{s_k}
        U=zeros(ComplexF64, Dl*d, nk)
        for k=1:nk
            o=(ps[k][il]-1)*Dl
            U[(o+1):(o+Dl), k] .= rv[k] .* conj.(λ[k])
        end
        G=U * reduce(hcat, (R[k][il] for k=1:nk))'

        if il == l
            #最后一个格点右键维数为 1，直接把投影后的目标态放进去
            Bd=reshape(Y - Δτ * G, Dl, d, Dr)
            Vr=domain(ψ.AL[l])[1]
            Bs[il]=TensorMap(Bd, Vl ⊗ phySpace ← Vr)
            break
        end

        #同一后缀的构型相干叠加：V = U S，S 为构型→后缀的指示矩阵
        suf=Dict{Vector{Int},Int}()
        sidx=[get!(suf, cbs[k][(2*il+1):end], length(suf)+1) for k=1:nk]
        Vm=zeros(ComplexF64, Dl*d, length(suf))
        for k=1:nk
            Vm[:, sidx[k]] .+= @view U[:, k]
        end
        YG=Y * G'
        ρ=Y * Y' - Δτ * (YG + YG') + Δτ^2 * (Vm * Vm')

        #(α,s) 行的扇区 = sector(α) ⊗ sector(s)；ρ 在扇区上块对角，逐块对角化
        lsec=dense_sectors(Vl)
        rowsec=[first(lsec[a] ⊗ psec[s]) for s=1:d for a=1:Dl]
        cands=Tuple{Float64,eltype(rowsec),Vector{ComplexF64}}[]
        rows=Dict{eltype(rowsec),Vector{Int}}()
        for c in unique(rowsec)
            ix=findall(==(c), rowsec)
            rows[c]=ix
            F=eigen(Hermitian(ρ[ix, ix]))
            for j in eachindex(F.values)
                push!(cands, (F.values[j], c, F.vectors[:, j]))
            end
        end
        sort!(cands; by=x -> -x[1])
        λmax=max(cands[1][1], eps())
        nkeep=min(Dr + Δb, maxdim, count(x -> x[1] > tol*λmax, cands))
        kept=cands[1:nkeep]

        #新的右键空间；稠密列按 sectors(Vr) 的顺序排列
        cnt=Dict{eltype(rowsec),Int}()
        for (_, c, _) in kept
            cnt[c]=get(cnt, c, 0)+1
        end
        Vr=Vect[FermionParity⊠U1Irrep]((c => nc for (c, nc) in cnt)...)
        Bm=zeros(ComplexF64, Dl*d, nkeep)
        j=0
        for c in sectors(Vr)
            for (_, c1, v) in kept
                c1 == c || continue
                j+=1
                Bm[rows[c], j] .= v
            end
        end
        Bd=reshape(Bm, Dl, d, nkeep)
        Bs[il]=TensorMap(Bd, Vl ⊗ phySpace ← Vr)

        #更新左环境：E ← Σ_s B_s† E AL_s，λ_k ← λ_k B_{s_k}
        Enew=zeros(ComplexF64, nkeep, size(ALd[il], 3))
        for s=1:d
            Enew+=Bd[:, s, :]' * E * ALd[il][:, s, :]
        end
        E=Enew
        for k=1:nk
            λ[k]=transpose(Bd[:, ps[k][il], :]) * λ[k]
        end
        Vl=Vr
    end
    return normalize!(FiniteMPS([Bs[il] for il=1:l]))
end
