#Construct Hamiltonian
function rs_Hamiltonian()
    ls=2*l
    enstring=zeros(Float64, 3*ls, 5)
    for i=1:(ls-1)
        enstring[(i-1)*3+1, :]=[-t i -i-1 0 0]
        enstring[(i-1)*3+2, :]=[-t i+1 -i 0 0]
        enstring[(i-1)*3+3, :]=[V i -i i+1 -i-1]
    end
    enstring[(ls-1)*3+1, :]=[t ls -1 0 0]
    enstring[(ls-1)*3+2, :]=[t 1 -ls 0 0]
    enstring[(ls-1)*3+3, :]=[V ls -ls 1 -1]
    return enstring
end

function ms_Hamiltonian()
    ls=2*l
    k2p = [x <= l ? 2x - 1 : 2 * (2l - x + 1) for x in 1:ls]
    enstring=zeros(Float64, ls^3 + ls, 5)
    for ik1=1:ls
        for ik2=1:ls
            for iq=1:ls
                enstring[(ik1-1)*ls^2+(ik2-1)*ls+iq, :]=[V/ls*cos(2*π/ls*(iq)) k2p[ik1] -k2p[mod(ik1+iq-1, ls)+1] k2p[ik2] -k2p[mod(ik2-iq-1, ls)+1]]
            end
        end
    end
    for i=1:l
        enstring[ls^3+(i-1)*2+1, :]=[-2*t*cos(2*π/ls*(i-1+1/2)) (i-1)*2+1 -(i-1)*2-1 0 0]
        enstring[ls^3+(i-1)*2+2, :]=[-2*t*cos(2*π/ls*(i-1+1/2)) (i-1)*2+2 -(i-1)*2-2 0 0]
    end
    return enstring
end

function Vk(s, m1, m2, m0)
    vk=0
    for k=0:2s
        if -k<=m0<=k
            vk+=(U0*(2k+1)-U1*k*(k+1)*(2k+1))*wigner3j(s, k, s, -m1, -m0, m1+m0)*wigner3j(s, k, s, -m2, m0, m2-m0)*wigner3j(s, k, s, -s, 0, s)^2
        end
    end
    return vk
end

function fsIsing_Hamiltonian()
    enstring=zeros(Float64, 0, 5)
    for im1=1:l
        for im2=1:l
            for m0=-l+1:l-1
                m1=-s+im1-1
                m2=-s+im2-1
                if -s<=m1+m0<=s&&-s<=m2-m0<=s
                    enstring=vcat(enstring, [(-1)^(2s+m0+m1+m2)*(2s+1)^2*Vk(s, m1, m2, m0) (im1-1)*2+1 -(m1+m0+s)*2-1 (im2-1)*2+2 -(m2-m0+s)*2-2])
                    enstring=vcat(enstring, [(-1)^(2s+m0+m1+m2)*(2s+1)^2*Vk(s, m1, m2, m0) (im1-1)*2+2 -(m1+m0+s)*2-2 (im2-1)*2+1 -(m2-m0+s)*2-1])
                end
            end
        end
    end
    for i=1:l
        enstring=vcat(enstring, [-h (i-1)*2+1 -(i-1)*2-2 0 0])
        enstring=vcat(enstring, [-h (i-1)*2+2 -(i-1)*2-1 0 0])
    end
    return enstring
end

#For the mc sampling
function mc_sample(m, config0, mcsample, α)
    Dall = Dict{Vector{Int},ComplexF64}()
    Dm = Dict{Vector{Int},Tuple{Float64,Float64,ComplexF64,ComplexF64,Vector{Any},Dict{Vector{Int},ComplexF64}}}()
    # 采样会访问 m.AL/m.AR 并移动 gauge center，之后 m[il] 会换成另一套规范，
    # 所以先固定一份张量：梯度和参数更新都必须基于这同一组张量
    mt=[copy(m[il]) for il=1:l]
    As=dense_slices(mt)

    # 规范变换只做一次，采样循环里不再反复访问 m.AL/m.AR
    ARd=dense_slices([m.AR[il] for il=1:l])
    ALd=dense_slices([m.AL[il] for il=1:l])

    updt=0
    totl=0
    for imc=1:mcsample
        if imc%2==1
            (config0, wt, updt1, totl1)=ldirect_sample(ARd, config0, α)
        else
            (config0, wt, updt1, totl1)=rdirect_sample(ALd, config0, α)
        end
        updt+=updt1
        totl+=totl1

        #measure observables here
        if !haskey(Dm, config0)
            co0=mps_amp(As, config0)
            (enm, Dc)=observable_config(enstring, As, co0, config0, Dall)
            gradm=measure_grad(As, co0, config0)
            Dm[config0] = (wt/mcsample, wt^2/mcsample, co0, enm, gradm, Dc)
        else
            Dm[config0] = (Dm[config0][1]+wt/mcsample, Dm[config0][2]+wt^2/mcsample, Dm[config0][3], Dm[config0][4], Dm[config0][5], Dm[config0][6])
        end
    end
    wm = 0
    wm2 = 0
    for cf in keys(Dm)
        wm += Dm[cf][1]
        wm2 += Dm[cf][2]
    end
    em = 0
    for cf in keys(Dm)
        em += Dm[cf][1]/wm * Dm[cf][4]
    end
    gm = Vector{Any}(undef, l)
    for il=1:l
        gm[il] = zeros(ComplexF64, slices_size(As[il]))
        for cf in keys(Dm)
            gm[il] += Dm[cf][1]/wm * Dm[cf][5][il]
        end
    end
    #ĥ(b) = conj( Σ_σ (w_σ/wm) * H_bσ / ψ(σ) ) ≈ ⟨b|H|ψ⟩/⟨ψ|ψ⟩（H 为实对称矩阵）
    #-E|ψ⟩ 部分不在这里减：它在全空间上非零，由 residual_expand 在 MPS 层面精确处理
    hb = Dict{Vector{Int},ComplexF64}()
    for cs in keys(Dm)
        for cb in keys(Dm[cs][6])
            hb[cb] = get(hb, cb, 0.0im) + Dm[cs][1] / wm * Dm[cs][6][cb] / Dm[cs][3]
        end
    end
    for cb in keys(hb)
        hb[cb] = conj(hb[cb])
    end

    println("average energy: ", em)
    println("flip ratio: ", updt/totl)
    println("effective sample ratio: ", wm^2/wm2)
    return wm, em, gm, Dm, mt, hb
end

#稠密切片对应的三维数组尺寸 (Dl, d, Dr)
slices_size(Asl) = (size(Asl[1], 1), length(Asl), size(Asl[1], 2))

#For energy gradient measurement
#As[il][c] 为稠密切片（见 dense_slices）
function measure_grad(As, co0, config0)
    vl = Vector{Any}(undef, l + 1)
    vr = Vector{Any}(undef, l + 1)
    vl[1] = ComplexF64[1]
    for il = 2:(l+1)
        vl[il] = As[il-1][s2p[config0[(il-2)*2+1], config0[(il-2)*2+2]]]' * vl[il-1]
    end
    vr[l+1] = ComplexF64[1]
    for il = l:-1:1
        vr[il] = vr[il+1] * As[il][s2p[config0[(il-1)*2+1], config0[(il-1)*2+2]]]'
    end

    mg=Vector{Any}(undef, l)
    for il=1:l
        mb=zeros(ComplexF64, slices_size(As[il]))
        mb[:, s2p[config0[(il-1)*2+1], config0[(il-1)*2+2]], :]=kron(vl[il], vr[il+1])
        mg[il]=mb/conj(co0)
    end
    return mg
end

#For configuration direct sampling
#采样用的稠密切片：As[il][c] 为第 il 个格点、物理态 c（与 s2p/p2s 编号一致）对应的 Dl×Dr 矩阵
function dense_slices(A)
    return [[Matrix{ComplexF64}(a[:, c, :]) for c=1:size(a, 2)] for a in (convert(Array, t) for t in A)]
end

#按 p^α 从单点条件概率 dmz 中抽样，返回 (物理态 c, 重加权因子 wb, 条件概率 p_c)
function sample_site(dmz, α)
    all(isfinite, dmz) || error("direct sampling: non-finite probabilities $dmz")
    psum=sum(dmz)
    psum > 0 || error("direct sampling: all probabilities vanish")
    # 相对偏离 1 过大说明张量不再满足正则条件（环境已归一化，psum 应≈1）
    abs(psum-1) < 1e-8 || @warn "direct sampling: conditional probabilities sum to $psum"
    pb=[p > 0 ? p^α : 0.0 for p in dmz]   # 显式跳过 0，α=0 时也不会抽到禁戒态
    cnow=sample(1:length(pb), Weights(pb))
    pc=dmz[cnow]/psum
    wb=pc/(pb[cnow]/sum(pb))
    return cnow, wb, pc
end

#从左到右直接采样；ARd = dense_slices(右正则张量)
#单条构型的左环境是秩 1 的行向量 L，p_c = ‖L A_c‖²（Σ_c A_c A_c† = I），结果天然非负
function ldirect_sample(ARd, config0, α)
    config1=zeros(Int, 2*l)
    updt=0
    wt=1.0
    L=ones(ComplexF64, 1, 1)
    for il=1:l
        LA=[L*A for A in ARd[il]]
        cnow, wb, pc = sample_site([sum(abs2, v) for v in LA], α)
        wt*=wb
        config1[2*il-1], config1[2*il] = p2s[cnow]
        updt+=(config1[2*il-1]!=config0[2*il-1] || config1[2*il]!=config0[2*il])
        L=LA[cnow]/sqrt(pc)   # 保持 ‖L‖=1，避免长链上概率连乘下溢
    end
    return config1, wt, updt, l
end

#从右到左直接采样；ALd = dense_slices(左正则张量)，p_c = ‖A_c R‖²（Σ_c A_c† A_c = I）
function rdirect_sample(ALd, config0, α)
    config1=zeros(Int, 2*l)
    updt=0
    wt=1.0
    R=ones(ComplexF64, 1, 1)
    for il=l:-1:1
        AR=[A*R for A in ALd[il]]
        cnow, wb, pc = sample_site([sum(abs2, v) for v in AR], α)
        wt*=wb
        config1[2*il-1], config1[2*il] = p2s[cnow]
        updt+=(config1[2*il-1]!=config0[2*il-1] || config1[2*il]!=config0[2*il])
        R=AR[cnow]/sqrt(pc)
    end
    return config1, wt, updt, l
end

# ---------------------------------------------------------------------------
# 电荷密度波结构因子 S(q) = (1/L)(⟨ρ_q†ρ_q⟩ - |⟨ρ_q⟩|²)，ρ_q = Σ_k c†_k c_{k+q}，L = 2l
# q = 2π·qidx/L；CDW 位于 q = π，即 qidx = l
# 动量模式编号与 ms_Hamiltonian 相同，算符串行格式 [系数, 算符…]：正数为产生算符，负数为湮灭算符
# ---------------------------------------------------------------------------

#ρ_q 与 ρ_q†ρ_q 的算符串，ρ_q†ρ_q = Σ_{k,k'} c†_{k+q} c_k c†_{k'} c_{k'+q}
function cdw_ostrings(qidx)
    L = 2*l
    k2p = [x <= l ? 2x - 1 : 2 * (2l - x + 1) for x in 1:L]
    shift(x) = mod(x + qidx - 1, L) + 1
    rho = zeros(Float64, L, 3)
    for ik = 1:L
        rho[ik, :] = [1.0 k2p[ik] -k2p[shift(ik)]]
    end
    rr = zeros(Float64, L^2, 5)
    for ik = 1:L, ik2 = 1:L
        rr[(ik-1)*L+ik2, :] = [1.0 k2p[shift(ik)] -k2p[ik] k2p[ik2] -k2p[shift(ik2)]]
    end
    return rho, rr
end

#Dm 为 mc_sample 返回的采样字典（权重 Dm[cf][1] 未归一化），mt 为同一时刻的 MPS 张量
function cdw_structure_factor(Dm, mt, qidx)
    As = dense_slices(mt)
    rho, rr = cdw_ostrings(qidx)
    Dall = Dict{Vector{Int},ComplexF64}()
    wsum = 0.0
    rho_avg = 0.0im
    rr_avg = 0.0im
    for (cf, v) in Dm
        w = v[1]
        (ro, _) = observable_config(rho, As, v[3], cf, Dall)
        (rro, _) = observable_config(rr, As, v[3], cf, Dall)
        wsum += w
        rho_avg += w * ro
        rr_avg += w * rro
    end
    rho_avg /= wsum
    rr_avg /= wsum
    return real(rr_avg - abs2(rho_avg)) / (2*l)
end

# ---------------------------------------------------------------------------
# 磁化（层赝自旋）结构因子 S(L) = (1/N_e) Σ_M ⟨V_LM† V_LM⟩，N_e = l（半满）
# 费米子模式：轨道 o = m+s+1 的两个层 f=1,2 分别为模式 2o-1, 2o；σ_1=+1, σ_2=-1
# V_LM = Σ_{m,f} (-1)^{s-m} (s L s; -m M m-M) σ_f c†_{m f} c_{m-M, f}（Wigner–Eckart，无额外归一化）
# 对旋转对称态 S(L) 与 M 无关；无关联（无序）时 S(L) 与 L 无关，故 R = 1 - S(1)/S(0) 对无序态趋于 0
# ---------------------------------------------------------------------------

mag_mode(m, f) = 2*(Int(round(m+s))+1) - (f==1 ? 1 : 0)   #轨道 m、层 f 对应的模式编号（1-based）

#V_LM 的单体项列表 [系数, 模式…]：返回 (m, f, w) 三元组，w = (-1)^{s-m}(s L s; -m M m-M) σ_f
function mag_terms(L, M)
    terms = Tuple{Float64,Int,Float64}[]   #(m, f, w) 的 w 不含 σ
    for im = 1:l
        m = im - 1 - s
        m2 = m - M
        (-s <= m2 <= s) || continue
        w = Float64((-1)^Int(round(s-m)) * wigner3j(s, L, s, -m, M, m2))   #wigner3j 返回精确有理根，需转为浮点
        for f = 1:2
            push!(terms, (m, f, w))
        end
    end
    return terms
end

#V_LM† V_LM 的算符串（4 个算符，按算符乘积顺序排列，与 observable_config 约定一致）
function mag_ostrings(L, M)
    terms = mag_terms(L, M)
    sig(f) = f == 1 ? 1.0 : -1.0
    ops = zeros(Float64, length(terms)^2, 5)
    k = 0
    for (ma, fa, wa) in terms, (mb, fb, wb) in terms
        k += 1
        #V† = Σ_a conj(w_a) σ_a c†_{m_a-M} c_{m_a};  V = Σ_b w_b σ_b c†_{m_b} c_{m_b-M}
        ops[k, :] = [wa*wb*sig(fa)*sig(fb), mag_mode(ma-M, fa), -mag_mode(ma, fa), mag_mode(mb, fb), -mag_mode(mb-M, fb)]
    end
    return ops
end

#单体算符 V_LM 的算符串，用于检验（作用在 Fock 态上的系数）
function mag_onebody_ostrings(L, M)
    terms = mag_terms(L, M)
    sig(f) = f == 1 ? 1.0 : -1.0
    ops = zeros(Float64, length(terms), 3)
    for (k, (m, f, w)) in enumerate(terms)
        ops[k, :] = [w*sig(f), mag_mode(m, f), -mag_mode(m-M, f)]
    end
    return ops
end

#Dm 为 mc_sample 返回的采样字典，mt 为同一时刻的 MPS 张量
function mag_structure_factor(Dm, mt, L)
    As = dense_slices(mt)
    wsum = sum(v[1] for v in values(Dm))
    S = 0.0
    for M = -L:L
        ops = mag_ostrings(L, M)
        Dall = Dict{Vector{Int},ComplexF64}()
        vv_avg = 0.0im
        for (cf, v) in Dm
            (vo, _) = observable_config(ops, As, v[3], cf, Dall)
            vv_avg += v[1] * vo
        end
        S += real(vv_avg / wsum)
    end
    return S / l
end
