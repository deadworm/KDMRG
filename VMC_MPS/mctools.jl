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

#For the mc sampling
function mc_sample(m, config0, mcsample, α)
    Dall = Dict{Vector{Int},ComplexF64}()
    Dm = Dict{Vector{Int},Tuple{Float64,ComplexF64,ComplexF64,Vector{Any},Float64}}()
    # 采样会访问 m.AL/m.AR 并移动 gauge center，之后 m[il] 会换成另一套规范，
    # 所以先固定一份张量：梯度和参数更新都必须基于这同一组张量
    mt=[copy(m[il]) for il=1:l]
    md=Vector{Any}(undef, l)
    for il=1:l
        md[il]=convert(Array, mt[il])
    end

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
            co0=mps_slice(md, config0)
            enm=observable_config(enstring, md, co0, config0, Dall)
            gradm=measure_grad(md, co0, config0)
            Dm[config0] = (wt/mcsample, co0, enm, gradm, wt^2/mcsample)
        else
            Dm[config0] = (Dm[config0][1]+wt/mcsample, Dm[config0][2], Dm[config0][3], Dm[config0][4], Dm[config0][5]+wt^2/mcsample)
        end
    end
    wm = 0
    wm2 = 0
    for cf in keys(Dm)
        wm += Dm[cf][1]
        wm2 += Dm[cf][5]
    end
    em = 0
    for cf in keys(Dm)
        em += Dm[cf][1] * Dm[cf][3] / wm
    end
    gm = Vector{Any}(undef, l)
    for il=1:l
        gm[il] = zeros(ComplexF64, size(md[il]))
        for cf in keys(Dm)
            gm[il] += Dm[cf][1] * Dm[cf][4][il] / wm
        end
    end

    println("average energy: ", em)
    println("flip ratio: ", updt/totl)
    println("effective sample ratio: ", wm^2/wm2)
    return wm, em, gm, Dm, mt
end

#For energy gradient measurement
function measure_grad(md, co0, config0)
    vl = Vector{Any}(undef, l + 1)
    vr = Vector{Any}(undef, l + 1)
    vl[1] = ComplexF64[1]
    for il = 2:(l+1)
        vl[il] = md[il-1][:, s2p[config0[(il-2)*2+1], config0[(il-2)*2+2]], :]' * vl[il-1]
    end
    vr[l+1] = ComplexF64[1]
    for il = l:-1:1
        vr[il] = vr[il+1] * md[il][:, s2p[config0[(il-1)*2+1], config0[(il-1)*2+2]], :]'
    end

    mg=Vector{Any}(undef, l)
    for il=1:l
        mb=zeros(ComplexF64, size(md[il]))
        mb[:, s2p[config0[(il-1)*2+1], config0[(il-1)*2+2]], :]=kron(vl[il],vr[il+1])
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
