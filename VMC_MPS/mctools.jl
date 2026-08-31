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
    Dm = Dict{Vector{Int},Tuple{Float64,ComplexF64,ComplexF64,Vector{TensorMap},Float64}}()

    updt=0
    totl=0
    for imc=1:mcsample
        if imc%2==1
            (config0, wt, updt1, totl1)=ldirect_sample(m, config0, α)
            updt+=updt1
            totl+=totl1
        else
            (config0, wt, updt1, totl1)=rdirect_sample(m, config0, α)
            updt+=updt1
            totl+=totl1
        end

        #measure observables here
        if !haskey(Dm, config0)
            mf=prod_mps(config0)
            co0=MPSKit.dot(mf, m)
            enm=observable_config(enstring, m, co0, config0, Dall)
            gradm=measure_grad(m, mf, co0)
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
    gm = Vector{TensorMap}(undef, l)
    for il=1:l
        gm[il] = zero(m[il])
    end
    for cf in keys(Dm)
        for il=1:l
            gm[il] += Dm[cf][1] * Dm[cf][4][il] / wm
        end
    end

    println("average energy: ", em)
    println("flip ratio: ", updt/totl)
    println("effective sample ratio: ", wm^2/wm2)
    return wm, em, gm, Dm
end

#For energy gradient measurement
function measure_grad(m, mf, co0)
    lt = Vector{TensorMap}(undef, l + 1)
    rt = Vector{TensorMap}(undef, l + 1)
    ltrt!(lt, rt, m, mf, l)

    mg=Vector{TensorMap}(undef, l)
    for il=1:l
        mg[il]=con_lrt(lt[il], rt[il+1], mf[il])/conj(co0)
    end
    return mg
end

#For configuration direct sampling
function ldirect_sample(m, config0, α)
    config1=zeros(Int, 2*l)
    mz=zeros(Float64, phySpace ← phySpace)
    updt=0
    totl=0
    wt=1
    mlast=0
    # direct sampling from left to right
    ileft=isometry(zb ← zb)
    for il=1:l
        if il!=1
            @planar ileft[r1; r2] := ileft[l1; l2] * mlast'[r1; l1 u] * mlast[l2 u; r2]
        end
        @planar mz[u2; u1] := ileft[l1; l2] * m.AR[il]'[r; l1 u1] * m.AR[il][l2 u2; r]
        dmz=real(vcat((diag(b) for (_, b) in blocks(mz))...))
        @assert all(isfinite, dmz)
        @assert minimum(dmz) >= -1e-12
        @assert sum(dmz) > 0
        pb = dmz .^ α
        @assert sum(pb) > 0
        cnow = sample(1:length(pb), Weights(pb))
        wb = (dmz[cnow]/sum(dmz))/(pb[cnow]/sum(pb))
        wt *= wb
        config1[2*il-1], config1[2*il] = p2s[cnow]
        if config1[2*il-1]==config0[2*il-1]&&config1[2*il]==config0[2*il]
            totl+=1
        else
            updt+=1
            totl+=1
        end
        ip=zeros(phylist[cnow] ← phySpace)
        if cnow == 3
            ip.data[2] = 1.0
        else
            ip.data[1] = 1.0
        end
        mlast = zeros(ComplexF64, codomain(m.AR[il], 1) ⊗ phylist[cnow] ← domain(m.AR[il]))
        @planar mlast[l, u; r] = ip[u; u0] * m.AR[il][l u0; r]
    end
    return config1, wt, updt, totl
end

function rdirect_sample(m, config0, α)
    config1=zeros(Int, 2*l)
    mz=zeros(Float64, phySpace ← phySpace)
    updt=0
    totl=0
    wt=1
    mlast=0
    # direct sampling from right to left
    iright=isometry(nb ← nb)
    for il=l:-1:1
        if il!=l
            @planar iright[l2; l1] := mlast'[r1; l1 u] * mlast[l2 u; r2] * iright[r2; r1]
        end
        @planar mz[u2; u1] := m.AL[il]'[r1; l u1] * m.AL[il][l u2; r2] * iright[r2; r1]
        dmz=real(vcat((diag(b) for (_, b) in blocks(mz))...))
        @assert all(isfinite, dmz)
        @assert minimum(dmz) >= -1e-12
        @assert sum(dmz) > 0
        pb = dmz .^ α
        @assert sum(pb) > 0
        cnow = sample(1:length(pb), Weights(pb))
        wb = (dmz[cnow]/sum(dmz))/(pb[cnow]/sum(pb))
        wt *= wb

        config1[2*il-1], config1[2*il] = p2s[cnow]
        if config1[2*il-1]==config0[2*il-1]&&config1[2*il]==config0[2*il]
            totl+=1
        else
            updt+=1
            totl+=1
        end
        ip=zeros(phylist[cnow] ← phySpace)
        if cnow == 3
            ip.data[2] = 1.0
        else
            ip.data[1] = 1.0
        end
        mlast = zeros(ComplexF64, codomain(m.AL[il], 1) ⊗ phylist[cnow] ← domain(m.AL[il]))
        @planar mlast[l, u; r] = ip[u; u0] * m.AL[il][l u0; r]
    end
    return config1, wt, updt, totl
end

