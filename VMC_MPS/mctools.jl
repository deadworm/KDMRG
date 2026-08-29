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
    enstring=zeros(Float64, ls^3 + ls, 5)
    for i1=1:ls
        if i1%2 == 1
            ik1 = i1÷2+1
        else
            ik1 = mod(-i1÷2, ls)+1
        end
        for i2=1:ls
            if i2%2 == 1
                ik2 = i2÷2+1
            else
                ik2 = mod(-i2÷2, ls)+1
            end
            for i=1:ls
                if i%2 == 1
                    iq = i÷2+1
                else
                    iq = mod(-i÷2, ls)+1
                end
                enstring[(ik1-1)*ls^2+(ik2-1)*ls+iq, :]=[V/ls*cos(2*π/ls*(iq)) ik1 -mod(ik1+iq-1, ls)-1 ik2 -mod(ik2-iq-1, ls)-1]
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
function mc_sample(m, config0, mcsample)
    Dall = Dict{Vector{Int},ComplexF64}()
    Dm = Dict{Vector{Int},Tuple{Int,ComplexF64,ComplexF64,Vector{TensorMap}}}()

    updt=0
    totl=0
    for imc=1:mcsample
        if imc%2==1
            (config0, updt1, totl1)=ldirect_sample(m, config0)
            updt+=updt1
            totl+=totl1
        else
            (config0, updt1, totl1)=rdirect_sample(m, config0)
            updt+=updt1
            totl+=totl1
        end

        #measure observables here
        if !haskey(Dm, config0)
            mf=prod_mps(config0)
            co0=MPSKit.dot(mf, m)
            enm=observable_config(enstring, m, co0, config0, Dall)
            gradm=measure_grad(m, mf)/conj(co0)
            Dm[config0] = (1, co0, enm, gradm)
        else
            Dm[config0] = (Dm[config0][1]+1, Dm[config0][2], Dm[config0][3], Dm[config0][4])
        end
    end

    em = 0
    for cf in keys(Dm)
        em += Dm[cf][1] * Dm[cf][3] / mcsample
    end
    gm = Vector{TensorMap}(undef, l)
    for il=1:l
        gm[il] = zero(m[il])
    end
    for cf in keys(Dm)
        for il=1:l
            gm[il] += Dm[cf][1] * Dm[cf][4][il] / mcsample
        end
    end

    println("average energy: ", em)
    println("update ratio: ", updt/totl)
    return em, gm, Dm
end

#For configuration direct sampling
function ldirect_sample(m, config0)
    config1=zeros(Int, 2*l)
    mz=zeros(Float64, phySpace ← phySpace)
    updt=0
    totl=0
    mlast=0
    # direct sampling from left to right
    ileft=isometry(zb ← zb)
    for il=1:l
        if il!=1
            @planar ileft[r1; r2] := ileft[l1; l2] * mlast'[r1; l1 u] * mlast[l2 u; r2]
        end
        @planar mz[u2; u1] := ileft[l1; l2] * m.AR[il]'[r; l1 u1] * m.AR[il][l2 u2; r]
        dmz=vcat((diag(b) for (_, b) in blocks(mz))...)
        cnow = sample(1:length(dmz), Weights(real(dmz))) - 1
        config1[2*il], config1[2*il-1] = divrem(cnow, 2)
        if config1[2*il-1]==config0[2*il-1]&&config1[2*il]==config0[2*il]
            totl+=1
        else
            updt+=1
            totl+=1
        end
        ip=zeros(phylist[cnow+1] ← phySpace)
        if cnow == 2
            ip.data[2] = 1.0
        else
            ip.data[1] = 1.0
        end
        mlast = zeros(ComplexF64, codomain(m.AR[il], 1) ⊗ phylist[cnow+1] ← domain(m.AR[il]))
        @planar mlast[l, u; r] = ip[u; u0] * m.AR[il][l u0; r]
    end
    return config1, updt, totl
end

function rdirect_sample(m, config0)
    config1=zeros(Int, 2*l)
    mz=zeros(Float64, phySpace ← phySpace)
    updt=0
    totl=0
    mlast=0
    # direct sampling from right to left
    iright=isometry(nb ← nb)
    for il=l:-1:1
        if il!=l
            @planar iright[l2; l1] := mlast'[r1; l1 u] * mlast[l2 u; r2] * iright[r2; r1]
        end
        @planar mz[u2; u1] := m.AL[il]'[r1; l u1] * m.AL[il][l u2; r2] * iright[r2; r1]
        dmz=vcat((diag(b) for (_, b) in blocks(mz))...)
        cnow = sample(1:length(dmz), Weights(real(dmz))) - 1
        config1[2*il], config1[2*il-1] = divrem(cnow, 2)
        if config1[2*il-1]==config0[2*il-1]&&config1[2*il]==config0[2*il]
            totl+=1
        else
            updt+=1
            totl+=1
        end
        ip=zeros(phylist[cnow+1] ← phySpace)
        if cnow == 2
            ip.data[2] = 1.0
        else
            ip.data[1] = 1.0
        end
        mlast = zeros(ComplexF64, codomain(m.AL[il], 1) ⊗ phylist[cnow+1] ← domain(m.AL[il]))
        @planar mlast[l, u; r] = ip[u; u0] * m.AL[il][l u0; r]
    end
    return config1, updt, totl
end

#For energy gradient measurement
function measure_grad(m, mf)
    lt = Vector{TensorMap}(undef, l + 1)
    rt = Vector{TensorMap}(undef, l + 1)
    ltrt!(lt, rt, m, mf, l)

    mg=Vector{TensorMap}(undef, l)
    for il=1:l
        mg[il]=con_lrt(lt[il], rt[il+1], mf[il])
    end
    return mg
end

