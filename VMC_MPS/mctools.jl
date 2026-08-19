#For the mc sampling
function mc_sample(m, config0, mcwarmup, mcsample)
    enlist=Vector{ComplexF64}(undef, mcsample)
    glist=Vector{TensorMap}(undef, mcsample * l)
    englist=Vector{TensorMap}(undef, mcsample * l)

    co0=get_coef(m, config0)
    w0=abs(co0)^2
    accp=0
    totl=0

    for imc=1:(mcwarmup+mcsample)
        for il=1:l
            config1=copy(config0)
            if config1[il]!=config1[mod(il, l)+1]
                config1[il]=1-config1[il]
                config1[mod(il, l)+1]=1-config1[mod(il, l)+1]
                co1=get_coef(m, config1)
                w1=abs(co1)^2
                if rand()<w1/w0
                    config0=config1
                    co0=co1
                    w0=w1
                    accp+=1
                    totl+=1
                else
                    totl+=1
                end
            else
                continue
            end
        end

        #measure observables here
        if imc>mcwarmup
            mcs=imc-mcwarmup
            enm=measure_en(l, m, co0, config0)
            enlist[mcs]=enm
            gradm=measure_grad(m, config0)/conj(co0)
            for il=1:l
                glist[(mcs-1)*l+il] = gradm[il]
                englist[(mcs-1)*l+il] = enm * gradm[il]
            end
        end
    end

    println("average energy: ", mean(enlist))
    println("acceptance ratio: ", accp/totl)
    return enlist, glist, englist
end

#For energy measurement
function measure_en(l, m, co, config)
    enstring=zeros(Float64, l, 5)
    for i=1:l
        enstring[i, :]=[-2 * t * cos(2*π/l*(i-1)) i -i 0 0]
        # enstring[(i-1)*2+1, :]=[-t mod(i, l)+1 -i 0 0]
        # enstring[(i-1)*2+2, :]=[-t i -mod(i, l)-1 0 0]
    end
    enconf=observable_config(enstring, m, co, config)
    return enconf
end

#For energy gradient measurement
function measure_grad(m, config)
    mf = prod_mps(config)
    lt = Vector{TensorMap}(undef, l + 1)
    rt = Vector{TensorMap}(undef, l + 1)
    ltrt!(lt, rt, m, mf, l)

    mg=Vector{TensorMap}(undef, l)
    for il=1:l
        mg[il]=con_lrt(lt[il], rt[il+1], mf[il])
    end 
    return mg
end

#For mps update
function mps_update(m, enlist, glist, englist, dm)
    nm=Vector{TensorMap}(undef, l)
    lmc=length(enlist)
    men=mean(enlist)
    mg=Vector{TensorMap}(undef, l)
    meng=Vector{TensorMap}(undef, l)
    for il=1:l
        mg[il]=sum(glist[il:l:((lmc-1)*l+il)])/lmc
        meng[il]=sum(englist[il:l:((lmc-1)*l+il)])/lmc
    end
    for il=1:l
        δm = meng[il] - men * mg[il]
        @show il, norm(real(δm))
        nm[il] = m[il] - dm * δm
    end
    return FiniteMPS([nm[il] for il in 1:l])
end