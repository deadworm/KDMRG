#For initial random mps
function rnd_mps(l, n, bond_dim)
    bond_list = []
    push!(bond_list, zb)
    for i = 1:(l-1)
        push!(bond_list, fuse(bond_list[end], phySpace))
    end
    push!(bond_list, nb)

    m = FiniteMPS([randn(ComplexF64, bond_list[i] ⊗ phySpace ← bond_list[i+1]) for i in 1:l])
    @show [dim(space(m[i], 3)) for i in 1:l]

    m = changebonds(m, SvdCut(trscheme=truncdim(bond_dim)))
    @show [dim(space(m[i], 3)) for i in 1:l]

    return m
end

#For initial direct-product mps
function prod_mps(config)
    bond_list = []
    push!(bond_list, zb)
    for i = 1:l
        push!(bond_list, fuse(bond_list[end], Vect[FermionParity ⊠ U1Irrep]((config[i], config[i]) => 1)))
    end
    mf = FiniteMPS([ones(Float64, bond_list[i] ⊗ phySpace ← bond_list[i+1]) for i in 1:l])
    return mf
end

#For the configuration coefficient of mps
function get_coef(m, config)
    mf = prod_mps(config)
    cm = MPSKit.dot(mf, m)
    return cm
end

#For computing observables [coefficient, site1, site2, ...] for each term
function observable_config(ostring, m, cm, config)
    oterms=size(ostring, 1)
    olength=size(ostring, 2)
    oc=ostring[:, 1]

    nconfig=[]
    osum=0
    for it=1:oterms
        nconfig=copy(config)
        for il=0:(olength-2)
            oind=Int(ostring[it, end-il])
            if oind==0
                continue
            elseif oind>0&&nconfig[oind]==0
                oc[it]=oc[it]*(-1)^sum(nconfig[1:(oind-1)])
                nconfig[oind]=1
            elseif oind<0&&nconfig[-oind]==1
                oc[it]=oc[it]*(-1)^sum(nconfig[1:(-oind-1)])
                nconfig[-oind]=0
            else
                oc[it]=0
                break
            end
        end
        if oc[it]!=0
            osum+=oc[it] * get_coef(m, nconfig)
        end
    end
    olocal=osum / cm

    return olocal
end