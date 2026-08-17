#For initial mps
function initial_mps(l, n, phySpace, bond_dim)
    bond_list = []
    push!(bond_list, Vect[FermionParity⊠U1Irrep]((0, 0) => 1))
    for i = 1:(l-1)
        push!(bond_list, fuse(bond_list[end], phySpace))
    end
    push!(bond_list, Vect[FermionParity⊠U1Irrep]((0, n) => 1))

    m = FiniteMPS([randn(Float64, bond_list[i] ⊗ phySpace ← bond_list[i+1]) for i in 1:l])
    @show [dim(space(m[i], 3)) for i in 1:l]

    m = changebonds(m, SvdCut(trscheme=truncdim(bond_dim)))
    @show [dim(space(m[i], 3)) for i in 1:l]

    return m
end

#For the configuration coefficient of mps
function get_config(m, config)
    bond_list = []
    push!(bond_list, Vect[FermionParity⊠U1Irrep]((0, 0) => 1))
    for i = 1:l
        push!(bond_list, fuse(bond_list[end], Vect[FermionParity⊠U1Irrep]((config[i], config[i]) => 1)))
    end
    mp = FiniteMPS([ones(Float64, bond_list[i] ⊗ phySpace ← bond_list[i+1]) for i in 1:l])
    cm = MPSKit.dot(mp, m)

    return cm
end

#For computing observables
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
            osum+=oc[it]*get_config(m, nconfig)
        end
    end
    olocal=osum/cm

    return olocal
end