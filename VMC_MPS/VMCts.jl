#VMC for t-V model
using Pkg
Pkg.activate(".")
using TensorKit, MPSKit
using LinearAlgebra, Random, StatsBase, IterativeSolvers
using CUDA
using Plots

include("mpstools.jl")
include("mctools.jl")
phySpace = Vect[FermionParity⊠U1Irrep]((0, 0)=>1, (1, 1)=>2, (0, 2)=>1)
s2p = Dict((0, 0)=>1, (1, 0)=>2, (0, 1)=>3, (1, 1)=>4)
p2s = [(0, 0), (1, 0), (0, 1), (1, 1)]
phylist = [Vect[FermionParity⊠U1Irrep]((0, 0)=>1), Vect[FermionParity⊠U1Irrep]((1, 1)=>1), Vect[FermionParity⊠U1Irrep]((1, 1)=>1), Vect[FermionParity⊠U1Irrep]((0, 2)=>1)]

t = 1
V = 1
l = 8
n = 8
bond_dim = 20
α = 0.2
dm = 0.01
gstep = 200
mcsample = 1000
enstring = ms_Hamiltonian();

entot=[]
zb=Vect[FermionParity⊠U1Irrep]((0, 0)=>1)
nb=Vect[FermionParity⊠U1Irrep]((n%2, n)=>1)
# config = [isodd(i) ? 0 : 1 for i in 1:(2*n)]
config = [i>n ? 0 : 1 for i in 1:(2*n)]
m = prod_mps(config)
al=ceil(Int, n/2)
ar=al+1

for ib=1:5
    for ig=1:10
        m = expand_active(m, al:ar; add=bond_dim)
    end
    @show [dim(space(m[i], 3)) for i in 1:(l-1)]
    m = changebonds(m, SvdCut(trscheme=truncdim(bond_dim)))

    for ig=1:gstep
        println("sweep: ", ig)
        (wm, em, gm, Dm, mt)=mc_sample(m, config, mcsample, α)
        m=normalize!(sr_update(mt, wm, em, gm, Dm, dm))
        push!(entot, em)
    end

    al=max(al-1, 1)
    ar=min(ar+1, l)
end

@show [dim(space(m[i], 3)) for i in 1:(l-1)]
plot(real(entot), title="Energy measurement", xlabel="SR steps", ylabel="Energy", label="VMC")
plot!(zeros(length(entot)) .- 8.083075391506974, label="Target energy")
entropy_list = mps_entropy(m, l)
plot(entropy_list, title="Entanglement entropy", xlabel="Site index", ylabel="Entropy")

# #Test code for energy measurement
# using MPSKitModels
# lat = FiniteChain(2*l)

# hopping = -t * (MPSKitModels.c_plusmin() + MPSKitModels.c_plusmin()');
# hubbard = V * MPSKitModels.c_number() ⊗ MPSKitModels.c_number();
# H = @mpoham sum(hopping{i,j} + hubbard{i,j} for (i, j) in nearest_neighbours(lat));
# H += @mpoham (-hopping{lat[2*l],lat[1]} + hubbard{lat[2*l],lat[1]});

# χ = 100#bond_dim÷2
# Vp = physicalspace(H)[1]
# V0 = Vect[FermionParity](0=>χ, 1=>χ)
# ψ = FiniteMPS(2*l, Vp, V0);
# ψ, _ = find_groundstate(ψ, H);
# ψ = changebonds(ψ, SvdCut(trscheme=truncerr(1e-1)))
# @show [dim(space(ψ[i], 3)) for i in 1:(2*l)]

# entropy_list = mps_entropy(ψ, 2*l)
# plot(entropy_list, title="Entanglement entropy", xlabel="Site index", ylabel="Entropy")
# m = changebonds(ψ, SvdCut(trscheme=truncdim(bond_dim)))