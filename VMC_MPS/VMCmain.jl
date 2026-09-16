#VMC for t-V model
using Pkg;
Pkg.activate(".")
using TensorKit, MPSKit
using LinearAlgebra, Random, StatsBase, IterativeSolvers
using CUDA
using Plots

include("contract.jl")
include("mpstools.jl")
include("mctools.jl")
phySpace = Vect[FermionParity⊠U1Irrep]((0,0)=>1,(1,1)=>2,(0,2)=>1)
s2p = Dict((0,0)=>1,(1,0)=>2,(0,1)=>3,(1,1)=>4)
p2s = [(0,0),(1,0),(0,1),(1,1)]
phylist = [Vect[FermionParity⊠U1Irrep]((0,0)=>1),Vect[FermionParity⊠U1Irrep]((1,1)=>1),Vect[FermionParity⊠U1Irrep]((1,1)=>1),Vect[FermionParity⊠U1Irrep]((0,2)=>1)]

srflag=true;
t = 1
V = 1
l = 8
n = 6
bond_dim = 20
zb=Vect[FermionParity⊠U1Irrep]((0,0)=>1)
nb=Vect[FermionParity⊠U1Irrep]((n%2,n)=>1)
m = rnd_mps(l, bond_dim);
enstring = ms_Hamiltonian();

α = 0.5
dm = 0.01
gstep = 200
mcsample = 1000

entot=[]
for ig=1:gstep
    println("sweep: ", ig)
    config = shuffle([i>n ? 0 : 1 for i in 1:(2*l)])
    if srflag
        (wm, em, gm, Dm)=mc_sample(m, config, mcsample, α)
        m=normalize!(sr_update(m, wm, em, gm, Dm, dm))
    else
        (wm, em, gm, Dm)=mc_sample(m, config, mcsample, α)
        m=normalize!(gn_update(m, wm, em, gm, Dm, dm))
    end
    push!(entot, em)
end

plot(real(entot), title="Energy measurement", xlabel="MC bin", ylabel="Energy")
entropy_list = mps_entropy(m, l)
plot(entropy_list, title="Entanglement entropy", xlabel="Site index", ylabel="Entropy")

#Test code for energy measurement
using MPSKitModels
lat = FiniteChain(2*l)

hopping = -t * (MPSKitModels.c_plusmin() + MPSKitModels.c_plusmin()');
hubbard = V * MPSKitModels.c_number() ⊗ MPSKitModels.c_number();
H = @mpoham sum(hopping{i,j} + hubbard{i,j} for (i, j) in nearest_neighbours(lat));
H += @mpoham (-hopping{lat[2*l],lat[1]} + hubbard{lat[2*l],lat[1]});

χ = bond_dim÷2
Vp = physicalspace(H)[1]
V0 = Vect[FermionParity](0=>χ, 1=>χ)
ψ = FiniteMPS(2*l, Vp, V0);
ψ, _ = find_groundstate(ψ, H);
@show [dim(space(ψ[i], 3)) for i in 1:2*l]

entropy_list = mps_entropy(ψ, 2*l)
plot(entropy_list, title="Entanglement entropy", xlabel="Site index", ylabel="Entropy")
m = changebonds(ψ, SvdCut(trscheme=truncdim(bond_dim)))
