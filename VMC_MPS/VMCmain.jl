#VMC for t-V model
using Pkg;Pkg.activate(".")
using TensorKit, MPSKit
using Random
using Statistics
using Plots
# using Combinatorics

include("contract.jl")
include("mpstools.jl")
include("mctools.jl")
phySpace = Vect[FermionParity ⊠ U1Irrep]((0, 0) => 1, (1, 1) => 1)

t = 1
l = 4
n = 2
bond_dim = 8
zb=Vect[FermionParity ⊠ U1Irrep]((0, 0) => 1)
nb=Vect[FermionParity ⊠ U1Irrep]((0, n) => 1)
m = rnd_mps(l, n, bond_dim)

dm = 0.01
gstep = 100
mstep = 100
mcwarmup = 100
mcsample = 100

entot=[]
for ig=1:(gstep+mstep)
    c0 = shuffle([i>n ? 0 : 1 for i in 1:l])
    (enlist, glist, englist)=mc_sample(m, c0, mcwarmup, mcsample);
    m = normalize!(mps_update(m, enlist, glist, englist, dm))
    if ig>gstep
        push!(entot, mean(enlist))
    end
end

plot(real(entot), title="Energy measurement", xlabel="MC bin", ylabel="Energy")


#Test code for l=4, n=2
configtest=[0 0 1 1;
    0 1 0 1;
    0 1 1 0;
    1 0 0 1;
    1 0 1 0;
    1 1 0 0];
normtest=0
cm = zeros(ComplexF64, 1, 6)
for id = 1:6
    config = configtest[id, :]
    cm[id] = get_coef(m, config)
    println("config: ", config, " coefficient: ", cm[id])
    normtest += abs(cm[id])^2
end
println("normtest: ", normtest)

ostringtest=[0.5 3 2 -3 -2;
    0.5 2 3 -3 -2;
    0.5 0 -3 1 0;
    0.5 3 0 -3 0;]
olocal=observable_config(ostringtest, m, cm[3], configtest[3, :])

#Test code for energy measurement
using MPSKitModels
lat = FiniteChain(l)

hopping = -t * (MPSKitModels.c_plusmin() + MPSKitModels.c_plusmin()')
H = @mpoham sum(hopping{i,j} for (i, j) in nearest_neighbours(lat))
H += @mpoham (hopping{lat[l],lat[1]});

χ = 8
Vp = physicalspace(H)[1]
V = Vect[FermionParity](0=>χ, 1=>χ)
ψ = FiniteMPS(l, Vp, V)
ψ, _ = find_groundstate(ψ, H);

expectation_value(ψ, H)
m = changebonds(ψ, SvdCut(trscheme=truncdim(bond_dim)))