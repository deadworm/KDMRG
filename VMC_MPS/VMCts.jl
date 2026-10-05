#VMC for t-V model
using Pkg
Pkg.activate(".")
Pkg.instantiate()
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
bond_dim = 100
dm = 0.01
α = 0.5
bsteps = 20
Δτ=0.01
Δb = 10
gsteps = 200
mcsample = 10000

enstring = ms_Hamiltonian();

entot=[]
zb=Vect[FermionParity⊠U1Irrep]((0, 0)=>1)
nb=Vect[FermionParity⊠U1Irrep]((n%2, n)=>1)
# config = [isodd(i) ? 0 : 1 for i in 1:(2*n)]
config = [i>n ? 0 : 1 for i in 1:(2*n)]
m = prod_mps(config)

for ig=1:gsteps
    println("sweep: ", ig)
    (wm, em, gm, Dm, mt, hb)=mc_sample(m, config, mcsample, α)
    if mod(ig, bsteps)==1
        # 一步虚时演化 + 扩维：每个键增加 Δb，但不超过 bond_dim
        global m = residual_expand(mt, hb, em, Δτ, Δb, bond_dim)
        @show [dim(space(m[i], 3)) for i in 1:(l-1)]
    else
        global m = normalize!(sr_update(mt, wm, em, gm, Dm, dm))
        push!(entot, em)
    end
end

@show [dim(space(m[i], 3)) for i in 1:(l-1)]
plot(real(entot), title="Energy measurement", xlabel="SR steps", ylabel="Energy", label="VMC")
plot!(zeros(length(entot)) .- 8.0815, label="Target energy")
entropy_list = mps_entropy(m, l)
plot(entropy_list, title="Entanglement entropy", xlabel="Site index", ylabel="Entropy")