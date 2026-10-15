#VMC for t-V model
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using TensorKit, MPSKit
using LinearAlgebra, Random, StatsBase, IterativeSolvers
using WignerSymbols
using CUDA
using Plots

include("mpstools.jl")
include("mctools.jl")
s2p = Dict((0, 0)=>1, (1, 0)=>2, (0, 1)=>3, (1, 1)=>4)
p2s = [(0, 0), (1, 0), (0, 1), (1, 1)]

l = 10
s = (l-1)/2
n = l
V0=4.75
V1=1
U1=(4s-1)*V1/s/(2s+1)^2
U0=(4s+1)*V0/(2s+1)^2+U1*s
hlist=3:0.1:4

Rblist=zeros(length(hlist), 1)
Rerlist=zeros(length(hlist), 1)
for idh=1:length(hlist)
    h=hlist[idh]
    bond_dim = 200
    dm = 0.01
    α = 0.5
    bsteps = 20
    Δτ = 0.001
    Δb = 10
    gsteps = 400
    mcsample = 10000

    enstring = fsIsing_Hamiltonian()

    entot=[]
    config = [isodd(i) ? 0 : 1 for i in 1:(2*n)]
    #目标 Lz 扇区由初始构型决定
    Lz = sum(orbital_m(il) * (config[2*il-1] + config[2*il]) for il in 1:l)
    zb=Vect[SymSector]((0, 0, 0)=>1)
    nb=Vect[SymSector]((n%2, n, Lz)=>1)
    m = prod_mps(config)

    for ig=1:gsteps
        println("sweep: ", ig)
        (wm, em, gm, Dm, mt, hb)=mc_sample(m, config, mcsample, α)
        if mod(ig, bsteps)==1
            # 一步虚时演化 + 扩维：每个键增加 Δb，但不超过 bond_dim
            global m = residual_expand(m, hb, em, Δτ, Δb, bond_dim)
            @show [dim(space(m[i], 3)) for i in 1:(l-1)]
        else
            global m = normalize!(sr_update(mt, wm, em, gm, Dm, dm))
            push!(entot, em)
        end
    end

    @show [dim(space(m[i], 3)) for i in 1:(l-1)]
    plot(real(entot), title="Energy measurement", xlabel="SR steps", ylabel="Energy", label="VMC")
    entropy_list = mps_entropy(m, l)
    plot(entropy_list, title="Entanglement entropy", xlabel="Site index", ylabel="Entropy")

    #能量与磁化结构因子：多批独立采样，取批间平均
    Rbh=[]
    Nbins=100
    for irbs=1:Nbins
        (wm, em, gm, Dm, mt, hb)=mc_sample(m, config, mcsample, α)
        global m = normalize!(sr_update(mt, wm, em, gm, Dm, dm))
        S0_b=mag_structure_factor(Dm, mt, 0)
        S1_b=mag_structure_factor(Dm, mt, 1)
        push!(Rbh, 1 - S1_b / S0_b)
    end
    Rblist[idh]=mean(Rbh)
    Rerlist[idh]=std(Rbh)/sqrt(Nbins)
end

# R10=Rblist
# Re10=Rerlist
# open("R8.txt", "w") do io
#     println(io, R8)
# end
# open("Re8.txt", "w") do io
#     println(io, Re8)
# end
# open("R10.txt", "w") do io
#     println(io, R10)
# end
# open("Re10.txt", "w") do io
#     println(io, Re10)
# end

entropy_list = mps_entropy(m, l)
plot(entropy_list, title="Entanglement entropy", xlabel="Site index", ylabel="Entropy")

plot!(hlist, R8, yerror=Re8, title="Correlation Ratio", xlabel="h", ylabel="R", label="N=8")
plot!(hlist, R10, yerror=Re10, title="Correlation Ratio", xlabel="h", ylabel="R", label="N=10")