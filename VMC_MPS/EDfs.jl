#ED for the bilayer Ising model in the LLL (same Hamiltonian as VMCfs.jl)
#扫描 h，计算基态的磁化结构因子 S(0), S(1) 与 correlation ratio R = 1 - S(1)/S(0)
using Pkg
Pkg.activate(".")
using LinearAlgebra, SparseArrays, KrylovKit
using Combinatorics
using WignerSymbols
using Printf

include("mctools.jl")

l = 12
s = (l-1)/2
n = l
V0=4.75
V1=1
U1=(4s-1)*V1/s/(2s+1)^2
U0=(4s+1)*V0/(2s+1)^2+U1*s
hlist = 3:0.1:4

orbital_m(il) = il - 1 - s
config_lz(cf) = sum(orbital_m(il) * (cf[2*il-1] + cf[2*il]) for il in 1:l)

#算符串作用在构型上：O|config0⟩ = Σ Dc[cf] |cf⟩，符号约定与 observable_config 相同
function apply_ostring(ostring, config0)
    Dc = Dict{Vector{Int},Float64}()
    for it=1:size(ostring, 1)
        c = ostring[it, 1]
        nconfig = copy(config0)
        for il=0:(size(ostring, 2)-2)
            oind = Int(ostring[it, end-il])
            if oind==0
                continue
            elseif oind>0 && nconfig[oind]==0
                c *= (-1)^sum(@view nconfig[1:(oind-1)])
                nconfig[oind] = 1
            elseif oind<0 && nconfig[-oind]==1
                c *= (-1)^sum(@view nconfig[1:(-oind-1)])
                nconfig[-oind] = 0
            else
                c = 0.0
                break
            end
        end
        c != 0 && (Dc[nconfig] = get(Dc, nconfig, 0.0) + c)
    end
    return Dc
end

#固定粒子数 n、固定 Lz 的构型基
function lz_basis(lz)
    basis = Vector{Int}[]
    for occ in combinations(1:(2*l), n)
        cf = zeros(Int, 2*l)
        cf[occ] .= 1
        config_lz(cf) == lz && push!(basis, cf)
    end
    return basis
end

function op_matrix(ostring, basis)
    idx = Dict(cf => i for (i, cf) in enumerate(basis))
    I, J, X = Int[], Int[], Float64[]
    for (j, cf) in enumerate(basis)
        for (cf1, c) in apply_ostring(ostring, cf)
            push!(I, idx[cf1])
            push!(J, j)
            push!(X, c)
        end
    end
    return sparse(I, J, X, length(basis), length(basis))
end

#S(L) = (1/l) Σ_M ⟨V_LM† V_LM⟩，与 mag_structure_factor 定义一致
mag_S(ψ, basis, L) = sum(real(dot(ψ, op_matrix(mag_ostrings(L, M), basis) * ψ)) for M=(-L):L) / l

lzs = sort(unique(config_lz(cf) for occ in combinations(1:(2*l), n) for cf in [begin c=zeros(Int, 2*l); c[occ].=1; c end]))
bases = Dict(lz => lz_basis(lz) for lz in lzs)
println("sector dims: ", [(lz, length(bases[lz])) for lz in lzs if lz >= 0])

results = []
for hv in hlist
    global h = hv
    enstring = fsIsing_Hamiltonian()
    #各 Lz 扇区的最低能量，确认基态所在扇区
    e_sec = Dict{Float64,Float64}()
    for lz in lzs
        lz < 0 && continue   #Lz → -Lz 对称
        H = op_matrix(enstring, bases[lz])
        e_sec[lz] = length(bases[lz]) <= 200 ? eigmin(Symmetric(Matrix(H))) :
                    first(eigsolve(H, 1, :SR; ishermitian=true)[1])
    end
    lz0 = argmin(e_sec)
    H = op_matrix(enstring, bases[0.0])
    herr = norm(H - H')
    vals, vecs, _ = eigsolve(H, 3, :SR; ishermitian=true, tol=1e-12)
    ψ = vecs[1]
    S0 = mag_S(ψ, bases[0.0], 0)
    S1 = mag_S(ψ, bases[0.0], 1)
    R = 1 - S1/S0
    push!(results, (hv, vals[1], vals[2]-vals[1], lz0, S0, S1, R))
    @printf("h=%.1f  E0=%.10f  gap=%.3e  GS Lz=%s  |H-H'|=%.1e  S0=%.6f  S1=%.6f  R=%.6f\n",
        hv, vals[1], vals[2]-vals[1], lz0, herr, S0, S1, R)
end

Rlist=zeros(length(hlist), 1)
for ih=1:length(hlist)
    Rlist[ih]=results[ih][7]
end
plot!(hlist, Rlist, title="Correlation Ratio", xlabel="h", ylabel="R", label="ED,N=10"; marker=:circle)

open("ED_R_l$(l).dat", "w") do io
    println(io, "# h E0 gap S0 S1 R   (l=$l, n=$n, V0=$V0, V1=$V1, Lz=0)")
    for (hv, e0, gap, _, S0, S1, R) in results
        @printf(io, "%.2f %.12f %.6e %.10f %.10f %.10f\n", hv, e0, gap, S0, S1, R)
    end
end
