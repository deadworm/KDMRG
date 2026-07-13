#dmrg test for Heisenberg model
using TensorKit
using KrylovKit
using Random
using Combinatorics

include("mpstools.jl")
include("pepstools.jl")

lx = 6
ly = 6
l = lx * ly

phySpace = Rep[U₁]((1 / 2) => 1, (-1 / 2) => 1)
config_sz = []
for ix = 1:lx
    for iy = 1:ly
        if (ix + iy) % 2 == 0
            push!(config_sz, Rep[U₁]((1 / 2) => 1))
        else
            push!(config_sz, Rep[U₁]((-1 / 2) => 1))
        end
    end
end

bond_dim = 4
bond_cut = 10
dt = 0.01
nsweep = 40

#generate qn randomly for PEPS
peps = PEPSinit(lx, ly, config_sz, phySpace)
rdPEPS!(lx, ly, peps, bond_dim, dt, nsweep)

#pick physical indices of PEPS
# config_sz[1]=config_sz[2]
# config_sz[2]=config_sz[3]
fpeps = Vector{TensorMap}(undef, l)
for ix = 1:lx
    for iy = 1:ly
        i = (ix - 1) * ly + iy
        ft = TensorMap(ones, config_sz[i] ← phySpace)
        @tensor fpeps[i][p l d; u r] := peps[i][p0 l d; u r] * ft[p; p0]
        # @show fpeps[i].data==peps[i].data
    end
end

#test PEPS contraction
# @time "PEPS contraction test" begin
mps_d = Vector{TensorMap}(undef, l)
mps_u = Vector{TensorMap}(undef, l)
for ix = 1:lx
    mps_d[ix] = fpeps[(ix-1)*ly+1]
    mps_u[(ly-1)*lx+ix] = permute(fpeps[(ix-1)*ly+ly], (1, 2, 4), (3, 5))
end
for iy = 2:ly
    peps_d = Vector{TensorMap}(undef, lx)
    peps_u = Vector{TensorMap}(undef, lx)
    for ix = 1:lx
        peps_d[ix] = fpeps[(ix-1)*ly+iy]
        peps_u[ix] = permute(fpeps[(ix-1)*ly+ly-iy+1], (1, 2, 4), (3, 5))
    end
    mps_d[(iy-1)*lx+1:(iy-1)*lx+lx] = SRCcontract(mps_d[(iy-2)*lx+1:(iy-2)*lx+lx], peps_d, lx, bond_cut)
    mps_u[(ly-iy)*lx+1:(ly-iy)*lx+lx] = SRCcontract(mps_u[(ly-iy+1)*lx+1:(ly-iy+1)*lx+lx], peps_u, lx, bond_cut)
end
for iy = 1:ly
    for ix = 1:lx
        mps_u[(ly-iy)*lx+ix] = permute(mps_u[(ly-iy)*lx+ix], (1, 2, 4), (3, 5))
    end
end
# end

ppeps = Vector{TensorMap}(undef, l)
lt = Vector{TensorMap}(undef, lx)
rt = Vector{TensorMap}(undef, lx)
for ix = 1:lx
    for iy = 1:ly
        if iy == 1
            PEPSlc!(lt, rt, mps_u[lx+1:lx+lx], 0, fpeps, iy, lx)
            if ix == 1
                pPEPSlc!(ppeps, 0, rt, mps_u[lx+1:lx+lx], 0, ix, iy, lx)
            elseif ix == lx
                pPEPSlc!(ppeps, lt, 0, mps_u[lx+1:lx+lx], 0, ix, iy, lx)
            else
                pPEPSlc!(ppeps, lt, rt, mps_u[lx+1:lx+lx], 0, ix, iy, lx)
            end
        elseif iy == ly
            PEPSlc!(lt, rt, 0, mps_d[(ly-2)*lx+1:(ly-2)*lx+lx], fpeps, iy, lx)
            if ix == 1
                pPEPSlc!(ppeps, 0, rt, 0, mps_d[(ly-2)*lx+1:(ly-2)*lx+lx], ix, iy, lx)
            elseif ix == lx
                pPEPSlc!(ppeps, lt, 0, 0, mps_d[(ly-2)*lx+1:(ly-2)*lx+lx], ix, iy, lx)
            else
                pPEPSlc!(ppeps, lt, rt, 0, mps_d[(ly-2)*lx+1:(ly-2)*lx+lx], ix, iy, lx)
            end
        else
            PEPSlc!(lt, rt, mps_u[(iy)*lx+1:(iy)*lx+lx], mps_d[(iy-2)*lx+1:(iy-2)*lx+lx], fpeps, iy, lx)
            if ix == 1
                pPEPSlc!(ppeps, 0, rt, mps_u[(iy)*lx+1:(iy)*lx+lx], mps_d[(iy-2)*lx+1:(iy-2)*lx+lx], ix, iy, lx)
            elseif ix == lx
                pPEPSlc!(ppeps, lt, 0, mps_u[(iy)*lx+1:(iy)*lx+lx], mps_d[(iy-2)*lx+1:(iy-2)*lx+lx], ix, iy, lx)
            else
                pPEPSlc!(ppeps, lt, rt, mps_u[(iy)*lx+1:(iy)*lx+lx], mps_d[(iy-2)*lx+1:(iy-2)*lx+lx], ix, iy, lx)
            end
        end
    end
end
norm(lt[lx])
norm(rt[1])
for i = 1:l
    println("Site $(i) fpepsnorm: ", norm(fpeps[i]))
    println("Site $(i) ppepsnorm: ", norm(ppeps[i]))
    # println("Site $(i) max: ", maximum(abs.(fpeps[i].data)))
end



