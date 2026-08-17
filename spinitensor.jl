#use itensor code for benchmark
using ITensors, ITensorMPS
@time begin
    let
        lx = 1
        ly = 4
        fm = 1
        afm = 1
        N = lx * ly
        sites = siteinds("S=1/2", N; conserve_qns=false)

        os = OpSum()
        os += -fm, "Sx", 1, "Sx", 2
        os += -fm, "Sy", 1, "Sy", 2
        os += -fm, "Sz", 1, "Sz", 2
        os += -fm, "Sx", 3, "Sx", 4
        os += -fm, "Sy", 3, "Sy", 4
        os += -fm, "Sz", 3, "Sz", 4
        os += afm, "Sx", 2, "Sx", 3
        os += afm, "Sy", 2, "Sy", 3
        os += afm, "Sz", 2, "Sz", 3

        H = MPO(os, sites, splitblocks=false)

        psi0 = random_mps(sites, linkdims=1)

        nsweeps = 10
        maxdim = [800]
        cutoff = [1e-14]

        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

        return
    end
end
