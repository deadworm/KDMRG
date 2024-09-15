#use itensor code for benchmark
using ITensors, ITensorMPS
@time begin
    let
        lx = 1
        ly = 4
        N = lx * ly
        sites = siteinds("Electron", N; conserve_qns=true)

        os = OpSum()
        # for j = 1:N-2
        #     os += "Sz", j, "Sz", j + 2
        #     os += 1 / 2, "S+", j, "S-", j + 2
        #     os += 1 / 2, "S-", j, "S+", j + 2
        # end
        # for j = 1:N-1
        #     os += "Sz", j, "Sz", j + 1
        #     os += 1 / 2, "S+", j, "S-", j + 1
        #     os += 1 / 2, "S-", j, "S+", j + 1
        # end

        # for ix = 1:lx-1
        #     for iy = 1:ly-1
        #         os += "Cdagup", (ix - 1) * ly + iy, "Cup", (ix - 1) * ly + iy + 1
        #         os += "Cdagup", (ix - 1) * ly + iy + 1, "Cup", (ix - 1) * ly + iy
        #         os += "Cdagdn", (ix - 1) * ly + iy, "Cdn", (ix - 1) * ly + iy + 1
        #         os += "Cdagdn", (ix - 1) * ly + iy + 1, "Cdn", (ix - 1) * ly + iy

        #         os += "Cdagup", (ix - 1) * ly + iy, "Cup", (ix) * ly + iy
        #         os += "Cdagup", (ix) * ly + iy, "Cup", (ix - 1) * ly + iy
        #         os += "Cdagdn", (ix - 1) * ly + iy, "Cdn", (ix) * ly + iy
        #         os += "Cdagdn", (ix) * ly + iy, "Cdn", (ix - 1) * ly + iy
        #     end
        #     os += "Cdagup", (ix - 1) * ly + ly, "Cup", (ix - 1) * ly + 1
        #     os += "Cdagup", (ix - 1) * ly + 1, "Cup", (ix - 1) * ly + ly
        #     os += "Cdagdn", (ix - 1) * ly + ly, "Cdn", (ix - 1) * ly + 1
        #     os += "Cdagdn", (ix - 1) * ly + 1, "Cdn", (ix - 1) * ly + ly
        # end
        # os += "Cdagup", (lx - 1) * ly + ly, "Cup", (lx - 1) * ly + 1
        # os += "Cdagup", (lx - 1) * ly + 1, "Cup", (lx - 1) * ly + ly
        # os += "Cdagdn", (lx - 1) * ly + ly, "Cdn", (lx - 1) * ly + 1
        # os += "Cdagdn", (lx - 1) * ly + 1, "Cdn", (lx - 1) * ly + ly
        # for ix = 1:lx-1
        #     os += "Cdagup", (ix - 1) * ly + ly, "Cup", (ix) * ly + ly
        #     os += "Cdagup", (ix) * ly + ly, "Cup", (ix - 1) * ly + ly
        #     os += "Cdagdn", (ix - 1) * ly + ly, "Cdn", (ix) * ly + ly
        #     os += "Cdagdn", (ix) * ly + ly, "Cdn", (ix - 1) * ly + ly
        # end
        for iy = 1:ly-1
            os += "Cdagup", (lx - 1) * ly + iy, "Cup", (lx - 1) * ly + iy + 1
            os += "Cdagup", (lx - 1) * ly + iy + 1, "Cup", (lx - 1) * ly + iy
            os += "Cdagdn", (lx - 1) * ly + iy, "Cdn", (lx - 1) * ly + iy + 1
            os += "Cdagdn", (lx - 1) * ly + iy + 1, "Cdn", (lx - 1) * ly + iy
        end

        for j = 1:N
            os += 8, "Nupdn", j
        end
        H = MPO(os, sites, splitblocks=true)

        # state = ["Up","Up","Up","Up","Up","Up","Up","Up"]
        state = [isodd(n) ? "Up" : "Dn" for n = 1:N]
        psi0 = MPS(sites, state)
        # psi0 = random_mps(sites, linkdims=1)

        nsweeps = 10
        maxdim = [800]
        cutoff = [1e-14]

        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

        return
    end
end
