#use itensor code for benchmark
using ITensors, ITensorMPS
@time begin
    let
        # lx = 4
        # ly = 4
        # N = lx * ly
        N = 2 * Nk
        sites = siteinds("Electron", N; conserve_qns=true)

        os = OpSum()
        tsla1 = Lamda[1, 2]
        tsla2 = Lamda[1, 33]
        for ik1 = 1:Nk
            for ik2 = 1:Nk
                for ia1 = 1:2
                    for ia2 = 1:2
                        for iq1 = 1:Nk
                            for iq2 = 1:Nk
                                for ib1 = 1:2
                                    for ib2 = 1:2
                                        if tsla1[(ik1-1)*2+ia1, (ik2-1)*2+ia2] != 0 && tsla2[(iq1-1)*2+ib1, (iq2-1)*2+ib2] != 0
                                            os += tsla1[(ik1-1)*2+ia1, (ik2-1)*2+ia2] * tsla2[(iq1-1)*2+ib1, (iq2-1)*2+ib2], "Cdagup", (ik1 - 1) * 2 + ia1, "Cup", (ik2 - 1) * 2 + ia2, "Cdagup", (iq1 - 1) * 2 + ib1, "Cup", (iq2 - 1) * 2 + ib2
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        # os = OpSum()
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
        # for iy = 1:ly-1
        #     os += "Cdagup", (lx - 1) * ly + iy, "Cup", (lx - 1) * ly + iy + 1
        #     os += "Cdagup", (lx - 1) * ly + iy + 1, "Cup", (lx - 1) * ly + iy
        #     os += "Cdagdn", (lx - 1) * ly + iy, "Cdn", (lx - 1) * ly + iy + 1
        #     os += "Cdagdn", (lx - 1) * ly + iy + 1, "Cdn", (lx - 1) * ly + iy
        # end

        # for j = 1:N
        #     os += 8, "Nupdn", j
        # end
        H = MPO(os, sites, splitblocks=true)

        # state = ["Up","Up","Up","Up","Up","Up","Up","Up"]
        state = ["Up"; ["Emp" for n = 1:N-1]]
        # state = [isodd(n) ? "Up" : "Dn" for n = 1:N]
        psi0 = MPS(sites, state)
        # psi0 = random_mps(sites, linkdims=1)

        nsweeps = 10
        maxdim = [800]
        cutoff = [1e-14]

        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

        return
    end
end
