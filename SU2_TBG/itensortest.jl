#use itensor code for benchmark
using ITensors, ITensorMPS
@time begin
    let
        # lx = 1
        # ly = 4
        # N = lx * ly
        include("readdata.jl")
        N = 2 * Nk
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
        # for iy = 1:ly-1
        #     os += "Cdagup", (lx - 1) * ly + iy, "Cup", (lx - 1) * ly + iy + 1
        #     os += "Cdagup", (lx - 1) * ly + iy + 1, "Cup", (lx - 1) * ly + iy
        #     os += "Cdagdn", (lx - 1) * ly + iy, "Cdn", (lx - 1) * ly + iy + 1
        #     os += "Cdagdn", (lx - 1) * ly + iy + 1, "Cdn", (lx - 1) * ly + iy
        # end

        # for j = 1:N
        #     os += 8, "Nupdn", j
        # end
        for inh=[1 6 7 8 31 36] #in eachindex(hlist)
            tsla1 = hlist[inh][1]#Lamda[1, 2]
            tsla2 = hlist[inh][1]'#Lamda[1, 33]
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
                                                os += tsla1[(ik1-1)*2+ia1, (ik2-1)*2+ia2] * tsla2[(iq1-1)*2+ib1, (iq2-1)*2+ib2], "Cdagup", (ik1 - 1) * 2 + ia1, "Cup", (ik2 - 1) * 2 + ia2, "Cdagdn", (iq1 - 1) * 2 + ib1, "Cdn", (iq2 - 1) * 2 + ib2
                                                os += tsla1[(ik1-1)*2+ia1, (ik2-1)*2+ia2] * tsla2[(iq1-1)*2+ib1, (iq2-1)*2+ib2], "Cdagdn", (ik1 - 1) * 2 + ia1, "Cdn", (ik2 - 1) * 2 + ia2, "Cdagup", (iq1 - 1) * 2 + ib1, "Cup", (iq2 - 1) * 2 + ib2
                                                os += tsla1[(ik1-1)*2+ia1, (ik2-1)*2+ia2] * tsla2[(iq1-1)*2+ib1, (iq2-1)*2+ib2], "Cdagdn", (ik1 - 1) * 2 + ia1, "Cdn", (ik2 - 1) * 2 + ia2, "Cdagdn", (iq1 - 1) * 2 + ib1, "Cdn", (iq2 - 1) * 2 + ib2
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        H = MPO(os, sites, splitblocks=true)

        # state = ["Up","Up","Up","Up","Up","Up","Up","Up"]
        Eki = zeros(1, Nk)
        for imp = 1:Nk
            state = [["Emp" for n = 1:((imp-1)*2)]; "Up"; ["Emp" for n = ((imp-1)*2+2):N]]
            # state = [isodd(n) ? "Up" : "Up" for n = 1:N]
            psi0 = MPS(sites, state)
            # psi0 = random_mps(sites, state, linkdims=10)

            nsweeps = 10
            maxdim = [800]
            cutoff = [1e-14]

            energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
            Eki[imp]=energy
        end
        return plot(Eki[1, 1:end], legend=false)
    end
end

imp=2
N=18
ts=[["Up" for n = 1:((imp-1)*2+1)]; "UpDn"; ["Up" for n = ((imp-1)*2+3):N]]
@show ts