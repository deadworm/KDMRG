#use itensor code for benchmark
using ITensors, ITensorMPS
@time begin
    let
        t = 1
        V = 1
        l = 8
        N = 4
        sites = siteinds("Fermion", l; conserve_qns=true)

        os = OpSum()
        for i = 1:l-1
            os += -t, "Cdag", i, "C", i+1
            os += -t, "Cdag", i+1, "C", i
            os += V, "N", i, "N", i+1
        end
        os += t, "Cdag", l, "C", 1
        os += t, "Cdag", 1, "C", l
        os += V, "N", l, "N", 1

        H = MPO(os, sites, splitblocks=true)

        states = [i<=N ? "1" : "0" for i in 1:l]
        psi0 = MPS(ComplexF64, sites, states)

        nsweeps = 10
        maxdim = [20]
        cutoff = [1e-14]

        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

        return
    end
end
