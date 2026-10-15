#use itensor code for benchmark
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using ITensors, ITensorMPS
@time begin
    let
        t = -1
        U = 1
        l = 2
        N = 3*l*l
        sites = siteinds("Electron", N; conserve_qns=false)

        os = OpSum()
        for ix = 1:l
            idx = mod(ix, l)+1
            for iy = 1:l
                idy = mod(iy, l)+1
                i = (ix-1)*l+iy
                iidx = (idx-1)*l+iy
                iidy = (ix-1)*l+idy
                iidxdy = (idx-1)*l+idy
                for inew in [iidx, iidy, iidxdy]
                    os += t, "Cdagup", (i-1)*3+1, "Cup", (inew-1)*3+1
                    os += t, "Cdagup", (i-1)*3+2, "Cup", (inew-1)*3+2
                    os += t, "Cdagup", (i-1)*3+3, "Cup", (inew-1)*3+3
                    os += t, "Cdagdn", (i-1)*3+1, "Cdn", (inew-1)*3+1
                    os += t, "Cdagdn", (i-1)*3+2, "Cdn", (inew-1)*3+2
                    os += t, "Cdagdn", (i-1)*3+3, "Cdn", (inew-1)*3+3

                    os += t, "Cdagup", (inew-1)*3+1, "Cup", (i-1)*3+1
                    os += t, "Cdagup", (inew-1)*3+2, "Cup", (i-1)*3+2
                    os += t, "Cdagup", (inew-1)*3+3, "Cup", (i-1)*3+3
                    os += t, "Cdagdn", (inew-1)*3+1, "Cdn", (i-1)*3+1
                    os += t, "Cdagdn", (inew-1)*3+2, "Cdn", (i-1)*3+2
                    os += t, "Cdagdn", (inew-1)*3+3, "Cdn", (i-1)*3+3
                end
            end
        end

        Mat=[0 1 -1; -1 0 1; 1 -1 0]
        for il=1:(l*l)
            for im1=1:3
                for in1=1:3
                    for im2=1:3
                        for in2=1:3
                            os += U*Mat[im1, in1]*Mat[im2, in2], "Cdagup", (il-1)*3+im1, "Cup", (il-1)*3+in1, "Cdagup", (il-1)*3+im2, "Cup", (il-1)*3+in2
                            os += U*Mat[im1, in1]*Mat[im2, in2], "Cdagup", (il-1)*3+im1, "Cup", (il-1)*3+in1, "Cdagdn", (il-1)*3+im2, "Cdn", (il-1)*3+in2
                            os += U*Mat[im1, in1]*Mat[im2, in2], "Cdagdn", (il-1)*3+im1, "Cdn", (il-1)*3+in1, "Cdagup", (il-1)*3+im2, "Cup", (il-1)*3+in2
                            os += U*Mat[im1, in1]*Mat[im2, in2], "Cdagdn", (il-1)*3+im1, "Cdn", (il-1)*3+in1, "Cdagdn", (il-1)*3+im2, "Cdn", (il-1)*3+in2
                        end
                    end
                end
            end
        end

        H = MPO(os, sites, splitblocks=true)
        states = [isodd(n) ? "Up" : "Dn" for n in 1:N]
        psi0 = random_mps(sites, linkdims=1)

        nsweeps = 20
        maxdim = [1000]
        cutoff = [1e-14]

        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

        return
    end
end
