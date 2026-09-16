#For initial random mps
function rnd_mps(l, bond_dim)
    B = Vector{typeof(zb)}(undef, l + 1)
    B[1] = zb
    # 构造中间虚键：B[i + 1] 对应第 i 个站点之后的键
    for i in 1:(l-1)
        qmin = max(0, n - 2 * (l - i))
        qmax = min(2 * i, n)
        charges = collect(qmin:qmax)
        nsectors = length(charges)
        # 将最大维数尽量平均分配给各个允许的电荷扇区
        base, remainder = divrem(bond_dim, nsectors)
        sector_dims = [
            (q % 2, q) => base + (j <= remainder ? 1 : 0)
            for (j, q) in enumerate(charges)
        ]
        B[i+1] = Vect[FermionParity⊠U1Irrep](sector_dims...)
    end
    # 固定右边界，从而保证整个 MPS 的总粒子数为 n
    B[end] = nb
    tensors = Vector{TensorMap}(undef, l)
    for i in 1:l
        tensors[i] = randn(ComplexF64, B[i] ⊗ phySpace ← B[i+1])
    end
    m = FiniteMPS([tensors[i] for i in 1:l])
    m = changebonds(m, RandExpand(trscheme=truncdim(bond_dim)))
    m = changebonds(m, SvdCut(trscheme=truncdim(bond_dim)))
    @show [dim(space(m[i], 3)) for i in 1:l]
    return m
end

#For initial direct-product mps
function prod_mps(config)
    B = Vector{typeof(zb)}(undef, l + 1)
    B[1] = zb
    tensors = Vector{TensorMap}(undef, l)
    for i = 1:l
        pindex=s2p[(config[2*i-1],config[2*i])]
        B[i+1] = fuse(B[i], phylist[pindex])
        tenow = zeros(Float64, B[i] ⊗ phySpace ← B[i+1])
        if pindex == 3
            tenow.data[2] = 1.0
        else
            tenow.data[1] = 1.0
        end
        tensors[i]=tenow
    end
    mf = FiniteMPS([tensors[i] for i in 1:l])
    return mf
end

#For computing observables [coefficient, site1, site2, ...] for each term
function observable_config(ostring, m, cm, config, Dall)
    oterms=size(ostring, 1)
    olength=size(ostring, 2)
    oc=ostring[:, 1]
    Dc = Dict{Vector{Int},ComplexF64}()
    for it=1:oterms
        nconfig=copy(config)
        for il=0:(olength-2)
            oind=Int(ostring[it, end-il])
            if oind==0
                continue
            elseif oind>0&&nconfig[oind]==0
                oc[it]=oc[it]*(-1)^sum(@view nconfig[1:(oind-1)])
                nconfig[oind]=1
            elseif oind<0&&nconfig[-oind]==1
                oc[it]=oc[it]*(-1)^sum(@view nconfig[1:(-oind-1)])
                nconfig[-oind]=0
            else
                oc[it]=0
                break
            end
        end
        if oc[it]!=0
            if !haskey(Dc, nconfig)
                Dc[nconfig]=oc[it]
            else
                Dc[nconfig]+=oc[it]
            end
        end
    end
    osum=0
    for cf in keys(Dc)
        if haskey(Dall, cf)
            osum+=Dc[cf] * Dall[cf]
        else
            mf=prod_mps(cf)
            Dall[cf] = MPSKit.dot(mf, m)
            osum+=Dc[cf] * Dall[cf]
        end
    end
    olocal=osum / cm
    return olocal
end

#For mps gradient-normalize update
function gn_update(m, wm, em, gm, Dm, dm)
    nm=Vector{TensorMap}(undef, l)
    configs=collect(keys(Dm))
    for il=1:l
        nm[il] = copy(m[il])
    end
    for cf in configs
        (cn, _, enm, gradm) = Dm[cf]
        for il=1:l
            δm = cn * (gradm[il]-gm[il]) * (enm-em) / wm
            nm[il] = nm[il] - dm * δm
        end
    end
    δm_norm=0
    for il=1:l
        δm_norm+=norm(nm[il]-m[il])
    end
    @show δm_norm
    return FiniteMPS([nm[il] for il in 1:l])
end

# A Hermitian matrix plus a scalar diagonal shift. Keeping the shift as an
# operator avoids allocating a second (potentially GPU-resident) Gram matrix.
struct SRShiftedHermitian{T,M<:AbstractMatrix{T},R<:Real} <: AbstractMatrix{T}
    gram::M
    shift::R
end

function SRShiftedHermitian(gram::M, shift::R) where {T,M<:AbstractMatrix{T},R<:Real}
    return SRShiftedHermitian{T,M,R}(gram, shift)
end

Base.size(A::SRShiftedHermitian) = size(A.gram)
Base.axes(A::SRShiftedHermitian) = axes(A.gram)

function LinearAlgebra.mul!(y::AbstractVector, A::SRShiftedHermitian, x::AbstractVector)
    mul!(y, A.gram, x)
    y .+= A.shift .* x
    return y
end

# Flatten the symmetry blocks of all log-derivative tensors into columns. The
# resulting dense matrix lets BLAS/cuBLAS replace O(nsamples^2 * nsites) many
# small TensorKit contractions with one large matrix multiplication.
function sr_pack_gradients(gm, Dm, wm, em, ::Type{RT}) where {RT<:AbstractFloat}
    isempty(Dm) && throw(ArgumentError("sr_update requires at least one sample"))
    isempty(gm) && throw(ArgumentError("sr_update requires at least one MPS tensor"))
    isfinite(wm) && wm > 0 || throw(ArgumentError("wm must be finite and positive"))

    configs = collect(keys(Dm))
    CT = Complex{RT}
    ST = sectortype(gm[1])
    layout = Vector{Vector{Tuple{ST,UnitRange{Int}}}}(undef, length(gm))
    offset = 0
    for il in eachindex(gm)
        site_layout = Tuple{ST,UnitRange{Int}}[]
        for sector in blocksectors(gm[il])
            block_length = length(block(gm[il], sector))
            range = (offset + 1):(offset + block_length)
            push!(site_layout, (sector, range))
            offset += block_length
        end
        layout[il] = site_layout
    end

    gradients = Matrix{CT}(undef, offset, length(configs))
    energies = Vector{CT}(undef, length(configs))
    for (ik, cf) in enumerate(configs)
        sample_weight = Dm[cf][1]
        isfinite(sample_weight) && sample_weight >= 0 ||
            throw(ArgumentError("sample weights must be finite and nonnegative"))
        scale = sqrt(RT(sample_weight / wm))
        energies[ik] = scale * (Dm[cf][3] - em)
        sample_gradient = Dm[cf][4]
        length(sample_gradient) == length(gm) ||
            throw(DimensionMismatch("sample and mean gradients have different lengths"))

        for il in eachindex(gm)
            for (sector, range) in layout[il]
                source = vec(block(sample_gradient[il], sector))
                mean_source = vec(block(gm[il], sector))
                length(source) == length(range) ||
                    throw(DimensionMismatch("incompatible TensorMap blocks in sample gradient"))
                destination = @view gradients[range, ik]
                @. destination = scale * (source - mean_source)
            end
        end
    end
    return gradients, energies, layout
end

function sr_solve(gradients, energies, regularization, reltol, maxiter, use_gpu)
    lk = length(energies)
    RT = real(eltype(gradients))
    mean_diagonal = real(sum(abs2, gradients)) / lk
    shift = RT(regularization) * max(mean_diagonal, eps(RT))

    if use_gpu
        device_gradients = CUDA.CuArray(gradients)
        device_energies = CUDA.CuArray(energies)
        gram = Hermitian(adjoint(device_gradients) * device_gradients)
        coefficients, history = cg(
            SRShiftedHermitian(gram, shift), device_energies;
            reltol=RT(reltol), maxiter=maxiter, log=true,
        )
        update = Array(device_gradients * coefficients)
    else
        gram = Hermitian(adjoint(gradients) * gradients)
        coefficients, history = cg(
            SRShiftedHermitian(gram, shift), energies;
            reltol=RT(reltol), maxiter=maxiter, log=true,
        )
        update = gradients * coefficients
    end
    return update, history, shift
end

function sr_cuda_enabled(backend)
    backend in (:auto, :cpu, :gpu) ||
        throw(ArgumentError("backend must be :auto, :cpu, or :gpu"))
    backend === :cpu && return false

    functional = CUDA.functional()
    backend === :gpu && !functional &&
        throw(ArgumentError("backend=:gpu requested, but CUDA is not functional"))
    return functional
end

"""
    sr_update(m, wm, em, gm, Dm, dt; backend=:auto, precision=Float64,
              regularization=1e-4, reltol=1e-5, maxiter=5length(Dm))

Apply a stochastic-reconfiguration update. With `backend=:auto`, the dense Gram
matrix construction, conjugate-gradient solve, and parameter update run on an
NVIDIA GPU when CUDA is functional, and otherwise fall back to CPU BLAS.

Set `backend=:gpu` to require CUDA or `backend=:cpu` to disable it. `Float64`
preserves the previous numerical precision; `precision=Float32` is usually
faster on consumer GPUs and uses half as much device memory.
"""
function sr_update(
    m, wm, em, gm, Dm, dt;
    backend=:auto,
    precision::Type{<:AbstractFloat}=Float64,
    regularization::Real=1e-4,
    reltol::Real=1e-5,
    maxiter::Integer=5length(Dm),
)
    regularization >= 0 || throw(ArgumentError("regularization must be nonnegative"))
    reltol > 0 || throw(ArgumentError("reltol must be positive"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    length(m) == length(gm) ||
        throw(DimensionMismatch("MPS and mean gradient have different lengths"))

    use_gpu = sr_cuda_enabled(backend)
    gradients, energies, layout = sr_pack_gradients(gm, Dm, wm, em, precision)
    update, cglog, shift = sr_solve(
        gradients, energies, regularization, reltol, Int(maxiter), use_gpu,
    )
    println(
        "SR backend: ", use_gpu ? "CUDA" : "CPU",
        ", samples: ", length(Dm),
        ", parameters: ", size(gradients, 1),
        ", diagonal shift: ", shift,
    )
    @show cglog

    nm = [copy(m[il]) for il in eachindex(m)]
    for il in eachindex(nm)
        for (sector, range) in layout[il]
            destination = vec(block(nm[il], sector))
            destination .-= dt .* @view(update[range])
        end
    end
    # ∂m_norm = sum(norm(nm[il] - m[il]) for il in eachindex(m))
    @show [norm(nm[il] - m[il]) for il in eachindex(m)]
    return FiniteMPS(nm)
end

#For computing the entanglement entropy
function mps_entropy(m, l)
    entropy_list = Vector{Float64}(undef, l-1)
    for b in 1:(l-1)
        λ = entanglement_spectrum(m, b)
        # collect singular values from all symmetry sectors
        λ_all = vcat(values(λ)...)
        # probabilities
        p = λ_all .^ 2
        p ./= sum(p)   # safe normalization
        # von Neumann entropy
        S = -sum(p .* log.(p .+ 1e-15))
        entropy_list[b] = S
    end
    return entropy_list
end
