#This is a test for multi-virtual legs mpo

Rm=Rep[U₁ × SU₂]
phySpace = Rm((0, 0) => 1, (1, 1/2) => 1, (2, 0) => 1)
dicop,dicqn,dicbf=op2ten(phySpace)

using MPSKitModels
using MPSKit
mpo=hubbard_model(ComplexF64, Irrep[U₁], Irrep[SU₂], FiniteChain(2);t=1, U=0, mu=0, n=1)
dmpo1=convert(MPSKitModels.SparseMPO, mpo)
ts=dmpo1.O

dmpo2=convert(MPSKitModels.DenseMPO, mpo)
a1=convert(Array,dmpo2.opp[1])
a2=convert(Array,dmpo2.opp[2])

blspace = Rep[U₁ × SU₂]((0, 0)=>1)
brspace = Rep[U₁ × SU₂]((0, 0)=>2, (1, 1/2)=>1, (-1, 1/2)=>1)
t1=TensorMap(convert(Array,dmpo2.opp[1])[1,:,:,:], blspace ⊗ phySpace, phySpace ⊗ brspace)
tt1=convert(Array,t1)

t2=TensorMap(convert(Array,dmpo2.opp[2])[:,:,:,2], brspace ⊗ phySpace, phySpace ⊗ blspace)
tt2=convert(Array,t2)

space(dmpo2.opp[1])
space(dmpo2.opp[2])

@tensor mpot[p1 p2; p1p p2p]:= dicop["cdf"][p1; p1p b] * dicop["c"][b p2; p2p]
# @tensor mpotp[p1 p2; p1p p2p]:= -dicop["cf"][b p1; p1p] * dicop["cd"][p2; p2p b]
# mpot'≈mpotp
tsh=mpot'
@tensor tsh=isomorphism(blspace * codomain(tsh), codomain(tsh))*tsh*isomorphism(domain(tsh), domain(tsh) * blspace)
U, S, V = tsvd(tsh, (1, 2, 4), (3, 5, 6), trunc=truncbelow(1e-10))
ts1=permute(U*S,(1,2),(3,4))
ts2=permute(V,(1,2),(3,4))
tts1=convert(Array,ts1)
tts2=convert(Array,ts2)

tt1=tt1[:,:,:,3:6]
tt2=tt2[3:6,:,:,:]
@tensor T1[1,2,4,3,5,6]:=tt1[1,2,3,-1]*tt2[-1,4,5,6]
@tensor T2[1,2,4,3,5,6]:=tts1[1,2,3,-1]*tts2[-1,4,5,6]
T1≈-T2/sqrt(2)

convert(Array,hopping)≈T2[1,:,:,:,:,1]

space(U)
space(S)
space(V)
using LinearAlgebra
sort(diag(convert(Array,S)))

using MPSKitModels

T=ComplexF64
particle_symmetry = U1Irrep
spin_symmetry = SU2Irrep
hopping = e_plusmin(T, particle_symmetry, spin_symmetry) + e_minplus(T, particle_symmetry, spin_symmetry)
interaction_term = e_number_updown(T, particle_symmetry, spin_symmetry)
N = e_number(T, particle_symmetry, spin_symmetry)

lattice=FiniteChain(3)
t=1;mu=1;U=1;

E=scalartype(hopping)
T=TensorMap
L=3

function _find_free_channel(data::Array{Union{E,T},3},loc)::Tuple{Int,Array{Union{E,T},3}} where {E<:Number,T<:AbstractTensorMap}
hit = findfirst(map(x -> _is_free_channel(data, loc, x), 2:(size(data, 2) - 1)))
#hit = findfirst(ismissing.(data[loc,1,2:end-1]));
if isnothing(hit)
ndata = fill!(Array{Union{E,T},3}(undef, size(data, 1), size(data, 2) + 1,
                  size(data, 2) + 1),
zero(E))
ndata[:, 1:(end - 1), 1:(end - 2)] .= data[:, :, 1:(end - 1)]
ndata[:, 1:(end - 2), end] .= data[:, 1:(end - 1), end]
ndata[:, end, end] .= data[:, end, end]
return size(data, 2), ndata
else
return hit + 1, data
end
end
_iszeronumber(x::Number) = iszero(x)
_iszeronumber(x::AbstractTensorMap) = false
function _is_free_channel(data, loc, channel)
    return all(_iszeronumber, data[mod1(loc, end), :, channel])
end

data = fill!(Array{Union{E,T},3}(undef, L, 2, 2), zero(E))
data[:, 1, 1] .= one(E)
data[:, end, end] .= one(E)
data::Array{Union{E,T},3}
opps=h31
opp1=opps.opps[1]
opp2=opps.opps[2]
op11=opp1.opp[1]
op12=opp1.opp[2]
op21=opp2.opp[1]
op22=opp2.opp[2]

op11==op21==t_mpo[1]
op12==op22==t_mpo[2]
opp1.opp==opp2.opp
opp1.inds
opp2.inds

using TupleTools
using MPSKit

lop=LocalOperator(hopping,FiniteChain(3)[1], FiniteChain(3)[2])









op1*op2
space(op1)
space(op2)
aop=convert(Array,op)

for opp in opps.opps
    linds = linearize_index.(opp.inds)
    mpo = opp.opp

    if length(mpo) == 1
        if data[linds[1], 1, end] == zero(E)
            data[mod1(linds[1], L), 1, end] = mpo[1]
        else
            data[mod1(linds[1], L), 1, end] += mpo[1]
        end
        continue
    end

    start, stop = first(linds), last(linds)
    hit, data = _find_free_channel(data, start)

    data[mod1(start, L), 1, hit] = mpo[1]
    for site in (start + 1):(stop - 1)
        mpo_ind = findfirst(linds .== site)
        o = isnothing(mpo_ind) ? one(E) : mpo[mpo_ind]

        if length(lattice(opps)) > 1 && _is_free_channel(data, site, hit)
            data[mod1(site, L), hit, hit] = o
        else
            nhit, data = _find_free_channel(data, site)
            data[mod1(site, L), hit, nhit] = o
            hit = nhit
        end
    end

    data[mod1(stop, L), hit, end] = mpo[end]
end

@mpoham h31=hopping{FiniteChain(3)[1], FiniteChain(3)[2]}

ops=h31.opp
ops1p=convert(Array,ops[1])
ops2p=convert(Array,ops[2])

ops1==ops1p
ops2==ops2p

+hopping{FiniteChain(3)[2], FiniteChain(3)[3]}
h32=sum(nearest_neighbours(lattice)) do (i, j)
 return @mpoham hopping{i,j}
end

(Irrep[SU₂ × U₁](1/2, 1)⊗Irrep[SU₂ × U₁](1/2, -1)...,)
u1=convert(Array,mpo[1])
bospace = Rep[SU₂ × U₁]((0, 0)=>1)
blspace = Rep[SU₂ × U₁]((1/2, 1)=>1)
brspace = Rep[SU₂ × U₁]((1, 2)=>1)
Rm=Rep[SU₂×U₁]
phySpace = Rm((0, 0) => 1, (1 / 2, 1) => 1, (0, 2) => 1)
MPO1 = TensorMap(ones, blspace ⊗ phySpace, phySpace ⊗ brspace)
# blocks(MPO1)[Irrep[SU₂ × U₁](1,2)].=-1
ampo=convert(Array,MPO1)
at=zeros(2,4,4,3)
at[1,:,:,1]=[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
at[2,:,:,1]=[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
at[1,:,:,2]=[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
at[2,:,:,2]=[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
at[1,:,:,3]=[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
at[2,:,:,3]=[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
# at[4,:,:,4]=[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
MPO2 = TensorMap(at, blspace ⊗ phySpace, phySpace ⊗ brspace)
MPO3 = TensorMap(at[1,:,:,1], bospace ⊗ phySpace, phySpace ⊗ bospace)