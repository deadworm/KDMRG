#dmrg test for Heisenberg model
using TensorKit, MPSKit
using KrylovKit
using Random
using Combinatorics

include("mpstools.jl")

l = 4
n = 2
bond_dim = 4
phySpace = Vect[FermionParity⊠U1Irrep]((0, 0) => 1, (1, 1) => 1)

m = initial_mps(l, n, phySpace, bond_dim);

configtest=[0 0 1 1;
0 1 0 1;
0 1 1 0;
1 0 0 1;
1 0 1 0;
1 1 0 0];
normtest=0;
cm = zeros(1, 6);
for id = 1:6
    config = configtest[id, :]
    cm[id] = get_config(m, config)
    println("config: ", config, " coefficient: ", cm[id])
    normtest += abs(cm[id])^2
end
println("normtest: ", normtest)

ostringtest=[0.5 3 2 -3 -2;
0.5 2 3 -3 -2;
0.5 0 -3 1 0;
0.5 3 0 -3 0;]

olocal=observable_config(ostringtest, m, cm[3], configtest[3, :])