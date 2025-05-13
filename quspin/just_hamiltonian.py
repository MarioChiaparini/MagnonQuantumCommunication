from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d


basis_1 = spin_basis_1d(L=1) # single spin-1/2 particle / qubit / two-level system
print(basis_1)
print(basis_1.blocks) # blocks of the basis
print(basis_1.Ns) # number of states in the basis
dynamic_terms = []

hx, hy, hz = 1.0, 1.0, 1.0

hx_list = [[hx, 0]]
hy_list = [[hy, 0]]  # use site 0 instead of 1
hz_list = [[hz, 0]]# coupling: hz multiplies spin labeled 2

# associate a Pauli matrix to each coupling list
static_terms = [['x',hx_list],
                ['y',hy_list],
                ['z',hz_list]]

H_1 = hamiltonian(static_terms, dynamic_terms, basis=basis_1)
print(H_1)
print(H_1.toarray()) # print the matrix representation of the Hamiltonian