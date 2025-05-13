from mpl_toolkits.mplot3d import Axes3D
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
from matplotlib import pyplot as plt
import numpy as np
from qutip import Qobj
from qutip import coherent_dm, squeeze, fock_dm, wigner, displace

from mpl_toolkits.mplot3d import Axes3D

Jxy, Jzz = 1.0, 1.0
# define coupling strengths
hx, hy, hz = 1.0, 1.0, 1.0

# associate each coupling to the spin labeled 0
hx_list = [[hx,0],] # coupling: hx multiplies spin labeled 0
hy_list = [[hy,0],] # coupling: hy multiplies spin labeled 0
hz_list = [[hz,0],] # coupling: hz multiplies spin labeled 0

# associate a Pauli matrix to each coupling list
static_terms = [['x',hx_list], # assign coupling list hx_list to Pauli operator sigma^x
    		    ['y',hy_list], # assign coupling list hy_list to Pauli operator sigma^y
    		    ['z',hz_list], # assign coupling list hz_list to Pauli operator sigma^z
    		   ]
dynamic_terms = []

# system size / number of spins
L = 10
basis_ising = spin_basis_1d(L=L,)
Jzz_list = [[Jzz, j,j+1] for j in range(L-1)] # L-1 bonds
hz_list = [[hz,j] for j in range(L)] # L sites
hx_list = [[hx,j] for j in range(L)] # L sites

H_terms = [['zz',Jzz_list],
		   ['x',hx_list],
		   ['z',hz_list],
		  ]


H_Ising = hamiltonian(H_terms,[], basis=basis_ising)
E, V = H_Ising.eigh()
# ground state
psi_GS = V[:,0]
E_GS = E[0]
print('\nGS energy = {}'.format(E_GS) )
# first-excited state
psi_ex_1 = V[:,1]
E_exc = E[1]
print('\nfirst excited state energy = {}\n'.format(E_exc) )

# compute histogram
DOS, energy =  np.histogram(E, bins=50) # number of bins can be adjusted

#density of state plot qutip

print(DOS)






