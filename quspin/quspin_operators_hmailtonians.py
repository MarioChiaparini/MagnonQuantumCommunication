from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
import numpy as np

#define models paramters 
L = 4
J = 1.0
g = 0.809
h = 0.9045
Omega = 4.5

def drive(t, Omega):
  return np.cos(Omega*t)

drive_args = [Omega]

drive_args

basis = spin_basis_1d(L=L, a=1, kblock=0, pblock=1)

basis.blocks

x_field = [[g, i] for i in range(L)]
z_field = [[h, i] for i in range(L)]
J_nn = [[J, i, (i + 1) % L] for i in range(L)]  # PBC

J_nn

static = [["zz", J_nn], ["z", z_field]]
dynamic = [["x", x_field, drive, drive_args]]

H = hamiltonian(static, dynamic, static_fmt="dia", dtype=np.float64, basis=basis)

print(H)
print(H.toarray())