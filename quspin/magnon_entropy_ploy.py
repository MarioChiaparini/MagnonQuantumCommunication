from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
from quspin.tools.measurements import ent_entropy
from quspin.operators import exp_op
import numpy as np
from scipy.sparse.linalg import expm
import matplotlib.pyplot as plt
# DOI 10.1007/s11128-011-0252-z
dynamic_terms = []
L = 4
Jzz = 1.0   # Ising interaction strength
hx = 1.0    # transverse field strength
hz = 1.0


Jzz_list = [[Jzz, j, j+1] for j in range(L-1)]  # nearest-neighbor ZZ interactions
hx_list = [[hx, j] for j in range(L)]             # single-spin X terms
hz_list = [[hz, j] for j in range(L)]  


static_terms = [["zz", Jzz_list],
                ["x", hx_list],
                ["z", hz_list]]

basis_ising = spin_basis_1d(L=L)

N_time_steps = 101
times, dt = np.linspace(0.0,5.0,N_time_steps, retstep=True)

H_Ising = hamiltonian(static_terms, dynamic_terms, basis=basis_ising)
psi_cat = np.zeros(basis_ising.Ns, dtype=np.complex128)

psi_cat[basis_ising.index('0'*L)] = 1.0
psi_cat[basis_ising.index('1'*L)] = 1.0

psi_cat /= np.linalg.norm(psi_cat)

N_timesteps = 101
times, dt = np.linspace(0.0, 5.0, N_timesteps, retstep=True)
U_dt = expm(-1j * dt * H_Ising.tocsc())


psi_cat_t = np.zeros((basis_ising.Ns, N_timesteps), dtype=np.complex128)
E_t = np.zeros(N_timesteps, dtype=np.float64)
Sent_t = np.zeros(N_timesteps, dtype=np.float64)

psi_cat_t[:, 0] = psi_cat.copy()
E_t[0] = H_Ising.expt_value(psi_cat).real
Sent_t[0] = basis_ising.ent_entropy(psi_cat)['Sent_A']

for j in np.arange(N_timesteps - 1):
    # Evolve the state one time step forward.
    psi_cat_t[:, j+1] = U_dt @ psi_cat_t[:, j]
    # Measure energy expectation value.
    E_t[j+1] = H_Ising.expt_value(psi_cat_t[:, j+1]).real
    # Compute entanglement entropy for, say, a half chain.
    Sent_t[j+1] = basis_ising.ent_entropy(psi_cat_t[:, j+1])['Sent_A']



plt.plot(times, E_t/basis_ising.L, label='$E(t)/L$')
#plt.plot(times, Sent_t, label='$S_\mathrm{ent}^\mathrm{vN}(t)/L_A$')
plt.xlabel('time $Jt$')
plt.legend()
plt.show()