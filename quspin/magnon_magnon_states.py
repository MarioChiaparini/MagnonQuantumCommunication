from quspin.basis import spin_basis_1d
import numpy as np

# 1. Define the 4-spin Hilbert space restricted to states with exactly 2 up-spins (magnons)
L = 4
basis_2mag = spin_basis_1d(L=L, Nup=2)

# 2. Initialize a state vector of zeros having dimension equal to the number of allowed basis states.
psi_W = np.zeros(basis_2mag.Ns, dtype=np.complex128)

# 3. Define the state coefficients for each basis state.
# Here we follow the form given in the paper, e.g.,
#    |W> = 1/(2√2) * [ √2|1100> + √2|1010> - |1001> + |0011> - √2|0101> + |0110> ]
# You can adjust the prefactors if your definition of spin "up" and "down" differs.
norm_factor = 1.0 / (2 * np.sqrt(2))
coefficients = {
    "1100": norm_factor * np.sqrt(2),
    "1010": norm_factor * np.sqrt(2),
    "1001": -norm_factor,
    "0011": norm_factor,
    "0101": -norm_factor * np.sqrt(2),
    "0110": norm_factor,
}

# 4. Fill the state vector using the `index` method from the basis.
for state, coeff in coefficients.items():
    # The `index` method maps the bitstring to the array index.
    idx = basis_2mag.index(state)
    psi_W[idx] = coeff

# 5. Normalize the state vector (if not already normalized).
psi_W /= np.linalg.norm(psi_W)

print("Representation of the two-magnon state |W> in the restricted basis:")
print(psi_W)
entropy_result = basis_2mag.ent_entropy(psi_W, sub_sys_A=[0,1])
print("Von Neumann entropy for subsystem A:", entropy_result['Sent_A'])