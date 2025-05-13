from qutip import *
import numpy as np
import matplotlib.pyplot as plt

# Parameters from your reference paper
ω_r = 2 * np.pi * 5.27e9       # Resonator frequency (5.27 GHz)
ω_m = 2 * np.pi * 5.405e9      # Magnon modes frequency (5.405 GHz)
g_mr = 2 * np.pi * 46e6        # Magnon-resonator coupling (46 MHz)
κ_m = 2 * np.pi * 0.98e6       # Magnon damping (~1 MHz)
g_mm = 2 * np.pi * 14.8e6  

# Drive pulse (Gaussian pulse on magnon 1)
def Ω(t, args):
    Ω0 = 2*np.pi*10e6
    t0, σ = 20e-9, 5e-9
    return Ω0 * np.exp(-(t - t0)**2 / (2 * σ**2))

N = 10 

a = tensor(destroy(N), qeye(N), qeye(N))   # Resonator mode
b1 = tensor(qeye(N), destroy(N), qeye(N))  # Magnon mode 1
b2 = tensor(qeye(N), qeye(N), destroy(N))

rho0 = ket2dm(tensor(basis(N,0), basis(N,0), basis(N,0)))

H = (ω_r * a.dag() * a +
     ω_m * (b1.dag() * b1 + b2.dag() * b2) +
     g_mr * (a.dag() * b1 + a * b1.dag() + a.dag() * b2 + a * b2.dag()))

H_drive = [b1 + b1.dag(), Ω]
H_total = [H, H_drive]
# Dissipation
c_ops = [np.sqrt(κ_m) * b1, np.sqrt(κ_m) * b2]
# Simulation times
tlist = np.linspace(0, 100e-9, 1000)
# Solve dynamics
result = mesolve(H_total, rho0, tlist, c_ops, [b1.dag()*b1, b2.dag()*b2, a.dag()*a])

plt.figure(figsize=(10, 6))
plt.plot(tlist*1e9, result.expect[0], label="Magnon mode 1")
plt.plot(tlist*1e9, result.expect[1], label="Magnon mode 2")
plt.plot(tlist*1e9, result.expect[2], label="Resonator")
plt.xlabel("Time (ns)")
plt.ylabel("Photon/Magnon number")
plt.legend()
plt.title("Coherent Magnon Dynamics")
plt.grid(True)
plt.tight_layout()
plt.show()