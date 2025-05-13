import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Define symbolic variables
a, b, w = sp.symbols('a b w')
a1, b1, pin, pout, wa, wb, v, k, d, s21, H = sp.symbols(
    'a_1 b_1 p_in p_out w_a w_b v k d s_21 H'
)
j = sp.I

# Define equations for the system
eq1 = sp.Eq(((w - wa + a1*j)*a)*j - (sp.sqrt(v)*pin)*j - v*a - sp.sqrt(k*v)*b - (d*b)*j, 0)
eq2 = sp.Eq(((w - wb + b1*j)*b)*j - (sp.sqrt(k)*pin)*j - k*b - sp.sqrt(k*v)*a - (d*a)*j, 0)
eq3 = sp.Eq(((w - wa + a1*j)*a)*j - (sp.sqrt(v)*pout)*j + v*a + sp.sqrt(k*v)*b - (d*b)*j, 0)
eq4 = sp.Eq(((w - wb + b1*j)*b)*j - (sp.sqrt(k)*pout)*j + k*b + sp.sqrt(k*v)*a - (d*a)*j, 0)
eq5 = sp.Eq(s21, (pout/pin) - 1)

# Solve for S21 parameter
sol1 = sp.solve(eq1, b)
solution1 = sp.solve(eq2.subs(b, sol1[0]), a)
sol2 = sp.solve(eq1, a)
solution2 = sp.solve(eq2.subs(a, sol2[0]), b)
new_eq = sp.Eq(eq1.lhs - eq3.lhs, 0)
solution3 = sp.solve(new_eq.subs({a: solution1[0], b: solution2[0]}), pout)
solution4 = sp.solve(eq5.subs(pout, solution3[0]), s21)
S = solution4[0]

# Parameters from the article's Table 1 (highest α case)
alpha_values = [
    (1.4e-5, 127.3e6),   # (α, g) in Hz
    (2.8e-2, 62.6e6)
]

# Convert GHz to Hz (1e9) and MHz to Hz (1e6)
wa_value = 5.33e9        # 5.33 GHz in Hz
beta = 4.69e-3           # Photon intrinsic damping (a1)
gamma_c = 24.997e6       # Photon extrinsic damping (v)

# Setup plots
plt.figure(figsize=(12, 8))

for idx, (alpha, g) in enumerate(alpha_values):
    # Convert parameters
    b1_val = alpha
    d_val = g / 1e9      # Convert g to GHz
    v_val = gamma_c / 1e9  # Convert to GHz

    # Substitute values into S
    s21_sub = S.subs({
        v: v_val,
        k: 0,            # Extrinsic magnon damping (γm=0)
        wa: wa_value / 1e9,  # Convert to GHz for calculation
        wb: wa_value / 1e9,  # At resonance ωm=ωc
        a1: beta,
        b1: b1_val,
        d: d_val
    })

    # Lambdify S21
    s21_func = sp.lambdify(w, sp.Abs(s21_sub), modules='numpy')

    # Frequency range around resonance
    freq = np.linspace(5.2, 5.4, 1000)  # GHz
    s21 = s21_func(freq)

    # Plot transmission spectra
    plt.subplot(2, 1, idx+1)
    plt.plot(freq, s21, label=f'α={alpha:.1e}')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('|S21|')
    plt.title(f'Transmission Spectrum at Resonance (α={alpha:.1e})')
    plt.grid(True)

plt.tight_layout()

# Contour plot for anti-crossing (example for α=2.8e-2)
alpha, g = alpha_values[1]
b1_val = alpha
d_val = g / 1e9

# Define H range to sweep through resonance
H_vals = np.linspace(1.1, 1.3, 100)  # Adjusted to cross ωc=5.33 GHz
wb_vals = 1.3828 + 3.2478 * H_vals   # Article's ωm vs H relation

# Frequency range
w_vals = np.linspace(5.2, 5.4, 500)  # GHz

# Compute S21 matrix
S_matrix = np.zeros((len(w_vals), len(H_vals)))
for i, H in enumerate(H_vals):
    wb = 1.3828 + 3.2478 * H
    s21_sub_h = S.subs({
        v: v_val,
        k: 0,
        wa: wa_value / 1e9,
        wb: wb,
        a1: beta,
        b1: b1_val,
        d: d_val
    })
    s21_func_h = sp.lambdify(w, sp.Abs(s21_sub_h), modules='numpy')
    S_matrix[:, i] = s21_func_h(w_vals)

# Plot contour
plt.figure(figsize=(10, 6))
plt.contourf(H_vals, w_vals, S_matrix, levels=50, cmap='jet')
plt.colorbar(label='|S21|')
plt.xlabel('H (arb. units)')
plt.ylabel('Frequency (GHz)')
plt.title('Anti-crossing Transition (α=2.8e-2)')
plt.show()