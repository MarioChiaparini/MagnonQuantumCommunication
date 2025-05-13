using SymPy, Plots

# Define symbolic variables
a, b, w = symbols("a b w")
a1, b1, pin, pout, wa, wb, v, k, d, s21 = symbols("a₁ b₁ p_{in} p_{out} w_a w_b v k d s_{21}")
j = im

# Define equations
eq1 = Eq(((w - wa + a1*j)*a)*j - sqrt(v)*pin*j - v*a - sqrt(k*v)*b - d*b*j, 0)
eq2 = Eq(((w - wb + b1*j)*b)*j - sqrt(k)*pin*j - k*b - sqrt(k*v)*a - d*a*j, 0)
eq3 = Eq(((w - wa + a1*j)*a)*j - sqrt(v)*pout*j + v*a + sqrt(k*v)*b - d*b*j, 0)
eq5 = Eq(s21, pout / pin - 1)

# Solve the coupled equations
sol1 = solve(eq1, b)[1]
a_sol = solve(subs(eq2, b => sol1), a)[1]
sol2 = solve(eq1, a)[1]
b_sol = solve(subs(eq2, a => sol2), b)[1]

# Substitute into difference equation
diff = eq1.lhs - eq3.lhs
new_eq = Eq(diff, 0)
new_eq_subs = subs(new_eq, Dict(a => a_sol, b => b_sol))
pout_sol = solve(new_eq_subs, pout)[1]
s21_sol = solve(subs(eq5, pout => pout_sol), s21)[1]

# Numeric parameter values
d_val = 0.0541
k_val = 0.0001
v_val = 0.04
wa_val = 5.275
b1_val = 0.0714
a1_val = 0.028

# Substitute numeric values into expression
s21_num_expr = subs(s21_sol, Dict(
    v => v_val,
    k => k_val,
    wa => wa_val,
    b1 => b1_val,
    a1 => a1_val,
    d => d_val
))

# Evaluation grid
H_vals = range(1.050, 1.350, length=101)
w_vals = range(4.5, 6.0, length=1001)

# Compute |S21| over the grid
Z = [abs(subs(s21_num_expr, Dict(w => wv, wb => 1.38284384 + 3.24789744 * H))) for wv in w_vals, H in H_vals]

# Plot
plot_title = "Contour Plot of s21"
contour = contourf(H_vals, w_vals, Z, cmap = :jet, levels = range(0, 1.5, length=501),
    xlabel = "H", ylabel = "w", colorbar_title = "|S₂₁|", title = plot_title)

# Save output to file
savefig("s21_contour_plot.png")