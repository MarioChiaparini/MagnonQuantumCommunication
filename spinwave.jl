using Sunny, GLMakie, StaticArrays

crystal_full = Crystal("YIG.cif"; symprec=0.01)
crystal = subcrystal(crystal_full, "Fe")  # focus only on magnetic Fe atoms
view_crystal(crystal)

sys = System(crystal, [1 => Moment(s=5/2, g=2)], :SUN)
J1 = 1.0  # nearest-neighbor exchange
set_exchange!(sys, J1, Bond(1, 2, [0, 0, 0]))

D = -0.001  # very small easy-axis anisotropy
set_onsite_coupling!(sys, S -> D * S[3]^2, 1)

dimers = [(1, 2), (3, 4)]
esys = Sunny.EntangledSystem(sys, dimers)

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys)

formfactors = [1 => FormFactor("Fe")]
measure = ssf_perp(sys; formfactors)
swt = SpinWaveTheory(sys; measure)

fwhm = 0.1  # smaller broadening for YIG magnon peaks
qpts = q_space_path(crystal, [[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0]], 400)
energies = range(0.0, 100.0, 400)  # YIG magnon range can go up to GHz → meV scale
res = intensities(swt, qpts; energies, kernel=gaussian(; fwhm))
plot_intensities(res)