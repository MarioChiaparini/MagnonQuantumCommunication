using Sunny, LinearAlgebra, GLMakie, StaticArrays

latvecs = lattice_vectors(1, 1, 2, 90, 90, 90)
positions = [[0, 0, 0], [0, 1/2 - 1e-4, 0]]

positions = [
    [0.0, 0.0, 0.0],
    [0.1, 0.0, 0.0],
    [0.2, 0.0, 0.0],
    [0.3, 0.0, 0.0],
    [0.4, 0.0, 0.0],
    [0.5, 0.0, 0.0],
    [0.6, 0.0, 0.0],
    [0.7, 0.0, 0.0],
    [0.8, 0.0, 0.0],
    [0.9, 0.0, 0.0],
]

crystal = Crystal(latvecs, positions)

view_crystal(crystal)