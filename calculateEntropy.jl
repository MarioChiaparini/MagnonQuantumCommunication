using LinearAlgebra
using Arpack
using Plots


function compute_rdmA(Psiab)
    @assert ndims(Psiab) == 2
    @assert size(Psiab, 1) < 2^13

end

function compute_entropy(dm)
    p = eigvals(dm) 
    entropy = 0.0
    #entropy = - p'*log2.(p) 
    for n in 1:size(p,1)
        if abs(p[n]) > 1e-12
            entropy = entropy - p[n]*log2(p[n])
        end
    end
    entropy
end


function normalizeVec(Vec)
    Vec = Vec/sqrt(abs(Vec'*Vec)) 
end

# here we implement the outer product
function separable(N)
    PsiInitial= reshape( (randn(2) + im*randn(2))*(randn(2) + im*randn(2))', 2^2 )
    for i in 3:N
        PsiInitial = reshape( PsiInitial*(randn(2) + im*randn(2))', 2^i )
    end
    return normalizeVec(PsiInitial)
end

function buildIsing(theta=pi/4) # Ising model with transverse magnetic field h (critical theta=pi/4 by default)
    X = [0. 1; 1 0]
    Z = [1. 0; 0 -1]
    E = diagm(0=>ones(2))
    XX = kron(X,X)
    HXX = XX
    HZ = kron(Z,E) + kron(E,Z)
    H2 = -(cos(theta)*XX + sin(theta)/2*HZ)
    return H2
end

function multiplyHPsi(Psi, H2) 
    L = length(Psi)  # Dimension of the vector space
    N = convert(Int64,log2(L))  # Number of spins
    HPsi = zeros(L)
    for n=1:N                              # This multiplies by the Hamiltonian H
        Psi = reshape(Psi,(4,2^(N-2)))
        HPsi = reshape(HPsi,(4,2^(N-2)))
        HPsi += H2*Psi                 
        Psi = reshape(Psi,(2,2^(N-1)))
        HPsi = reshape(HPsi,(2,2^(N-1)))
        Psi = permutedims(Psi,(2,1))
        HPsi = permutedims(HPsi,(2,1))
    end
    Psi = reshape(Psi,L)
    HPsi = reshape(HPsi, L)
    return HPsi
end

N=2*5
PsiInitial = separable(N);

nA = Int(N/2) # number spins/qubits in A
nB = Int(N/2) # number spins/qubits in B

PsiAB = reshape(PsiInitial, (2^(nA),2^(nB)))

rdmA = compute_rdmA(PsiAB)

isapprox(compute_entropy(rdmA), 0, atol=10^-14)

theta = pi/4 # critical magnetic field
H2 = buildIsing(theta) # Hamiltonian
D,U = eigen(H2)
shiftE = D[end] 
H2 = H2 - shiftE*diagm(0=>ones(4))


init_steps = 1
entropy = zeros(1)
normdiff = zeros(1)
initial_step=1
Nsteps = 0;

N = 2 * 5
Psi = separable(N);
deltat = 0.0001

initial_step += Nsteps
Nsteps = 80000  # change this number if the entanglement entropy has not reached its plateau 
final_step = initial_step + Nsteps-1


for n=initial_step:final_step
    Psi = Psi .- im * deltat .* multiplyHPsi(Psi,H2)
    if n%400==0
    normdiff = [normdiff; abs(1 - sqrt(abs(Psi'*Psi)))]
    #Psi = normalizeVec(Psinew)
    PsiAB = reshape(Psi, (2^(Int(N/2)),2^(Int(N/2))))
    rdmA = compute_rdmA(PsiAB)
    entropynew = compute_entropy(rdmA)
    entropy = [entropy; entropynew]
    end
    if n%4000==0 print(n-initial_step+1, ":", Nsteps, " ") end
end


tarray = collect(range(0, stop = deltat * Nsteps, length = length(normdiff)))
standard = [10^-5 for i in 1:length(tarray)]
plot(tarray, standard, label = "goal", xlabel = "t", ylabel = "normalization loss")
plot!(tarray, normdiff, label = "numerical", markershape = :hexagon, markerstrokealpha = 0,
    markersize = 2)