from __future__ import print_function, division
import sys,os
quspin_path = os.path.join(os.getcwd(),"../../")
sys.path.insert(0,quspin_path)
from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis import spin_basis_1d # Hilbert space spin basis
from quspin.tools.measurements import *
import numpy as np
import matplotlib.pyplot as plt

L=12 # size
#coupling strength
J=1.0 #spin-spin coupling
h=0.8945 
g=0.945 
#create site-coupling lists 
J_zz=[[J,i,(i+1)%L] for i in range(L)] # PBC
x_field=[[h,i] for i in range(L)]
z_field=[[g,i] for i in range(L)]
#static and dynamic lists
static_1=[["x",x_field],["z",z_field]]
static_2=[["zz",J_zz],["x",x_field],["z",z_field]]

dynamic=[]

basis = spin_basis_1d(L=L, kblock=0, pblock=1)

H1=hamiltonian(static_1,dynamic,basis=basis,dtype=np.float64)
H2=hamiltonian(static_2,dynamic,basis=basis,dtype=np.float64)

E1,V1=H1.eigh()
psi1=V1[:,14] # pick any state as initial state
E2,V2=H2.eigh()


Sent=ent_entropy(psi1,basis,chain_subsys=[1,3,6,7,11])
print(Sent['Sent_A'])

