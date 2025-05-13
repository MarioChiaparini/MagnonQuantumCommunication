from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.operators import exp_op
from quspin.tools.measurements import ent_entropy

import numpy as np
import os
import pickle as cPickle

import sys


def Hamiltonian(L, J, hz, fun=None, fun_args=[]):
    basis = spin_basis_1d(L=L, kblock=0, pblock=1, pauli=False)
    zz_int = [[-J,i,(i+1)%L] for i in range(L)] # PBC
    x_field=[[-1.0,i] for i in range(L)]
    z_field = [[-hz,i] for i in range(L)]

    static = [["zz",zz_int],["z",z_field]]
    dynamic = [["x",x_field,fun,fun_args]]

    kwargs = {'dtype':np.float64,
              'basis':basis,
              'check_symm':False,
              'check_herm':False,
              'check_pcon':False}
    H = hamiltonian(static,dynamic,**kwargs)
    return H

def Unitaries(delta_time,L,J,hz,action_min,var_max,var_min,state_i,save=False,save_str=''):

    b = 0
    lin_fun = lambda t: b
    H = hamiltonian(L=L, fun=lin_fun, **{'J': J, 'hz': hz})

    n = int((var_max - var_min) / action_min)
    for i in range(n+1):
        b = state_i[0]+1*action_min
        expm_dict = np.asarray([exp_op(H,)]) 
    
    b =+2.0
    _,V_target = H.eigh()
    
    if save:
        
        str1=os.getcwd()
        str2=str1.split('\\')
        n=len(str2)
        my_dir = str2[n-1]
		
        # create directory if non-existant
		
        save_dir = my_dir+"/unitaries"
		
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        save_dir="unitaries/"
        
        dataname  = save_dir + "unitaries_L={}".format(L)+save_str+'.pkl'
        cPickle.dump(expm_dict, open(dataname, "wb" ) )
        
        dataname  = save_dir + "target_basis_L={}".format(L)+save_str+'.pkl'
        cPickle.dump(V_target, open(dataname, "wb" ) )
    else:
        return expm_dict



def MB_observables(psi,times,protocol,pos_actions,h_field,L,J=1.0,hx_i=-2.0,hx_f=2.0,hz=1.0,fin_vals=False,bang=True):
    
    str1=os.getcwd()
    str2=str1.split('\\')
    n=len(str2)
    my_dir = str2[n-1]

	# define Hamiltonian
    b=hx_f  
    lin_fun = lambda t: b
    # define Hamiltonian
    H = Hamiltonian(L,fun=lin_fun,**{'J':J,'hz':hz})
    if bang:
        exp_dict_dataname = my_dir+"/unitaries/unitaries_L={}_bang".format(L)+'.pkl'
        V_target_dataname = my_dir+"/unitaries/target_basis_L={}_bang".format(L)+'.pkl'
    else:
        exp_dict_dataname = my_dir+"/unitaries/unitaries_L={}_cont".format(L)+'.pkl'
        V_target_dataname = my_dir+"/unitaries/target_basis_L={}_cont".format(L)+'.pkl'
    
    subsys=[i for i in range(L//2)]

    # preallocate variables
    Fidelity,E,delta_E,Sd,Sent=[],[],[],[],[]

    expm_dict=cPickle.load(open( exp_dict_dataname, "rb" ))
    V_target=cPickle.load(open( V_target_dataname, "rb" ))

    i=0
    while True:
        Fidelity.append( abs(psi.conj().dot(V_target[:,0]))**2 )
        EGS = H.eigsh(k=1,which='SA',maxiter=1E10,return_eigenvectors=False)[0]
        E.append( H.matrix_ele(psi,psi).real/L - EGS/L)
        delta_E.append( np.sqrt( (H(time=0)*H).matrix_ele(psi,psi) - H.matrix_ele(psi,psi)**2 + 1E-21*1j).real/L   )
        pn = abs( V_target.conj().T.dot(psi) )**2.0 + np.finfo(psi[0].dtype).eps
        Sd.append( -pn.dot(np.log(pn))/L )
        # entanglement entropy
        if L==1:
            Sent.append(0.0)
        else:
            Sent.append( ent_entropy(psi,H.basis,chain_subsys=subsys)['Sent'] )
        if i == len(protocol):
            break
        else:
            b=protocol[i] # --> induces a change in H
            #psi = exp_H.dot(psi)
            psi = expm_dict[int(np.rint((b - min(h_field))/min(pos_actions)))].dot(psi)
            i+=1


    if fin_vals:
        return Fidelity[-1],E[-1],delta_E[-1],Sd[-1],Sent[-1]
    else:
        return Fidelity,E,delta_E,Sd,Sent
    
sys.stdout.flush()
str1=os.getcwd()
str2=str1.split('\\')
n=len(str2)
my_dir = str2[n-1]

max_t_steps_vec=np.array([4,8,10,12,20,30,40,48,50,60,70,80]) 


L = int(sys.argv[4]) # system size

if L==1:
	J=0
else:
	J = 1.0 # zz interaction

hz = 1.0 # hz field
hx_i= -2.0 # initial hx coupling
hx_f= +2.0 # final hx coupling

b=hx_i
lin_fun = lambda t: b

H_params = {'J':J,'hz':hz}
H = Hamiltonian.Hamiltonian(L,fun=lin_fun,**H_params)