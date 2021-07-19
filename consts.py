import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import mpmath as mpm
#import scipy.linalg as slin

# total lattice number (each lattice constains 2 sublattices A and B)
N=2**10

# Gaussian wavepacket parameters
# center
xc=N

#width
sgm=30

#parameters for coefficients v, w, u
D0=2
d0=0.8
J=-1
#parameters of linear part of Hamiltonian
omega=0.03
omegaF=0.05

T=2*np.pi/omega
Q=2**14
tTot=3*T
dt=tTot/Q
#coef functions
def u(t):
    return D0 * np.cos(omega * t)


def v(t):
    return (J + d0 * np.sin(omega * t))


def w(t):
    return (J - d0 * np.sin(omega * t))

#init Gaussian part of the wavefunction
wvFcnt = [np.exp(-(n - xc) ** 2 / (4 * sgm ** 2)) for n in range(0, 2 * N)]
# bt=0.1
#
# wvFcnt=[float(mpm.sech(bt*(n-xc))**2) for n in range(0,2*N)]

# lower band
t0=0
k0=0
AVal=(v(t0)+w(t0))*np.cos(k0)
BVal=(v(t0)-w(t0))*np.sin(k0)
CVal=u(t0)
EVAL=np.sqrt(AVal**2+BVal**2+CVal**2)
# v1=[iB-A, C+E] without normalization const
u1A=1j*BVal-AVal
u1B=CVal+EVAL
for j in range(0,N):
    wvFcnt[2*j]*=u1A
    wvFcnt[2*j+1]*=u1B
#FT coef
for n in range(0,2*N):
    wvFcnt[n]*=np.exp(-1j*k0*n)


# # normalization const
c02Tmp = 0
for elem in wvFcnt:
    c02Tmp += np.abs(elem) ** 2
C0 = np.sqrt(c02Tmp)

psi0=[]
for n in range(0,2*N):
    psi0.append(wvFcnt[n]/C0)



