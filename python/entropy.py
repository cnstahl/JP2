import math
import sys 
import numpy as np
import random
from cmath import rect

def anti_com(a,b):
    return np.dot(a,b)+np.dot(b,a)

def Fermi(n,j):
    F = 0
    for i in np.arange(0,j):
        if (n&2**i) != 0: F += 1
    return F

def psi_i(N,j):
    assert N%1 == 0, "N must be an integer"
    assert j%1 == 0, "j must be an integer"
    assert j<N, "j must be less than N"
    assert N>=0, "N must be greater than or equal to 0"
    assert j>=0, "j must be greater than or equal to 0"
    k = np.floor(j/2)
    psi = np.zeros((2**N,2**N))
    for n, row in enumerate(psi):
        for m, elem in enumerate(row):
            if (n-m == 2**j) and (n&2**j != 0) and (m&2**j == 0): 
                psi[m,n] = (-1)**Fermi(n,j)
    return psi

def psi_bar_i(N,j):
    assert N%1 == 0, "N must be an integer"
    assert j%1 == 0, "j must be an integer"
    assert j<N, "j must be less than N"
    assert j>=0, "j must be greater than or equal to 0"
    psi = np.zeros((2**N,2**N))
    for m, row in enumerate(psi):
        for n, elem in enumerate(row):
            if   (m-n == 2**j) and (n&2**j == 0) and (m&2**j != 0): 
                psi[m,n] = (-1)**Fermi(n,j)
    return psi

def psis(N):
    psi = np.zeros((N,2**N,2**N))
    for i in range(N):
        psi[i] = psi_i(N,i)
    return psi

def psi_bars(N):
    psi = np.zeros((N,2**N,2**N))
    for i in range(N):
        psi[i] = psi_bar_i(N,i)
    return psi

def check_psis(N):
    assert N%1 == 0, "N must be an integer"
    assert N>=0, "N must be greater than or equal to 0"
    psi = psis(N)
    psi_bar = psi_bars(N)
    for i in np.arange(N):
        for j in np.arange(N):
            acom = anti_com(psi[j], psi[i])
            if np.any(acom): return(0, "{psi(%d,%d), psi(%d,%d)} \ne 0" % (N,j,N,i))
            acom = anti_com(psi_bar[i], psi_bar[j])
            if np.any(acom): return(0, "{psi_bar(%d,%d), psi_bar(%d,%d)} \ne 0" % (N,j,N,i))
            acom = anti_com(psi[j], psi_bar[i])
            if i == j:
                if np.any(acom - np.eye(2**N)): return(0, "{psi(%d,%d), psi_bar(%d,%d)} \ne 1" % (N,j,N,i))
            else: 
                if np.any(acom): return(0, "{psi(%d,%d), psi_bar(%d,%d)} \ne 0" % (N,j,N,i))
    return(1, "good")

def product(*args):
    prod = np.eye(np.shape(args[0])[0])
    for arg in args:
        prod = np.dot(prod,arg)
    return prod

def get_C(N,J):
    len = N*(N-1)*(N-2)
    Cs = np.zeros((N,N,N),dtype='complex128')
    for k in range(N):
        for j in range(k):
            for i in range(j):
                mod = np.sqrt(abs(random.gauss(0,np.sqrt(2*J)/N)))
                phase = 2*np.pi*random.random()
                Cs[i,j,k] = rect(mod,  phase)
    return Cs

def hamiltonian(N, J=1, C=None):
    Q    = np.zeros((2**N, 2**N),dtype='complex128')
    Qbar = np.zeros((2**N, 2**N),dtype='complex128')
    if C is None: C = get_C(N,J)
    # print(type(C))
    psi = psis(N)
    psi_bar = psi_bars(N)
    #  we want 0\le i<j<k<N
    for k in range(N):
        for j in range(k):
            for i in range(j):
                # mod = np.sqrt(abs(random.gauss(0,np.sqrt(2*J)/N)))
                # phase = 2*np.pi*random.random()
                Q +=                 C[i,j,k] *product(psi[i],     psi[j],     psi[k])
                Qbar += np.conjugate(C[i,j,k])*product(psi_bar[i], psi_bar[j], psi_bar[k])
    return -anti_com(Q, Qbar)

def damiltonian(N, J):
    Q    = np.zeros((2**N, 2**N),dtype='complex128')
    Qbar = np.zeros((2**N, 2**N),dtype='complex128')
    psi = psis(N)
    psi_bar = psi_bars(N)
    #  we want 0\le i<j<k<N
    for k in range(N):
        for j in range(k):
            for i in range(j):
                Q +=    product(psi[i],     psi[j],     psi[k])
                Qbar += product(psi_bar[i], psi_bar[j], psi_bar[k])
    return anti_com(Q, Qbar)


def set_zeros(w):
    tol = 1e-10
    w.real[abs(w.real) < tol] = 0.0
    if np.iscomplexobj(w): w.imag[abs(w.imag) < tol] = 0.0
    if np.all(w.imag == 0): return w.real
    return w

def single_trace(rho):
    N = len(rho)
    assert N%2 == 0, "\psi has to have even length: %s" % N
    M = int(N/2)
    r = np.zeros((M,M), dtype='complex')
    for idx in range(M):
        for jdx in range(M):
            r[idx, jdx] = rho[idx, jdx] + rho[idx + N/2, jdx + N/2]
    return r

def partial_trace(rho, traces=1):
    N = len(rho)
    assert N >= 2**traces, "Too many particles to trace over"
    for i in range(traces):
        rho = single_trace(rho)
    return rho

def eigen(H):
    w, v = np.linalg.eigh(H)
    #w = set_zeros(w)
    #v = set_zeros(v)
    return w, v
    p = w.argsort()
    return w[p], v.T[p].T

def reconstruct(w,v):
    return np.dot(np.dot(v, np.diag(w)), v.conj().T)

def fentropy(M):
    w, v = np.linalg.eig(M)
    S = 0
    set_zeros(w)
    assert np.all(np.imag(w) == 0)
    for i in np.real(w):
        if (i != 0): S += i*np.log(i)
    return -S

def density_mat(psi):
    return np.outer(psi, psi.conj())


N = 10
entropy = np.zeros(N+1)
I = 10
for i in range(I):
    print(i, "start")
    H = hamiltonian(N,1)
    w, v = eigen(H)
    w = set_zeros(w)
    ngs = sum(w==0)
    tot_gs = np.copy(v[:,0])
    for i in range(1,ngs):
        tot_gs += v[:,i]
    tot_gs = tot_gs/math.sqrt(ngs)
    rho = density_mat(tot_gs)
    entropies = np.zeros(N+1)
    rhoi = np.copy(rho)
    for idx, val in enumerate(entropies):
        if idx != 0: rhoi = single_trace(rhoi)
        entropies[idx] = fentropy(rhoi)
    entropy = entropy + entropies
    print(i, "end")
np.savetxt("../data/N%savg_entropy.txt" % N, np.transpose(entropy/I))
