import numpy as np

def set_zeros(w):
    tol = 1e-14
    w.real[abs(w.real) < tol] = 0.0
    if np.iscomplexobj(w): w.imag[abs(w.imag) < tol] = 0.0
    if np.all(w.imag == 0): return w.real
    return w

def partial_trace(M, psi):
    N = len(psi)-M
    rho = np.zeros((M,M))
    for i in range(M):
        for j in range(M):
            rho[i,j] = sum(psi[i:i+N+1]*psi[j:j+N+1].conjugate())
    return rho/M

def eigen(H):
    w, v = np.linalg.eigh(H)
    #w = set_zeros(w)
    #v = set_zeros(v)
    return w, v
    p = w.argsort()
    return w[p], v.T[p].T

def reconstruct(w,v):
    return np.dot(np.dot(v, np.diag(w)), v.conj().T)

def entropy(M):
    w, v = np.linalg.eig(M)
    rho = 0
    for i in w:
        if (i != 0): rho += i*np.log(i)
    return rho