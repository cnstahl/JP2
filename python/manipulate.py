import numpy as np

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
    r = np.zeros((M,M))
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

def entropy(M):
    w, v = np.linalg.eig(M)
    S = 0
    for i in w:
        if (i != 0): S += i*np.log(i)
    return -S

def density_mat(psi):
    return np.outer(psi, psi.conj())