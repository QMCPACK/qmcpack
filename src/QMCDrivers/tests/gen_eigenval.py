import numpy as np
import scipy.linalg

# Compute values for test_LinearMethod.cpp, the solveGeneralizedEigenvalues test case.

# Overlap
S = [[1.0, 0.1], [0.15, 1.0]]

# Hamiltonian
H = [[-4.0, 0.2], [0.3, 1.0]]


def compute_eigen_generalized(S, H):
    print("Eigenvalues from generalized eigenvalue problem for S and H")
    out = scipy.linalg.eig(H, S, left=False, right=True)
    evals = out[0]
    r_evec = out[1]
    print("Eigenvalues")
    print(evals)
    print("Raw Eigenvectors")
    print(r_evec)
    print("Scaled Eigenvectors")
    print(r_evec[:, 0] / r_evec[0, 0])
    print(r_evec[:, 1] / r_evec[1, 1])


def compute_eigen_inv(S, H):
    print("Eigenvalues from S^-1*H")
    r2 = np.dot(np.linalg.inv(S), H)
    out = scipy.linalg.eig(r2)
    evals = out[0]
    r_evec = out[1]
    print("Eigenvalues")
    print(evals)
    print("Raw Eigenvectors")
    print(r_evec)
    print("Scaled Eigenvectors")
    print(r_evec[:, 0] / r_evec[0, 0])
    print(r_evec[:, 1] / r_evec[1, 1])


if __name__ == "__main__":
    compute_eigen_generalized(S, H)
    # compute_eigen_inv(S,H)
