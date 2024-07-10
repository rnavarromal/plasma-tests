import numpy as np
import matplotlib.pyplot as plt


def potential_relaxation_method_loop(f, phi, h, tol=1e-5, maxiter=1000):
    """Solves the poisson equation `Del phi = f` through the relaxation
    method. The potential phi[i,j] at the node [i,j] is calculated
    explicitly within a for-loop."""
    phi = np.copy(phi)

    for iteration in range(maxiter):
        max_diff = 0.0
        for i in range(1, N-1):
            for j in range(1, N-1):
                phi_new = 0.25 * (phi[i-1, j] + phi[i+1, j] + phi[i, j-1] + phi[i, j+1] - h**2 * f[i, j])

                max_diff = max(max_diff, np.abs(phi[i,j]-phi_new))
                phi[i,j] = phi_new

        if max_diff < tol:
            print(f"Converged after {iteration+1} iterations")
            break
    else:
        print("Did not converge within the maximum number of iterations")

    return phi


def potential_relaxation_method_vectorized(f, phi, h, tol=1e-5, maxiter=1000):
    """Solves the poisson equation `Del phi = f` through the relaxation
    method. The potential phi is calculated through vectorization of the relaxation method."""
    phi = np.copy(phi)

    for iteration in range(maxiter):
        phi_new = 0.25 * (phi[:-2, 1:-1] + phi[2:,1:-1] + phi[1:-1, :-2] + phi[1:-1, 2:] - h**2 * f[1:-1, 1:-1])

        # Check for convergence
        max_diff = np.max(np.abs(phi[1:-1, 1:-1] - phi_new))
        phi[1:-1, 1:-1] = phi_new

        if max_diff < tol:
            print(f"Converged after {iteration+1} iterations")
            break
    else:
        print("Did not converge within the maximum number of iterations")

    return phi




def potential_relaxation_method_SOR(f, phi, h, w=1.5, tol=1e-5, maxiter=1000):
    """Solves the poisson equation `Del phi = f` through a succesive
    over-relaxation method. The potential phi is calculated through an
    explicit for-loop.
    """

    phi = np.copy(phi)

    for iteration in range(maxiter):
        max_diff = 0.0
        for i in range(1, N-1):
            for j in range(1, N-1):
                phi_new = (1-w)*phi[i,j] + w * 0.25 * (phi[i-1, j] + phi[i+1, j] + phi[i, j-1] + phi[i, j+1] - h**2 * f[i, j])

                max_diff = max(max_diff, np.abs(phi[i,j]-phi_new))
                phi[i,j] = phi_new

        if max_diff < tol:
            print(f"Converged after {iteration+1} iterations")
            break
    else:
        print("Did not converge within the maximum number of iterations")

    return phi


def potential_relaxation_method_SOR_vectorized(f, phi, h, w=1.5, tol=1e-5, maxiter=1000):
    """Solves the poisson equation `Del phi = f` through a succesive
    over-relaxation method. The potential phi is calculated through vectorization of the relaxation method."""
    phi = np.copy(phi)

    for iteration in range(maxiter):
        phi_new = (1-w)*phi[1:-1, 1:-1] + w * 0.25 * (phi[:-2, 1:-1] + phi[2:,1:-1] + phi[1:-1, :-2] + phi[1:-1, 2:] - h**2 * f[1:-1, 1:-1])

        # Check for convergence
        max_diff = np.max(np.abs(phi[1:-1, 1:-1] - phi_new))
        phi[1:-1, 1:-1] = phi_new

        if max_diff < tol:
            print(f"Converged after {iteration+1} iterations")
            break
    else:
        print("Did not converge within the maximum number of iterations")

    return phi


N = 50  # Number of grid points
h = 1.0 / (N - 1)  # Grid spacing
phi = np.zeros((N, N))  # Potential array initialized to zero
f = np.zeros((N, N))  # Source term array

# Example source term: point source at the center
f[N//2, N//2] = 100.0

# Solve the Poisson equation
phi = potential_relaxation_method_vectorized(f, phi, h, tol=1e-5)

# Visualize the solution

plt.imshow(phi, extent=[0, 1, 0, 1], origin='lower', cmap='viridis')
plt.colorbar(label='Potential $\phi$')
plt.title('Solution of Poisson Equation using Relaxation Method')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


plt.plot(phi, 'o')
plt.show()
