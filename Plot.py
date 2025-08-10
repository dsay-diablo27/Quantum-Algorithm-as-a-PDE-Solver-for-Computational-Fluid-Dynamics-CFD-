import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt

# Simulate Burgers' equation using Cole-Hopf transformation and Crank-Nicolson scheme
def simulate_burgers_cole_hopf(nu, L, Nx, T, dt, snapshot_times=[]):
    dx = L / (Nx - 1)
    x = np.linspace(0, L, Nx)
    
    # Initial condition: Riemann step (u=1 for x<=0.5, u=0 otherwise)
    u0 = np.where(x <= 0.5, 1.0, 0.0)

    # Compute integral for Cole-Hopf initial transformation
    integral = np.zeros(Nx)
    for i in range(1, Nx):
        # Trapezoidal rule for integration
        integral[i] = integral[i-1] + 0.5 * (u0[i] + u0[i-1]) * dx
    phi = np.exp(-integral / (2 * nu))
    
    # Crank-Nicolson matrices for time-stepping the heat equation
    r = nu * dt / dx**2
    A = np.eye(Nx)
    B = np.eye(Nx)
    for i in range(1, Nx-1):
        # Fill in tridiagonal entries for interior points
        A[i, i-1] = -r/2
        A[i, i]   = 1 + r
        A[i, i+1] = -r/2
        B[i, i-1] = r/2
        B[i, i]   = 1 - r
        B[i, i+1] = r/2
    # Neumann boundary conditions for phi (zero derivative at boundaries)
    A[0,0] = A[-1,-1] = 1
    A[0,1] = A[-1,-2] = 0
    B[0,0] = B[-1,-1] = 1
    B[0,1] = B[-1,-2] = 0
    
    Nt = int(T/dt)
    times = np.arange(0, T+dt/2, dt)
    snapshot_indices = [int(round(t/dt)) for t in snapshot_times]
    
    phi_hist = []
    phi_hist.append(phi.copy())
    for n in range(Nt):
        b = B @ phi
        # Apply boundary conditions for phi at each step
        b[0] = phi[0]
        b[-1] = phi[-1]
        phi_new = solve(A, b)
        phi = phi_new
        # Save snapshots at requested times
        if (n+1) in snapshot_indices:
            phi_hist.append(phi.copy())
    
    u_hist = []

    for idx, phi_snap in enumerate(phi_hist):
        # Recover u from phi using Cole-Hopf formula
        lnphi = np.log(np.maximum(phi_snap, 1e-14))
        dlnphi_dx = np.zeros_like(lnphi)
        # Use central differences for interior, one-sided for boundaries
        dlnphi_dx[1:-1] = (lnphi[2:] - lnphi[:-2])/(2*dx)
        dlnphi_dx[0] = (-3*lnphi[0] + 4*lnphi[1] - lnphi[2])/(2*dx)
        dlnphi_dx[-1] = (3*lnphi[-1] - 4*lnphi[-2] + lnphi[-3])/(2*dx)
        u_hist.append(-2*nu*dlnphi_dx)
    return x, snapshot_times, np.array(u_hist)
