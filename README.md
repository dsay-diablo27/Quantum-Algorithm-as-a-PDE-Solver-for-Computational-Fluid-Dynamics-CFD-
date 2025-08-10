# Quantum Algorithm as a PDE Solver for Computational Fluid Dynamics (CFD)

Welcome to the submission for **Project 3: Quantum Algorithm as a PDE Solver for CFD** under the [WISER Quantum Program 2025](https://www.thewiser.org/about-wiser). This repository presents a  quantum-algorithm approach to solving the **Burgers' Equation** via quantum algorithms.

---

## 🚀 Project Overview

Classical Computational Fluid Dynamics (CFD) solvers face increasing demands for resolution and computational power, especially for nonlinear partial differential equations (PDEs). Quantum algorithms, co-designed with classical techniques, offer state compression and speed-ups that could revolutionize CFD.

This project is an open challenge (in partnership with [BQP](https://www.bqpsim.com/)) to implement and benchmark quantum-enhanced PDE solvers for the 1D **Burgers' Equation**:

$$
\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} = \nu \frac{\partial^2 u}{\partial x^2}
$$

**Boundary conditions:**
- u(0, t) = U_L  (left),  u(1, t) = U_R  (right)

**Initial condition:**
- \( u(x, 0) = 1 \) for \( x < 0.5 \), otherwise \( u(x, 0) = 0 \)

---

## 💡 Motivation

- The **Burgers’ equation** serves as an ideal benchmark: retains the nonlinear and viscous core of Navier–Stokes, remains analytically tractable, and is a proving ground for quantum solvers.
- Quantum-classical hybrid methods validated here can be translated to full Navier–Stokes workloads. So, its Important to first Simulate Burger's Equation.

---

## 🧠 Core Algorithms Used:
- **Hydrodynamic Schrödinger Equation (HSE):** Recasts incompressible flow as quantum wavefunction dynamics ([Meng & Yang, 2023](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.033182)).


---

## 📝 Project Deliverables
- Two Jupyter notebook one covers the classical approach and other covers the Quantum Transformation based Algorithm, the framework uses the Research paper mentioned above.
- A Report file approx 5 page with analysis.
---


Note: The classical solution uses the Cole–Hopf transformation.
