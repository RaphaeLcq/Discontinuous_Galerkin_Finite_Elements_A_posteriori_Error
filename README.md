# A Posteriori Error Estimation for Adaptive Mesh Refinement in Discontinuous Galerkin Methods

## Overview
This project implements an **a posteriori error estimator** for **Discontinuous Galerkin (DG)** methods, with the goal of enabling **adaptive mesh refinement (AMR)** for the diffusion operator. It is based on the Symmetric Interior Penalty method.

The work is based on the estimator proposed by Ainsworth (2007), extended to handle **heterogeneous diffusion problems**, and validated through numerical experiments.

This is a research prototype implemented with Octave GNU.

The PDF report sums up the theory used and the numerical results.
---

## Motivation
In numerical PDEs, **a priori estimates** predict convergence rates, but they do not provide the actual error of a computed solution.

In practice:
- the exact solution is unknown
- error is **non-uniformly distributed**
- uniform mesh refinement is inefficient

**A posteriori error estimation** addresses this by:
- estimating the error from the computed solution itself
- identifying regions where the solution is inaccurate
- guiding **local mesh refinement**

---

## Problem Setting

We consider diffusion problems of the form:

- Diffusion PDEs (e.g. Poisson)
- Possibly **heterogeneous diffusion coefficients**
- Solved using **Discontinuous Galerkin (DG)** SIP or SWIP methods

The numerical error is measured in an energy norm (broken flux norm), denoted:

- global error: error_h
- local error: error_T on each mesh element T  

---

## Method

### A Posteriori Error Estimator

An estimator eta_h is constructed such that:

- it depends only on:
  - the computed solution
  - the mesh
  - the problem data
- it is as much **equivalent to the true error asymptotically** as possible

This ensures:
- reliability: eta_h ≥ error_h  
- efficiency: eta_h ≈ error_h as h → 0  

---

### Error Decomposition

The error is split into two components:

- **Conforming error (eta_CF)**  
- **Non-conforming error (eta_NC)**  

This reflects the DG structure:
- discontinuities across elements introduce a non-conforming contribution  
- approximation inside elements contributes to conforming error

---

### Adaptive Mesh Refinement (AMR)

The estimator is used to:

1. Compute **local error indicators** eta_T
2. Identify elements with largest error  
3. Refine the mesh **locally**  

Goal:
> concentrate degrees of freedom where they are needed

This leads to:
- improved accuracy for fixed computational cost  
- better handling of singularities and irregular solutions

---

## Numerical Results

### Key Observations

- The estimator shows **good effectivity** (eta_h ≈ error_h)  
- Local error maps correctly identify high-error regions  
- Adaptive refinement reduces error more efficiently than uniform refinement  

For smooth problems:
- improvement is moderate (expected)

For singular problems (e.g. L-shaped domain):
- **significant gain in convergence rate**  
- local refinement captures singularities effectively 

---

## Applications

- Adaptive finite element methods (hp-AMR)
- Industrial simulations
- Multiphysics and heterogeneous media problems

---

## Implementation

- Prototype developed in:
  - MATLAB, Octave

- Designed for future integration into **TrioCFD (CEA)** :

---

## Repository Structure
- DEBUG : debugging files, used for sanity check
- EstimPosterioriAinsworth : error estimation files
    - EstimationConforme : Conform error estimation files
    - EstimationNonConforme : Non conform error estimation files
- FonctionsProblemes : Problem functions for different benchmarks
    - HarmonicSquare : Files for Delta u = 0, u = r^alpha sin(alpha theta) on the border (alpha is an integer)
    - Lshape : Files for Delta u = 0, u = r^alpha sin(alpha theta) on the border (alpha = 2/3)
    - LshapeNeumann : Files for Delta u = 0, grad_n u = grad_n (r^alpha sin(alpha_theta)) on the border, (alpha = 2/3)
    - SquareNeumann : Files for Delta u =  sin(lambda * pi .* x) .* sin(lambda * pi .* y), grad_n u =  grad_n(sin(lambda * pi .* x) .* sin(lambda * pi .* y)) on the border (lambda is an integer)
    - SquarePoisson : Files for Delta u =  sin(lambda * pi .* x) .* sin(lambda * pi .* y), grad_n u = 0 on the border (lambda is an integer)
    - SWIP_Square : Files for div( kappa grad u ) =  0, u = r^alpha( a cos(alpha theta) + b sin(alpha theta)) on the border (alpha > 0), kappa = kappa(T).
    - TopNeumann : Files for TP3 of https://ecoles-cea-edf-inria.fr/en/schools/ecole-analyse-numerique-2022/.
- MatGlo : Global matrices files
- MatGloGradient_SWIP : Global matrices files used for SWIP method
- MatLoc : Local matrices files
- MeshUtils : Meshing files
- Poisson : Files embedding the different computations for the solution
- Solution : Files used to compute to build the system Ax = b using different methods
- SWIP : Files used for the Symmetric Weighted Interior Penalty method


---

## How to Use a ready to run problem

0. Create a .msh file using the .geo files already created.
Use GMESH export : Version 2, ASCII, checking "Save all elements" only. 
The name of the geo file is : "ProblemName.geo".
The .msh files should be named : ProblemName_h{X}.msh, where X = 1, ..., 5. X = 1 being the first non refined file, X > 1 are a sequence refined meshes.
Note the program expects that the greater X is, the more the mesh is refined.

1. Run TrouveChemin.m to add the project folder in the path of Octave.

2.  Open mainPoissonDG.m and select the problem you want to run.

3. Choose what vizualisation you want the the "GLOBAL visu".

4. Choose the number of simulations you want to run. 
The simulations will run for X = m0 to X = m1 if m0 > m1.

5. Refine mesh based on eta_T.

5. Repeat.

---

## References

- M. Ainsworth (2007) — A posteriori error estimation for DG methods  
- Di Pietro & Ern (2012) — DG theory  
- Ern & Vohralík — flux-based estimators  

---

## Author

Raphaël Lecoq

Erell Jamelot

Andrew Peitavy

Thanks to Melissa Mroueh for her help and some parts of the script

---

## Notes
This repository is a **research prototype**, intended for:
- numerical analysis validation
- methodological exploration

Not production-ready.

# Copyright 

Copyright (c), CEA
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
