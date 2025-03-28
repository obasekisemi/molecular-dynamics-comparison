# Molecular Dynamics Simulation Validation

This repository contains implementations of Lennard-Jones molecular dynamics simulations in both Python and LAMMPS. The purpose is to validate the custom Python implementation against the established LAMMPS package.

## Contents

- **python-md/**: Custom Python implementation of molecular dynamics for Lennard-Jones particles
- **lammps-md/**: LAMMPS input scripts for equivalent simulation

## Validation Results

Comparison of the two implementations shows excellent agreement:
- Temperature: < 1% difference
- Total Energy: ~2% difference
- Potential Energy: ~15% difference
- Pressure: ~14% difference

The Radial Distribution Functions demonstrate that both implementations correctly capture the structural properties of a Lennard-Jones fluid, confirming the validity of the custom Python implementation.

## Units

All simulations use reduced Lennard-Jones units:
- Energy: ε = 1.0
- Length: σ = 1.0
- Pressure: ε/σ³

## Parameters

- 100 particles
- Cubic box of side length 10σ
- Cutoff distance: 2.5σ
- Timestep: 1.0E-4
- 1000 simulation steps
- Target temperature: 300
