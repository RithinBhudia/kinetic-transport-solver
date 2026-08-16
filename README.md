# Kinetic Transport Solver

A numerical solver written in C for a two-dimensional kinetic transport problem in position–velocity space.

The project explores finite-difference methods for advection and diffusion, with configurable initial conditions, boundary conditions and physical parameters.

## Features

* First-order upwind discretisation for transport in position and velocity
* CFL-based automatic timestep selection
* Implicit finite-difference treatment of velocity diffusion
* Tridiagonal linear systems solved using LAPACK
* Configurable Gaussian initial conditions and inflow boundary conditions
* Configurable spatially varying force
* Simulation parameters loaded from a configuration file
* Results written to CSV for further analysis and visualisation

## Project Structure

```text
kinetic-transport-solver/
├── src/
│   ├── main.c
│   ├── config.c
│   ├── solver.c
│   └── kinetic.h
├── examples/
│   └── beam.cfg
├── scripts/
│   └── plot_solution.py
├── output/
├── Makefile
├── README.md
└── .gitignore
```

## Numerical Method

The distribution function is represented on a two-dimensional grid in position (x) and velocity (v).

Transport in both dimensions is approximated using a first-order upwind finite-difference scheme. The direction of the difference is selected according to the sign of the corresponding transport velocity.

The timestep is automatically chosen using a CFL-type stability restriction based on the maximum velocity and force.

Velocity-space diffusion is treated implicitly. At each spatial grid point this produces a tridiagonal linear system, which is solved using the LAPACK `DGTSV` routine.

## Configuration

Simulation parameters are supplied through a configuration file using

```text
key = value
```

syntax.

Parameters include:

* number of spatial and velocity grid points
* spatial and velocity domain limits
* final simulation time
* CFL number
* diffusion coefficient
* force parameters
* Gaussian initial-condition parameters
* left and right inflow conditions
* output filename

An example configuration is provided in:

```text
examples/beam.cfg
```

## Building

The project requires:

* a C compiler such as GCC
* BLAS
* LAPACK
* `make`

Compile with:

```bash
make
```

This creates the executable:

```text
solver
```

To remove generated build files:

```bash
make clean
```

## Running

Run using the example configuration:

```bash
./solver examples/beam.cfg
```

If no configuration file is supplied, the program uses the default example configuration.

Simulation results are written as CSV data containing:

```text
x,v,f
```

where `f` is the computed distribution function at each point in position–velocity space.

## Motivation

This project was developed to explore numerical methods used in scientific computing and to practise implementing finite-difference algorithms efficiently in C.

It demonstrates the use of explicit and implicit numerical methods, stability constraints, dynamic memory management, configuration parsing and external numerical libraries.
