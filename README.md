# Numerical Analysis for Partial Differential Equations

This repository contains **Python implementations of numerical methods for solving Partial Differential Equations (PDEs)**. The scripts were developed as part of the course **Numerical Analysis for Partial Differential Equations** at **TU Delft**.

The main objective of this repository is to demonstrate the **numerical discretisation and solution of PDEs**.

---

## Repository Structure

```text
Numerical-Analysis-for-PDE/
├── assignment1/
│   ├── Problem2/
│   ├── Problem3/
│   └── Problem4/
├── assignment2/
│   ├── Problem1/
│   ├── Problem2/
│   ├── Problem3/
│   ├── Problem4/
│   └── Problem5/
├── utils/
│   ├── SparseMatrices.py
│   └── meshGrid.py
├── assignments.py
├── main.py
├── __init__.py
├── README.md
└── LICENSE
```

## Contents
Assignments

Each assignment folder contains individual Python scripts for each problem defined in the course assignments.

    - assignment1/
    Contains numerical solutions for Assignment 1, with one subfolder per problem.

    - assignment2/
    Contains numerical solutions for Assignment 2, with one subfolder per problem.

Each problem folder includes a Python script implementing the corresponding numerical method and simulation.

## Centralized python scripts

- `assignments.py`: Module used to run individual assignment problems from a single interface.

It defines two main functions:

    problems_assignment_1(current_directory, problem_to_run)

    problems_assignment_2(current_directory, problem_to_run)

Each function executes the selected problem script using Python’s runpy module.
Example usage

    from assignments import problems_assignment_1, problems_assignment_2
    import os

    current_directory = os.getcwd()

    ## Run Problem 2 from Assignment 1

    problems_assignment_1(current_directory, problem_to_run=2)

    ## Run Problem 4 from Assignment 2**

    problems_assignment_2(current_directory, problem_to_run=4)

Only valid problem numbers defined in each assignment can be executed. If an invalid number is provided, an error message is printed.

- `main.py`: Entry point for running selected numerical experiments and simulations.

## Utility Module

The utils/ folder contains reusable helper modules shared across multiple assignments:

    SparseMatrices.py
    Functions for assembling and manipulating sparse matrices arising from the discretisation of PDE operators.

    meshGrid.py
    Utilities for generating computational grids and meshes, including structured mesh creation and coordinate handling.

## Topics Covered

Depending on the assignment, the repository includes numerical methods related to:

    Finite difference discretisation

    Sparse matrix assembly

    Time integration schemes

    Stability and convergence analysis

    Boundary and initial value problems for PDEs

## Intended Audience

    Students enrolled in numerical analysis or scientific computing courses

    Researchers seeking reference implementations of basic PDE solvers

    Self-learners studying numerical methods for Partial Differential Equations

## Notes

    The code is intended primarily for educational purposes.

    Scripts reflect the structure and requirements of the original course assignments.

    Mathematical problem descriptions may be expanded in future updates.

## License

This repository is distributed under the GNU General Public License v3.0 (GPL-3.0).

