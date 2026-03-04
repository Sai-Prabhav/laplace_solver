# Laplace Solver

A numerical solver for the Laplace equation using the deal.II finite element library.

## Description

This project solves the Laplace equation using the Finite Element Method with the deal.II library:

$$\nabla^2 u = 0$$

## Building and Running

```bash
cmake .
make run
```

## Example Problems

### Problem 1 (2D)

Solves Laplace's equation on a rectangular domain $[0, \pi] \times [0, 2\pi]$ with boundary conditions:

- $u(0, y) = 0$
- $u(w, y) = 0$
- $u(x, 0) = 0$
- $u(x, h) = \sin(x)$

The exact solution is:

$$u(x,y) = \frac{\sinh(y)}{\sinh(h)} \sin(x)$$

### Visualization

![Heatmap of solution](https://github.com/Sai-Prabhav/laplace_solver/blob/main/problem-1_4mil.png)

### Errors

| Cycle |    L2 Error |    H1 Error |
| :---: | ----------: | ----------: |
|   1   |    0.280613 |    0.829229 |
|   2   |    0.099025 |    0.433508 |
|   3   |   0.0276992 |    0.216964 |
|   4   |  0.00713487 |    0.108348 |
|   5   |  0.00179731 |   0.0541507 |
|   6   | 0.000450184 |   0.0270723 |
|   7   | 0.000112602 |   0.0135357 |
|   8   | 2.81711e-05 |  0.00676782 |
|   9   | 7.08543e-06 |  0.00338391 |
|  10   | 1.77896e-06 |  0.00169196 |
|  11   | 4.49083e-07 | 0.000845995 |

### Problem 2 (3D)

Solves Laplace's equation on a unit cube $[0,1]^3$ with solution $u(x,y,z) = \sin(\pi x)\,\sin(\pi y)\,\sinh(\sqrt{2}\,\pi z)$.

## Errors

| Cycle |    L2 Error | H1 Error |
| :---: | ----------: | -------: |
|   1   |     2.03909 |  23.9678 |
|   2   |    0.537871 |  12.3252 |
|   3   |    0.136677 |  6.20055 |
|   4   |   0.0343238 |   3.1048 |
|   5   |  0.00859095 |  1.55296 |
|   6   |  0.00214837 |  0.77655 |
|   7   | 0.000537227 | 0.388284 |

### Visualization

![Contour Plot o the solution](https://github.com/Sai-Prabhav/laplace_solver/blob/main/problem-2_2m.png)

![Animation of the cross sectional view of the solution](https://github.com/Sai-Prabhav/laplace_solver/blob/main/problem-2.gif)

In both the above problems we can clearly observe that L2 error is of order $O(h^2)$ where H1 error is of the order $O(h)$.
