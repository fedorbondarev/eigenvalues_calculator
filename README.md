# Eigenvalues Calculator using Double Step QR Algorithm

This C program calculates the eigenvalues of a matrix using the Double Step QR Algorithm.

## Overview

The Double Step QR Algorithm is an iterative method used to compute the eigenvalues of a matrix. It iteratively transforms the matrix into a quasi-triangular form, converging towards the eigenvalues. This implementation utilizes the Double Step QR Algorithm to efficiently calculate the eigenvalues of a given matrix.

## Usage

Input and output files are passed as first and second arguments.

## Matrix Format

The input file contains the matrix in a specific format:

```
n
a11 a12 a13 ... a1n
a21 a22 a23 ... a2n
.   .   .       .
.   .   .       .
.   .   .       .
an1 an2 an3 ... ann
```

Where:
- `n` is the size of the square matrix (number of rows/columns).
- `aij` represents the elements of the matrix.

For example:

```
3
2.3 1.2 0.5
1.2 3.4 1.1
0.5 1.1 2.9
```