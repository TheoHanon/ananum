#ifndef _FEM_MATRICES_H_
#define _FEM_MATRICES_H_

#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmshc.h>
#include <string.h>
#include <stdint.h>

int assemble_system(Matrix** K, Matrix** M, size_t** clamped_nodes, size_t* n_clamped_nodes, size_t** symmetry_nodes, size_t* n_symmetry_nodes, double** coordinates, double E, double nu, double rho);

void visualize_in_gmsh(double* SOL, int n_nodes);

#endif
