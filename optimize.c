#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "elasticity.h"
#include <gmshc.h>
#include "optimize.h"
#include "math.h"
#include "lu.h"
#include "design.h"
#include "eigen.h"



double f_target(double r1, double r2, double e, double l, double meshSizeFactor) {

    double E = 0.7e11;  // Young's modulus for Aluminum
    double nu = 0.3;    // Poisson coefficient
    double rho = 3000;
	// Initialize Gmsh and create geometry
    int ierr;
    gmshInitialize(0, NULL, 0,0, &ierr);
    designTuningFork(r1, r2, e, l, meshSizeFactor, NULL);


    Matrix *K, *M;
    size_t* boundary_nodes;
    size_t n_boundary_nodes;
    double * coord;
    assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);
    
    Matrix *K_new;
    Matrix *M_new;
    remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);
    inverse(K_new, M_new);

    double * v = malloc(M_new->m * sizeof(double));
    double lambda, freq;
    lambda = power_iteration(M_new, v);
    freq = 1./(2*M_PI*sqrt(lambda));
    printf("freq = %f\n", freq);
    gmshFinalize(&ierr);
    free(v);
    free_matrix (K);
    free_matrix (M);
    free_matrix (K_new);
    free_matrix (M_new);
    free(boundary_nodes);
    free(coord);
    

    return freq;

}