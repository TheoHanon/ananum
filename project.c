#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include "matrix.h"
#include "elasticity.h"
#include "math.h"
#include "lu.h"
#include "design.h"
#include "eigen.h"
#include "optimize.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

int main (int argc, char *argv[]) {

  if (argc < 2){
    printf("Usage: \n"
			"./project <k> <out>\n" 
			"---------------------------- \n\n"
			"- k is the number of frequencies to compute. \n "
			"- out is the output file to write the frequencies. \n "
      "\n");
		return -1;
  }

    // Define physical constants
  double E = 0.7e11; // Young's modulus for Aluminum
  double nu = 0.3;   // Poisson coefficient
  double rho = 3000; // Density of Aluminum

  // Initialize Gmsh and create geometry
  int ierr;
  double r1 = 6e-3;
  double r2 = 11e-3;
  double e  = 19e-3;
  double l  = 82e-3;
  double c  = 19e-3;

  
  gmshInitialize(argc, argv, 0, 0, &ierr);
  designTuningForkSym2(r1, r2, e, l, c, 0.3, NULL);
  
  // Number of vibration modes to find
  int k = atoi(argv[1]);


  // Assemble the 2 matrices of the linear elasticity problem:
  // M is the mass matrix
  // K is the stiffness matrix
  Matrix *K, *M;
  BandMatrix *bK, *bM, *cbK;
  size_t* boundary_nodes;
  size_t n_boundary_nodes;
  size_t* symmetry_nodes;
  size_t n_symmetry_nodes;
  double *coord;
  assemble_system(&K, &M, &boundary_nodes, &n_boundary_nodes, &symmetry_nodes, &n_symmetry_nodes, &coord, E, nu, rho);

  // Remove lines from matrix that are boundary
  remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, symmetry_nodes, n_symmetry_nodes, &bK, &bM, &cbK, coord);


  // Power iteration + deflation to find k largest eigenvalues
  double *v = malloc(bM->m * sizeof(double));
  double lambda, freq;
  FILE *file = fopen(argv[2], "w"); // open file to write frequencies

  lu_band(bK);
  /*
  for(int i = 0; i < k; i++){
    lambda = eigen(bM, bK, cbK, v);

    freq = 1. / (2 * M_PI * sqrt(lambda));

    fprintf(file, "%.9lf ", freq);

    printf("lambda = %.9e, f = %.3lf\n", lambda, freq);

    if(i != k-1) deflate(cbK,bM,v,lambda);
  }*/
  lambda = eigen(bM, bK, cbK, v);

  freq = 1. / (2 * M_PI * sqrt(lambda));

  fprintf(file, "%.9lf ", freq);

  printf("lambda = %.9e, f = %.3lf\n", lambda, freq);


  fclose(file);
  free_matrix(K);
  free_matrix(M);
  free_band_matrix(bK);
  free_band_matrix(bM);
  free_band_matrix(cbK);
  free(boundary_nodes);
  free(coord);
  free(v);
  return 0;
}
