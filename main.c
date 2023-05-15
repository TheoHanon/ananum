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


int main(int argc, char *argv[])
{

  if (argc < 5)
  {
    printf("Usage: \n"
           "./main <r1> <r2> <e> <l> <MeshSizeFactor> \n"
           "---------------------------- \n\n"
           "r1 = inner radius (half-distance between prongs)\n"
           "r2 = outer radius (half-width of fork)\n"
           "e  = length of handle\n"
           "l  = length of prongs\n"
           "meshSizeFactor = meshSize / width of prongs\n"
           "\n");
    return -1;
  }

  // Define physical constants
  double E = 0.7e11; // Young's modulus for Aluminum
  double nu = 0.3;   // Poisson coefficient
  double rho = 3000; // Density of Aluminum

  // Initialize Gmsh and create geometry
  int ierr;
  double r1, r2, e, l, meshSizeFactor;

   r1 = atof(argv[1]);
   r2 = atof(argv[2]);
   e = atof(argv[3]);
   l = atof(argv[4]);
   meshSizeFactor = atof(argv[5]); 

  gmshInitialize(argc, argv, 0, 0, &ierr);
  designTuningFork(r1, r2, e, l, meshSizeFactor,  NULL);
  

  // Assemble the 2 matrices of the linear elasticity problem:
  // M is the mass matrix
  // K is the stiffness matrix
  Matrix *K, *M;
  BandMatrix *bK, *bM, *cbK;
  size_t *boundary_nodes;
  size_t n_boundary_nodes;
  double *coord;
  assemble_system(&K, &M, &boundary_nodes, &n_boundary_nodes,&coord, E, nu, rho);

  
  // Remove lines from matrix that are boundary
  remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &bK, &bM, &cbK, coord);

  


  // Power iteration + deflation to find k largest eigenvalues
  
  double *v = malloc(bM->m * sizeof(double));
  double lambda, freq;
  FILE *file = fopen("out.txt", "w"); // open file to write frequencies

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