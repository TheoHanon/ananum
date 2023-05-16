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
  int ** rmbl = remove_bnd_lines(K, &M, boundary_nodes, n_boundary_nodes, symmetry_nodes, n_symmetry_nodes, &bK, &bM, &cbK, coord);
  int * remaining = rmbl[0];
  int * map       = rmbl[1];
  free(rmbl);

  // Power iteration + deflation to find k largest eigenvalues
  double *v = malloc(bM->m * sizeof(double));
  double lambda, freq;
  FILE *file = fopen(argv[2], "w"); // open file to write frequencies

  lu_band(bK);

  lambda = eigen(bM, bK, cbK, v);

  freq = 1. / (2 * M_PI * sqrt(lambda));

  fprintf(file, "%.9lf ", freq);

  printf("lambda = %.9e, f = %.3lf\n", lambda, freq);

  
  for(int i = 1; i < k; i++){
    deflate(cbK,M,v,lambda);

    lambda = eigen_full(M, bK, cbK, v);

    freq = 1. / (2 * M_PI * sqrt(lambda));

    fprintf(file, "%.9lf ", freq);

    printf("lambda = %.9e, f = %.3lf\n", lambda, freq);

  }
  fclose(file);

  file = fopen("displacement.txt", "w");
  size_t n_new = M->n;
  size_t n_rem = K->n - n_new;
  double x, y, dis_x, dis_y;
  for(size_t i = 0; i < n_rem; i++){
    if(i+1 == n_rem || remaining[i+1] != remaining[i]+1){
      //printf("Condition rem[i] = %i\n",remaining[i]);
      for(size_t j = 0; j < n_new; j++){
        if(map[j] == remaining[i]+1){
          x = coord[remaining[i]];
          y = coord[remaining[i]+1];
          dis_x = 0.0;
          dis_y = v[j];
          map[j] = -1;
          break;
        }
      }
    } else {
      x = coord[remaining[i]];
      //printf("rem[i] = %i\n",remaining[i]+1);
      y = coord[remaining[i]+1];
      dis_x = 0.0;
      dis_y = 0.0;
      i++;
    }
    fprintf(file, "%.9lf; %.9lf; %.9lf; %.9lf\n", x, y, dis_x, dis_y);
  }
  for(size_t i = 0; i < n_new; i++){
    if(map[i] == -1){
      continue;
    }
    x = coord[map[i]];
    y = coord[map[i]+1];
    dis_x = v[i];
    i++;
    dis_y = v[i];
    fprintf(file, "%.9lf; %.9lf; %.9lf; %.9lf\n", x, y, dis_x, dis_y);
  }
  fclose(file);

  free(map);
  free(remaining);
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
