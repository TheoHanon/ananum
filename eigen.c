#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "eigen.h"
#include "lu.h"
#include <stdio.h>
#include <string.h>

#define max(a, b) \
  ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a, b) \
  ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

// Remove lines and columns corresponding to boundary nodes from K and M
int** remove_bnd_lines(Matrix *K, Matrix **M, size_t *clamped_nodes, size_t n_clamped_nodes, size_t *symmetry_nodes, size_t n_symmetry_nodes, BandMatrix **K_new, BandMatrix **M_new, BandMatrix **cK, double* coord){
  
  size_t n = K->n;
  size_t n_new;
  int *map;
  int *remaining;
  int perm[n/2]; // = (int *) malloc(sizeof(int)*n/2);
  int band_size;

  for (int i = 0 ; i < (int) n/2; i++) perm[i] = i;
  compute_permutation(perm, coord, (int) n /2);
  
  if(n_symmetry_nodes == 0){
    n_new = K->n - 2*n_clamped_nodes;
  
    map = malloc(n_new * sizeof(int));
    remaining = malloc(2*n_clamped_nodes*sizeof(int));
    int i_bnd = 0;
    int i_new = 0;
    int i_rem = 0;
    int to_write = 0;

    for (int i=0; i< (int) n/2; i++){
      to_write = 1;
      for (i_bnd = 0; i_bnd < n_clamped_nodes; i_bnd++){
        if (perm[i] == clamped_nodes[i_bnd]) {
          to_write = 0;
          remaining[2*i_rem  ] = 2*perm[i];
          remaining[2*i_rem+1] = 2*perm[i]+1;
          i_rem++;
          break;
        } 
      }
      if (to_write){
        map[2*i_new  ] = 2*perm[i];
        map[2*i_new+1] = 2*perm[i]+1;
        i_new++;}
    }
  }

  else {
    size_t i_clm=0,i_sym=0;
    size_t nbr_comm = 0;
    while(i_clm < n_clamped_nodes && i_sym < n_symmetry_nodes){

      if (clamped_nodes[i_clm] < symmetry_nodes[i_sym]) { 
        i_clm++;
      } else if(symmetry_nodes[i_sym] < clamped_nodes[i_clm]){
        i_sym++;
      } else {
        symmetry_nodes[i_sym] = SIZE_MAX;
        i_clm++;
        i_sym++;
        nbr_comm++;
      }
    }

    n_new = K->n - 2*n_clamped_nodes - n_symmetry_nodes + nbr_comm;

    map = malloc(n_new * sizeof(int));
    remaining = malloc((K->n - n_new)*sizeof(int));
    int i_new = 0;
    int i_rem = 0;
    int to_write = 0;


    for (int i=0; i< (int) n/2; i++){
      to_write = 2;
      for (i_clm = 0; i_clm < n_clamped_nodes; i_clm++){
        if (perm[i] == clamped_nodes[i_clm]) {
          to_write = 0;
          remaining[i_rem] = 2*perm[i];
          i_rem++;
          remaining[i_rem] = 2*perm[i]+1;
          i_rem++;
          break;
        }
      }
      if(to_write == 2){
        for (i_sym = 0; i_sym < n_symmetry_nodes; i_sym++){
          if (perm[i] == symmetry_nodes[i_sym]) {
            to_write = 1;
            remaining[i_rem] = 2*perm[i];
            i_rem++;
            break;
          } 
        }
      }
      if (to_write == 1){
        map[i_new] = 2*perm[i]+1;
        i_new++;}
      if (to_write == 2){
        map[i_new] = 2*perm[i];
        i_new++;
        map[i_new] = 2*perm[i]+1;
        i_new++;}
    }
  }

  band_size = 0;
  for (int i = 0 ; i < n_new; i++) {
    for (int j = i ; j < n_new; j++) {
      if (fabs(K -> a[map[i]][map[j]]) > 0.0 || fabs((*M) -> a[map[i]][map[j]]) > 0.0) band_size = max(band_size, abs(i - j));
    }
  }

  printf("Band size = %d\n", band_size);

  *K_new = allocate_band_matrix(n_new, band_size);
  *M_new = allocate_band_matrix(n_new, band_size);
  *cK    = allocate_band_matrix(n_new, band_size);

  for (int i = 0 ; i < n_new; i++){
    for (int j = max(0, i-band_size) ; j <= min(n_new-1, i+band_size); j++) {
        (*K_new) -> a[i][j] = K -> a[map[i]][map[j]];
        (*cK)    -> a[i][j] = K -> a[map[i]][map[j]];
        (*M_new) -> a[i][j] = (*M) -> a[map[i]][map[j]];
    }
  }
  free_matrix(*M);
  *M = allocate_matrix(n_new,n_new);
  for (int i = 0 ; i < n_new; i++){
    for (int j = max(0, i-band_size) ; j <= min(n_new-1, i+band_size); j++) {
        (*M) -> a[i][j] = (*M_new) -> a[i][j];
    }
  }
  int ** res = malloc(2*sizeof(int *));
  res[0] = remaining;
  res[1] = map;
  return res;
}

typedef struct
{
  int i;       // index
  double x, y; // coordinates
} Node;

// Comparateur
int cmp(const void *a, const void *b)
{

  Node *na = (Node *)a;
  Node *nb = (Node *)b;
  double diff = na->y - nb->y;
  if (diff > 0)
    return 1;
  else if (diff < 0)
    return -1;
  diff = na->x - nb->x;
  if (diff > 0)
    return 1;
  else if (diff < 0)
    return -1;
  return 0;
}

int compute_permutation(int *perm, double *coord, int n_nodes)
{

  // // We assume perm is allocated but not initialized
  // for(int i = 0; i < n_nodes; i++)
  // 	perm[i] = i;

  // qsort_r(perm, n_nodes, sizeof(int), coord, cmp);

  // Create Node structs
  Node *nodes = malloc(n_nodes * sizeof(Node));
  for (int i = 0; i < n_nodes; i++)
  {
    nodes[i].i = i;
    nodes[i].x = coord[2 * i];
    nodes[i].y = coord[2 * i + 1];
  }

  // Sort nodes
  qsort(nodes, n_nodes, sizeof(Node), cmp);

  // Fetch permutation (we assume perm is allocated)
  for (int i = 0; i < n_nodes; i++)
    perm[i] = nodes[i].i;

  return 0;
}

// Matrix times vector product
int matrix_times_vector(Matrix *A, double *b, double *res)
{
  int n = A->n;
  for (int i = 0; i < n; i++)
  {
    double b_i = 0.;
    for (int j = 0; j < n; j++)
    {
      b_i += A->a[i][j] * b[j];
    }
    res[i] = b_i;
  }
  return 0;
}

int band_matrix_times_vector(BandMatrix *A, double *v, double *res)
{
  double v_i;
  int b = A->k;
  int n = A->m;

  for (int i = 0; i < n; i++)
  {
    v_i = 0.;
    for (int j = max(0, i - b); j <= min(n - 1, i + b); j++)
    {
      v_i += A->a[i][j] * v[j];
    }
    res[i] = v_i;
  }
  return 0;
}

// inner product of 2 vectors
double inner(double *a, double *b, int n)
{
  double res = 0.;
  for (int i = 0; i < n; i++)
    res += a[i] * b[i];
  return res;
}

// normalize vector
double normalize(double *a, int n)
{
  double norm = 0.;
  for (int i = 0; i < n; i++)
    norm += a[i] * a[i];
  norm = sqrt(norm);
  for (int i = 0; i < n; i++)
    a[i] /= norm;
  return norm;
}

void copy_matrix(Matrix *A, Matrix *cA)
{
  int m, n;

  m = A->m;
  n = A->n;
  memcpy(cA->data,A->data,m*n*sizeof(double));
  /*for (int i = 0; i < m; i++)
  {
    for (int j = max(0, i - k); j <= min(m - 1, i + k); j++)
    {
      cA->a[i][j] = A->a[i][j];
    }
  }*/
  return;
}


void copy_band_matrix(BandMatrix *A, BandMatrix *cA)
{
  int m, k;

  m = A->m;
  k = A->k;
  memcpy(cA->data,A->data,m*(2*k+1)*sizeof(double));
  /*for (int i = 0; i < m; i++)
  {
    for (int j = max(0, i - k); j <= min(m - 1, i + k); j++)
    {
      cA->a[i][j] = A->a[i][j];
    }
  }*/
  return;
}

void shift(BandMatrix *M, BandMatrix *K, double lambda)
{
  int n, k;

  n = M->m;
  k = M->k;

  for (int j = 0; j < n; j++)
  {
    M->a[j][j] -= lambda * K->a[j][j];
    for (int l = j + 1; l <= min(j + k, n - 1); l++)
    {
      M->a[j][l] -= lambda * K->a[j][l];
      M->a[l][j] -= lambda * K->a[j][l];
    }
  }
  return;
}

void shift_full(Matrix *M, BandMatrix *K, double lambda)
{
  int n, k;

  n = K->m;
  k = K->k;

  for (int j = 0; j < n; j++)
  {
    M->a[j][j] -= lambda * K->a[j][j];
    for (int l = j + 1; l <= min(j + k, n - 1); l++)
    {
      M->a[j][l] -= lambda * K->a[j][l];
      M->a[l][j] -= lambda * K->a[l][j];
    }
  }
  return;
}



double get_eigenvalue(BandMatrix *M, BandMatrix *K, double * eigen_vector){
    
    int n = M->m;
    double res[n], vMv, vKv;
 

    band_matrix_times_vector(M, eigen_vector, res);
    vMv = inner(eigen_vector, res, n);

    band_matrix_times_vector(K, eigen_vector, res);
    vKv = inner(eigen_vector, res, n);

    return vMv/vKv;
}


double get_eigenvalue_full(Matrix *M, BandMatrix *K, double * eigen_vector){
    
    int n = M->m;
    double res[n], vMv, vKv;
 

    matrix_times_vector(M, eigen_vector, res);
    vMv = inner(eigen_vector, res, n);

    band_matrix_times_vector(K, eigen_vector, res);
    vKv = inner(eigen_vector, res, n);

    return vMv/vKv;
}


double eigen(BandMatrix *bM, BandMatrix *bK, BandMatrix *ccK, double *v)
{
  int n = bM->m, k;
  double Mv[n], lambda, lambda_prev, diff, rtol;
  BandMatrix *cM;
  // Initialization

  lambda = 1;
  rtol = 1e-5;
  k = bM->k;
  cM = allocate_band_matrix(n, k);

  //lu_band(bK);
  for (int i = 0; i < n; i++)
    v[i] = 0;
  v[0] = v[1] = 1;
  normalize(v, n);

  // Power Iteration
  for (int it = 0; it < 1e4; it++)
  {
    lambda_prev = lambda;
    band_matrix_times_vector(bM, v, Mv);
    solve_band(bK, Mv);
    normalize(Mv, n);
    band_matrix_times_vector(bM, Mv, v);
    //printf("\n===== Iteration %d =====\n", it+1);
    lambda = inner(Mv, v, n);
    //printf("λ = %.9e\n", lambda);
    diff = fabs((lambda - lambda_prev) / lambda_prev);
    for (int i = 0; i < n; i++)
      v[i] = Mv[i];
    if (diff < rtol)
      break;
  }


  //printf("λ = %.9e\n", lambda);

  lambda = get_eigenvalue(bM, ccK, v);

  rtol = 1e-9;
  

  // Rayleigh iteration

  for (int it = 0; it < 1e4; it++)
  {
    lambda_prev = lambda;
    copy_band_matrix(bM, cM);
    shift(cM, ccK, lambda);
    lu_band(cM);
    band_matrix_times_vector(ccK, v, Mv);
    solve_band(cM, Mv);
    normalize(Mv, n);
    //printf("\n===== Iteration %d =====\n", it+1);
    lambda = get_eigenvalue(bM, ccK, Mv);
    //printf("λ = %.9e\n", lambda);
    diff = fabs((lambda - lambda_prev) / lambda_prev);
    for (int i = 0; i < n; i++) v[i] = Mv[i];
    if (diff < rtol)
      break;
  }

  free_band_matrix(cM);
  return lambda;
}



double eigen_full(Matrix *M, BandMatrix *bK, BandMatrix *ccK, double *v)
{
  int n = M->m;
  double Mv[n], lambda, lambda_prev, diff, rtol;
  Matrix *cM;
  // Initialization

  lambda = 1;
  rtol = 1e-5;
  cM = allocate_matrix(n, n);

  //lu_band(bK);
  for (int i = 0; i < n; i++)
    v[i] = 0;
  v[0] = v[1] = 1;
  normalize(v, n);

  // Power Iteration
  for (int it = 0; it < 1e4; it++)
  {
    lambda_prev = lambda;
    matrix_times_vector(M, v, Mv);
    solve_band(bK, Mv);
    normalize(Mv, n);
    matrix_times_vector(M, Mv, v);
    //printf("\n===== Iteration %d =====\n", it+1);
    lambda = inner(Mv, v, n);
    //printf("λ = %.9e\n", lambda);
    diff = fabs((lambda - lambda_prev) / lambda_prev);
    for (int i = 0; i < n; i++)
      v[i] = Mv[i];
    if (diff < rtol)
      break;
  }


  lambda = get_eigenvalue_full(M, ccK, v);
  //printf("λ = %.9e\n", lambda);
  rtol = 1e-9;
  

  // Rayleigh iteration

  for (int it = 0; it < 1e4; it++)
  {
    lambda_prev = lambda;
    copy_matrix(M, cM);
    shift_full(cM, ccK, lambda);
    lu(cM);
    band_matrix_times_vector(ccK, v, Mv);
    solve(cM, Mv);
    normalize(Mv, n);
    //printf("\n===== Iteration %d =====\n", it+1);
    lambda = get_eigenvalue_full(M, ccK, Mv);
    //printf("λ = %.9e\n", lambda);
    diff = fabs((lambda - lambda_prev) / lambda_prev);
    for (int i = 0; i < n; i++) v[i] = Mv[i];
    if (diff < rtol)
      break;
  }

  free_matrix(cM);
  return lambda;
}



double band_power_iteration(BandMatrix *bM, BandMatrix *bK, BandMatrix *ccK, double *v)
{
  int n = bM->m;
  lu_band(bK);

  double Mv[n];

  for (int i = 0; i < n; i++)
    v[i] = 0;
  v[0] = v[1] = 1;
  normalize(v, n);

  double lambda = 1, lambda_prev, diff;
  double rtol = 1e-9; // relative tolerance on lambda

  for (int it = 0; it < 1e4; it++)
  {
    lambda_prev = lambda;
    band_matrix_times_vector(bM, v, Mv);
    solve_band(bK, Mv);
    normalize(Mv, n);
    band_matrix_times_vector(bM, Mv, v);
    lambda = inner(Mv, v, n);
    diff = fabs((lambda - lambda_prev) / lambda_prev);
    for (int i = 0; i < n; i++)
      v[i] = Mv[i];
    if (diff < rtol)
      break;
  }
  band_matrix_times_vector(ccK, v, Mv);
  lambda /= inner(v, Mv, n);

  return lambda;
}

double power_iteration(Matrix *A, double *v)
{
  int n = A->n;

  // initial guess
  double Av[n];
  for (int i = 0; i < n; i++)
    v[i] = 0;
  v[0] = v[1] = 1;
  normalize(v, n);

  double lambda = 1, lambda_prev, diff;
  double rtol = 1e-9; // relative tolerance on lambda

  for (int it = 0; it < 1e4; it++)
  {
    lambda_prev = lambda;
    // printf("\n===== Iteration %d =====\n", it+1);
    matrix_times_vector(A, v, Av);
    lambda = inner(v, Av, n);
    // printf("λ = %.9e\n", lambda);
    diff = fabs((lambda - lambda_prev) / lambda_prev);
    for (int i = 0; i < n; i++)
      v[i] = Av[i];
    normalize(v, n);
    if (diff < rtol)
      break;
  }

  return lambda;
}


void deflate(BandMatrix *K, Matrix *M, double *eigenvect, double lambda){
  int m = K->m;
  double Kv[m];
  band_matrix_times_vector(K,eigenvect,Kv);
  //double fac = inner(eigenvect,eigenvect,m);
  for(int i = 0; i < m; i++){
    for (int j = 0; j < m; j++){
      M -> a [i][j] -= lambda * Kv[i] * eigenvect[j]; // / fac;
    }
  }
}