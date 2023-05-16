#include <stdint.h>

void remove_bnd_lines(Matrix *K, Matrix *M, size_t *clamped_nodes, size_t n_clamped_nodes, size_t *symmetry_nodes, size_t n_symmetry_nodes, BandMatrix **K_new, BandMatrix **M_new, BandMatrix **cK, double* coord);

// inner product of 2 vectors
double inner(double* a, double* b, int n);

// normalize vector 
double normalize(double* a, int n);

double power_iteration(Matrix *A, double *v);

int compute_permutation(int * perm, double * coord, int n_nodes) ;
void boundary(size_t* boundary_nodes, size_t n_boundary_nodes, int* perm, double* coord, int n_nodes, int* perm_size);

double band_power_iteration(BandMatrix *bM, BandMatrix *bK, BandMatrix* ccK, double *v);
int cmp(const void * a, const void * b);

double eigen(BandMatrix *bM, BandMatrix *bK, BandMatrix *ccK, double *v);
void shift(BandMatrix *M, BandMatrix *K, double lambda);
void copy_band_matrix(BandMatrix *A, BandMatrix *cA);
double get_eigenvalue(BandMatrix *M, BandMatrix *K, double * eigen_vector);

void deflate(BandMatrix *K, BandMatrix *M, double *eigenvect, double lambda);