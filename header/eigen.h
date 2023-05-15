void remove_bnd_lines (Matrix *K, Matrix *M, size_t *bnd_nodes, size_t n_bnd_nodes, BandMatrix **K_new, BandMatrix **M_new, BandMatrix **cK, double* coord);

// inner product of 2 vectors
double inner(double* a, double* b, int n);

// normalize vector 
double normalize(double* a, int n);

double power_iteration(Matrix *A, double *v);

int compute_permutation(int * perm, double * coord, int n_nodes) ;
void boundary(size_t* boundary_nodes, size_t n_boundary_nodes, int* perm, double* coord, int n_nodes, int* perm_size);

double band_power_iteration(BandMatrix *bM, BandMatrix *bK, BandMatrix* ccK, double *v);
int cmp(const void * a, const void * b);