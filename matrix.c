#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "lu.h"

#include "matrix.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

	 

Matrix * allocate_matrix(int m, int n) {
	Matrix * mat = (Matrix*) malloc(sizeof(Matrix));
	mat->m = m, mat->n = n;
	mat->data = (double*) malloc(m*n*sizeof(double));
	mat->a = (double**) malloc(m*sizeof(double*));
	for(int i = 0; i < m; i++)
		mat->a[i] = mat->data+i*n;
	return mat;
}


BandMatrix * allocate_band_matrix(int m, int k) {
	// To make things simple, we allocate 2k+1 for each row
	// (even though we need less than that for the k first and last rows).
	// We want to have a[i][i] == data[k + i*(2*k+1)]

	BandMatrix * mat = (BandMatrix*) malloc(sizeof(BandMatrix));
	mat->m = m, mat->k = k;
	mat->data = (double*) calloc(m*(2*k+1), sizeof(double));
	mat->a = (double**) malloc(m*sizeof(double*));
	// This is the tricky part :-)
	// We want a[i][i] == data[(2*k+1)*i + k]
	// which means a[i] + i == data + (2*k+1)*i + k
	// and therefore a[i] == data + k + 2*k*i
	for(int i = 0; i < m; i++)
		mat->a[i] = mat -> data + k + 2*k*i;
	
	return mat;
}

int inverse_band(BandMatrix* K , BandMatrix* M) {
	int m = M -> m;
	int b = M -> k;
	
	int err = lu_band(K);
	double* temp = malloc(m * sizeof(double));
	for (int k = 0; k < m; k++) {
		for (int i = max(0, k-b); i <= min(m-1, i+b); i++) temp[i] = M -> a[i][k];
		solve_band(K,temp);
		for (int i = max(0, k-b); i <= min(m-1, i+b); i++) M -> a[i][k] = temp[i];
	}
	free(temp);
	return err;

}

int solve_band(BandMatrix * LU, double * y) {
	int m = LU->m, b = LU->k;
	double * x = y;

	// Résolution de L*x = y par substitution avant
	for(int i = 0; i < m; i++) {
		for(int j = max(0, i-b); j < i; j++) {
			x[i] -= LU->a[i][j] * x[j];
		}
	}

	// Résolution de U*x = L^{-1} y par substitution arrière
	for(int i = m-1; i >= 0; i--) {
		for(int j = i+1; j <= min(m-1, i+b); j++)
			x[i] -= LU->a[i][j] * x[j];
		x[i] /= LU->a[i][i];
	}

	return 0;
}

int lu_band(BandMatrix * A) {
	int m = A->m, b = A->k;
	BandMatrix * LU = A; // factorization is done in-place
	for(int k = 0; k < m-1; k++) {
		if(fabs(LU->a[k][k]) < EPS) return -1;
		for(int j = k+1; j <= min(m-1, k+b); j++) {
			LU->a[j][k] /= LU->a[k][k];
			for(int l = k+1; l <= min(m-1, k+b); l++)
				LU->a[j][l] -= LU->a[j][k] * LU->a[k][l];
		}
	}
	return 0;
}


void free_matrix(Matrix * mat) {
	free(mat->a);
	free(mat->data);
	free(mat);
}

void print_vector(double * v, int n) {
	for(int i = 0; i < n; i++)
		printf("%.3e ", v[i]);
	printf("\n");
}

void print_matrix(Matrix * A) {
	for(int i = 0; i < A->m; i++)
		print_vector(A->a[i], A->n);
}

int is_symmetric(Matrix * K) {
	int symmetric = 1;
	for(int i = 0; i < K->m; i++)
		for(int j = i+1; j < K->n; j++)
			if(fabs((K->a[i][j] - K->a[j][i]) / K->a[i][j]) > 1e-12) {
				printf("%d %d\n", i, j);
				printf("%lf %lf\n", K->a[i][j], K->a[j][i]);
				symmetric = 0;
			}
	return symmetric;
}


int inverse(Matrix* K, Matrix* M){
	int m = K -> m;
	int err = lu(K);
	double* temp = malloc(m * sizeof(double));
	for (int k = 0; k < m; k++) {
		for (int i = 0; i < m; i++) temp[i] = M -> a[i][k];
		solve(K,temp);
		for (int i = 0; i < m; i++) M -> a[i][k] = temp[i];
	}
	free(temp);
	return err;
}

void free_band_matrix(BandMatrix * mat) {
	free(mat->a);
	free(mat->data);
	free(mat);
}
