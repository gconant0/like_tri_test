//c++ include file for the library lbprsh_lapack
//Computes the inverse of a matrix--requires LU factored matrix
//from dgetrf
extern "C" int dgetri_(long int *n, double *a, long int *lda, 
		       long int *ipiv, double *work, long int *lwork, 
		       long int *info);

//Computes the LU factorization of a matrix
extern "C" int dgetrf_ (long int *m, long int *n, double *a, long int *lda, 
			long int *ipiv, long int *info);

//Computes the eigen values and vectors of a matrix 
extern "C" int dgeev_(char *jobvl, char *jobvr, long int *n, 
		  double *a, long int *lda, double *wr, 
		  double *wi, double *vl, long int *ldvl, 
		  double *vr, long int *ldvr, double *work, 
		  long int *lwork, long int *info, long int *jobvl_len, 
		  long int *jobvr_len);

//Performs a matrix/matrix multiplication
extern "C" int dgemm_ (char *transa, char *transb, long int *m, long int *n, 
		       long int *k, double *alpha, double *a, 
		       long int *lda, double *b, long int *ldb, 
		       double *beta, double *c__, long int *ldc, 
		       long int *transa_len, long int *transb_len);




