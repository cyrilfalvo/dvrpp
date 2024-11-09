
extern "C" void dlaruv_(int iseed[4], const int &n, double* x );

extern "C" void dlarnv_(const int &idist, int iseed[4], const int &n, double* x );

extern "C" double dlaran_(int iseed[4]);

extern "C" void dgesvd_(
  const char &jobu,
  const char &jobvt,
  const int &m,
  const int &n,
  double *a,
  const int &lda,
  double *s,
  double *u,
  const int &ldu,
  double *vt,
  const int &ldvt,
  double* work,
  const int &lwork,
  int &info
  );

extern "C" void dsysv_(
  const char &job,
  const int &n,
  const int &nrhs,
  double *a,
  const int &lda,
  int *ipiv,
  double *b,
  const int &ldb,
  double* work,
  const int &lwork,
  int &info
  );

extern "C" void dsyev_(
  const char &job,
  const char &uplo,
  const int &n,
  double *a,
  const int &lda,
  double *w,
  double *work,
  const int &lwork,
  int &info
  );

extern "C" void dgemm_(
	const char &transa,
	const char &transb,
	const int &m,
	const int &n,
	const int &k,
	const double &alpha,
	double* A,
	const int &lda,
	double* B,
	const int &ldb,
	const double &beta,
	double* C,
	const int &ldc
	);

extern "C" double ddot_(
	const int &n,
	double* dx,
	const int &incx,
	double* dy,
	const int &incy
	);


