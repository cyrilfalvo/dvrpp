extern "C" void dsaupd_(int *ido, char *bmat, int *n, char *which,
			int *nev, double *tol, double *resid, int *ncv,
			double *v, int *ldv, int *iparam, int *ipntr,
			double *workd, double *workl, int *lworkl,
			int *info);

extern "C" void dseupd_(int *rvec, char *All, int *select, double *d,
			double *v1, int *ldv1, double *sigma, 
			char *bmat, int *n, char *which, int *nev,
			double *tol, double *resid, int *ncv, double *v2,
			int *ldv2, int *iparam, int *ipntr, double *workd,
			double *workl, int *lworkl, int *ierr);





