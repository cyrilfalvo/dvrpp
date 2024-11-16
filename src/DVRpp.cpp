
#include "../include/DVRpp.h"
#include "../include/BasisHerm.h"
#include "../include/BasisSine.h"
#include "../include/BasisExp.h"
#include <cmath>
#include <cstdlib>
#include "../include/lapack.h"
#ifdef _ARPACK_
  #include "../include/arpack.h"
#endif
#include <cstdio>
#include <ctime>

DVR::DVR() {
  Basis = 0;
  Dim = 0;
  Size = 0;
  Vi = 0;
  EVectors = 0;
  EEnergies = 0;
  NGrid = 0;
  Pop = 0;
  Fact = 0;
  flag = 1;
}

DVR::~DVR() {
  if(Basis!=0) {
    for(int i=0;i<Dim;i++) { if(Basis[i] != 0 ) { delete Basis[i]; } }
    delete [] Basis;
    Basis = 0;
  }
  if (Basis!=0) { delete [] Basis; }  
  if (Vi != 0) { delete [] Vi; }
  if ( EVectors != 0) { delete [] EVectors; }
  if ( EEnergies != 0) { delete [] EEnergies; }	
  if ( NGrid != 0) { delete [] NGrid; }
  if ( Pop != 0 ) { delete [] Pop; }
  if ( Fact != 0 ) { delete [] Fact; }
  Dim = 0;
  Size = 0;
  Vi = 0;
  Basis = 0;
  NGrid = 0;
  Pop = 0;
  Fact = 0;
  flag = 1;
}

int DVR::GetSize() {
  return Size;
}

void DVR::SetDimension(int DimV) {
  Dim = DimV;
  Basis = new BasisSet*[Dim];
  NGrid = new int[Dim];
  for(int i=0;i<Dim;i++) {
    Basis[i] = 0;
    NGrid[i] = 0;
  }
}

int DVR::AddHermBasis(int mode, int N,double Fact, double Delta, double Xeq) {
  if ( (mode>=Dim) && (mode<0) ) { return EXIT_FAILURE; }
  if (Basis[mode] != 0 ) { return EXIT_FAILURE; }
  NGrid[mode] = N;
  Basis[mode] = new BasisHerm();
  Basis[mode]->SetFact(Fact);
  Basis[mode]->SetBasisSet(N,Delta,Xeq);
  return EXIT_SUCCESS;
}

int DVR::AddSineBasis(int mode, int N,double Fact, double Delta, double Xeq) {
  if ( (mode>=Dim) && (mode<0) ) { return EXIT_FAILURE; }
  if (Basis[mode] != 0 ) { return EXIT_FAILURE; }
  NGrid[mode] = N;
  Basis[mode] = new BasisSine();
  Basis[mode]->SetFact(Fact);
  Basis[mode]->SetBasisSet(N,Delta,Xeq);
  return EXIT_SUCCESS;
}

int DVR::AddExpBasis(int mode, int N,double Fact, double Delta, double Xeq) {
  if ( (mode>=Dim) && (mode<0) ) { return EXIT_FAILURE; }
  if (Basis[mode] != 0 ) { return EXIT_FAILURE; }
  NGrid[mode] = N;
  Basis[mode] = new BasisExp();
  Basis[mode]->SetFact(Fact);
  Basis[mode]->SetBasisSet(N,Delta,Xeq);
  return EXIT_SUCCESS;
}


int DVR::PrepareDvr() {
  for(int i=0;i<Dim;i++) { 
    if (Basis[i]==0) { return EXIT_FAILURE; }
  }
  Size = 1;
  for(int i=0;i<Dim;i++) {
    Size *= NGrid[i];
  }
  
  Vi = new double[Size];
  for(int i=0;i<Size;i++) { Vi[i] = 0.; }
  
  Pop = new int[Size*Dim];
  for(int s1=0;s1<Size;s1++) {
    int index1 = s1;
    for(int d=0;d<Dim;d++) {
      Pop[s1*Dim+d] = index1%NGrid[d];
      index1/= NGrid[d];
    } 
  }

  Fact = new int[Dim];
  Fact[0] = 1;
  for(int d=1;d<Dim;d++) { Fact[d] = Fact[d-1]*NGrid[d-1]; }

  flag = 0; // everything is ok
  
  return EXIT_SUCCESS;
}

void DVR::LoadPotential(int index, double V) {
  if(flag==1) { return ; }
  Vi[index] = V;
}
 
void DVR::LoadPotential(double* V) {
  if(flag==1) { return ; }
  for(int i=0;i<Size;i++) {
    Vi[i] = V[i];
  }
}

void DVR::GetNodeCoord(int index, double* R) {
  if(flag==1) { return ; }
  for(int i=0;i<Dim;i++) {
    R[i] = Basis[i]->Xi[index % NGrid[i]];
    index /= NGrid[i];
  }
}

void DVR::MultHamiltonian(double* d1, double* d2) {
  double* Tq;
  for(int i=0;i<Size;i++) {
    d2[i] = Vi[i]*d1[i];
  }
  for(int b=0;b<Dim;b++) {
    Tq = Basis[b]->Tq;
    int fact = Fact[b];
    for(int s1=0;s1<Size;s1++) {
      int vv1 = Pop[s1*Dim+b];
      for(int vv2 = 0; vv2<NGrid[b];vv2++) {
	int s2 = s1 + (vv2-vv1)*fact;
	d2[s2] += Tq[vv1+NGrid[b]*vv2]*d1[s1];
      }
    }
   }
}

int DVR::SolveHamiltonianLapack() { 
  if(flag==1) { return EXIT_FAILURE; }

  if(EVectors != 0 ) { delete [] EVectors; }
  if(EEnergies != 0 ) { delete [] EEnergies; } 

  //printf("Eigenvectors storage will takes %.6f Mbytes\n", (double)Size*Size*sizeof(double)/1024/1024);
  
  EEnergies = new double[Size];
  EVectors = new double[Size*Size];
  
  for(int i=0;i<Size;i++) { EEnergies[i] = 0.; }
  for(int i=0;i<Size*Size;i++) { EVectors[i] = 0.; }
  
  double* Tq;
  int* pop = new int[Dim];
  
  for(int b=0;b<Dim;b++) {
    Tq = Basis[b]->Tq;
    int fact = Fact[b];
    for(int s1=0;s1<Size;s1++) {
      int vv1=Pop[s1*Dim+b];
      for(int vv2 = 0; vv2<NGrid[b];vv2++) {
	int s2 = s1 + (vv2-vv1)*fact;
	EVectors[s1+Size*s2] += Tq[vv1+NGrid[b]*vv2];
      }
    }
  }
  delete [] pop;
  
  for(int i=0;i<Size;i++) {
    EVectors[i+Size*i] += Vi[i];
  }
  
  int lwork;
  int info;
  double* work;
  
  lwork=-1;
  work = new double[2];  
  dsyev_('V','U',Size,EVectors,Size,EEnergies,work,lwork,info);
  lwork = (int) work[0];
  delete [] work;
  work = new double[lwork];
  dsyev_('V','U',Size,EVectors,Size,EEnergies,work,lwork,info);
  delete [] work;
  
  if(info==0) { return EXIT_SUCCESS; } 
  else { printf("Error in the Lapack diagonalization\n"); return EXIT_FAILURE; }
}


#ifdef _ARPACK_
int DVR::SolveHamiltonianArpack(int nev, int ncv) {
  if(flag==1) { return EXIT_FAILURE; }
  
  if(EVectors != 0 ) { delete [] EVectors; }
  if(EEnergies != 0 ) { delete [] EEnergies; } 
  
  clock_t clock_start,clock_end;
  
  clock_start = clock();
  
  //printf("Eigenvectors storage will takes %.6f Mbytes\n", (double)Size*nev*sizeof(double)/1024/1024);
  
  EEnergies = new double[nev];
  EVectors = new double[Size*nev];
  
  for(int i=0;i<nev;i++) { EEnergies[i] = 0.; }
  for(int i=0;i<Size*nev;i++) { EVectors[i] = 0.; }
  
  if (ncv>Size) { ncv = Size; }  
  
  int ido = 0;
  char bmat[2] = "I";
  char which[3] = "SA";
  char HowMny[4] = "All";
  double tol = 0.0;
  double* resid = new double[Size];
  double* v = new double[Size*ncv];
  int iparam[11];
  int ipntr[11];
  iparam[0] = 1;
  iparam[2] = 300;
  iparam[6] = 1;
  double* workd = new double[3*Size];
  int lworkl = ncv*(ncv+8);
  double* workl = new double[lworkl];
  int info = 0;
  int rvec = 1;
  int* selec = new int[ncv];
  for(int i=0;i<ncv;i++) { selec[i] = 0; }
  double* d = new double[nev];
  double sigma;
  int ierr;
  
  dsaupd_(&ido,bmat,&Size,which,&nev,&tol,resid,&ncv,v,&Size,iparam,ipntr,workd,workl,&lworkl,&info);
  while((ido==-1)||(ido==1)) {
    MultHamiltonian(workd+ipntr[0]-1,workd+ipntr[1]-1);  
    dsaupd_(&ido,bmat,&Size,which,&nev,&tol,resid,&ncv,v,&Size,iparam,ipntr,workd,workl,&lworkl,&info);
  }
  if(info != 0 ) { printf("error code %d with procedure dsaupd_\n",info); return EXIT_FAILURE; }
  dseupd_(&rvec,HowMny,selec,EEnergies,EVectors,&Size,&sigma,bmat,&Size,which,&nev,&tol,resid,&ncv,v,&Size,iparam,ipntr,workd,workl,&lworkl,&ierr);
  if(ierr != 0 ) { printf("error code %d with procedure dseupd_\n",ierr); return EXIT_FAILURE; }
  
  delete [] resid;
  delete [] v;
  delete [] workd;
  delete [] workl;
  delete [] selec;
  delete [] d;
  
  clock_end = clock();
  
  double delay = (double) (clock_end-clock_start)/(double) CLOCKS_PER_SEC;
  
  printf("Arpack diagonalization performed in %12.3f seconds\n",delay);
  

  return EXIT_SUCCESS;
}
#endif

double DVR::GetEnergy(int i) {
  if (EEnergies == 0) { return 0.;}
  return EEnergies[i];
}

int DVR::GetVector(int i, double* Psi) {
  if (EVectors == 0) { return EXIT_FAILURE; } 
  double* Pt = EVectors + i*Size;
  for(int j=0;j<Size;j++) {
    Psi[j] = Pt[j];
  }
  return EXIT_SUCCESS;
}

int DVR::GetPop(int i, int* ni) {
  if(Pop==0) { return EXIT_FAILURE; }
  if(i>=Size) { return EXIT_FAILURE; }
  for(int j=0;j<Dim;j++) {
    ni[j] = Pop[i*Dim+j];
  }
  return EXIT_SUCCESS;
}

/*
int DVR::GetElmQOpt(int PowMax, int i, double* AvgOpt) {
  if(AvgOpt==0) { return EXIT_FAILURE; }
  if(EVectors == 0 ) { return EXIT_FAILURE; }
  double* Pt1 = EVectors + i*Size;
  double fact;
  for(int d=0;d<PowMax;d++) {
    AvgOpt[d] = 0.;
  }
  for(int k=0;k<Size;k++) {
    fact = Pt1[0]*Pt1[0];
    for(int d=0;d<PowMax;d++) {
      fact*=Basis[0]->Xi[k];
      AvgOpt[d]+=fact;
    }   
    Pt1++;
  }
  return EXIT_SUCCESS;
}
*/

double DVR::GetElmDiagOperatorSingle(int i, int j, double* Ai) {
  if (Ai == 0) { return 0.; }
  if (EVectors == 0) { return 0.; }
  double* Pt1 = EVectors + i*Size;
  double* Pt2 = EVectors + j*Size; 
  double* Pt3 = Ai;
  double Elm = 0;
  for(int k=0;k<Size;k++) {
    Elm += Pt1[0]*Pt2[0]*Pt3[0];
    Pt1++;
    Pt2++;
    Pt3++;
  }
  return Elm;
}


int DVR::GetElmDiagOperator(int i, int j,double* Ai, int rank, double* Op) {
  if (Ai == 0) { return EXIT_FAILURE; }
  if (Op == 0) { return EXIT_FAILURE; }
  if (EVectors == 0) { return EXIT_FAILURE; }
  for(int l=0;l<rank;l++) { Op[l] = 0.; }
  double* Pt1 = EVectors + i*Size;
  double* Pt2 = EVectors + j*Size; 
  double* Pt3 = Ai;
  double* Opp;
  for(int k=0;k<Size;k++) {
    Opp = Op;
    for(int l=0;l<rank;l++) {
      Opp[0] += Pt1[0]*Pt2[0]*Pt3[0];
      Opp++;
      Pt3++;
    }
    Pt1++;
    Pt2++;
  }
  return EXIT_SUCCESS;
}

