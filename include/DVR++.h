
//###########################################//
//#                                         #//
//# Class to compute DVR of a nD potential  #//
//#                                         #//
//# units in atomic units                   #//
//#                                         #//
//###########################################// 


#ifndef _CLASS_DVR_
#define _CLASS_DVR_

#include "BasisSet.h"

class DVR {

public :

  DVR();
  ~DVR();

  void SetDimension(int DimV);
  int  AddHermBasis(int mode,int N,double Freq, double Delta, double Xeq);
  int  AddSineBasis(int mode,int N,double Freq, double Delta, double Xeq);
  int  AddExpBasis(int mode,int N,double Freq, double Delta, double Xeq);
  int  PrepareDvr();
  int  GetNodeCoord(int index, double* R);
  int  LoadPotential(double* V);
  int  LoadPotential(int index, double V);
  int  LoadDiagOperator(int index, double A, int rank);
	
#ifdef _ARPACK_
  int SolveHamiltonianArpack(int nev, int ncv);
#endif

  int SolveHamiltonianLapack();
	
  double GetEnergy(int i);
  double GetElmDiagOperatorSingle(int i, int j, double* Ai);
  int    GetElmDiagOperator(int i, int j, double* Ai, int rank,double* Op);
  int    GetElmQOpt(int PowMax,int i, double* AvgOpt);
  int    GetPop(int i,int* ni);
  int    Size;
  int    GetVector(int i, double* Psi);
	
private:

	BasisSet** Basis;
	
	int Dim;
	int* NGrid;
	int MaxGrid;
	
	double* Vi;
  
	double* EVectors;
	double* EEnergies;

	int MultHamiltonian(double* d1, double* d2);
	
	int* Pop;
	int* Fact;
};    

#endif
