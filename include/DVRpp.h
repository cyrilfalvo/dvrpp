
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
  int  AddHermBasis(int mode,int N,double Fact, double Delta, double Xeq);
  int  AddSineBasis(int mode,int N,double Fact, double Delta, double Xeq);
  int  AddExpBasis(int mode,int N,double Fact, double Delta, double Xeq);
  int  PrepareDvr();
  void  GetNodeCoord(int index, double* R);
  void  LoadPotential(double* V);
  void  LoadPotential(int index, double V);
	
#ifdef _ARPACK_
  int SolveHamiltonianArpack(int nev, int ncv);
#endif

  int SolveHamiltonianLapack();
	
  double GetEnergy(int i);
  double GetElmDiagOperatorSingle(int i, int j, double* Ai);
  int    GetElmDiagOperator(int i, int j, double* Ai, int rank,double* Op);
  int    GetElmQOpt(int PowMax,int i, double* AvgOpt);
  int    GetPop(int i,int* ni);
  int    GetVector(int i, double* Psi);
  int    GetSize();
	
private:

	BasisSet** Basis;
	
	int Dim;
        int Size;
	int* NGrid;
	int MaxGrid;
	
	double* Vi;
  
	double* EVectors;
	double* EEnergies;
	
    	void MultHamiltonian(double* d1, double* d2);
	
	int* Pop;
	int* Fact;

	int flag;   // this flag is used to check that everything is ok 

};    

#endif
