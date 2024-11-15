
#include "../include/DVRpp.h"
#include <cstdio>

extern "C" {
  
  // C wrapper for C++ DVR++ contructor
  DVR* DVRpp_create() {
    return new DVR();
  }

  // C wrapper for C++ DVR+++ destructor
  void DVRpp_destroy(DVR* dvr) {
    if(dvr != 0 )  {
      delete dvr;
    }
  }

  void DVRpp_SetDimension(DVR* dvr, int DimV) { 
    dvr->SetDimension(DimV);
  }
  
  void DVRpp_AddHermBasis(DVR* dvr, int mode, int N, double Fact, double Delta, double Xeq, int info) {
    info = dvr->AddHermBasis(mode,N,Fact,Delta,Xeq);
  }
  
  void DVRpp_AddSineBasis(DVR* dvr, int mode, int N, double Fact, double Delta, double Xeq, int info) {
    info = dvr->AddSineBasis(mode,N,Fact,Delta,Xeq);
  }
  
  void DVRpp_AddExpBasis(DVR* dvr, int mode, int N, double Fact, double Delta, double Xeq, int info) {
    info = dvr->AddExpBasis(mode,N,Fact,Delta,Xeq);
  }
  
  void DVRpp_PrepareDvr(DVR* dvr, int info) {
    info = dvr->PrepareDvr();
  }

  void DVRpp_GetNodeCoord(DVR* dvr, int index, double* R, int dim) {
    dvr->GetNodeCoord(index, R);
  }
  
  void DVRpp_LoadPotential_array(DVR* dvr, double* V) {
    dvr->LoadPotential(V);
  }
  
  void DVRpp_LoadPotential(DVR* dvr, int index, double V) {
    dvr->LoadPotential(index,V);
  }

#ifdef _ARPACK_
  void DVRpp_SolveHamiltonianArpack(DVR* dvr, int nev, int ncv, int info) {
    info = dvr->SolveHamiltonianArpack(nev,ncv);
  }
#endif
  
  void DVRpp_SolveHamiltonianLapack(DVR* dvr, int info) {
    info = dvr->SolveHamiltonianLapack();
  }

  double DVRpp_GetEnergy(DVR* dvr, int i) {
    return dvr->GetEnergy(i);
  }

  double DVRpp_GetElmDiagOperatorSingle(DVR* dvr, int i, int j, double* Ai) {
    return dvr->GetElmDiagOperatorSingle(i,j,Ai);
  }

  void DVRpp_GetElmDiagOperator(DVR* dvr, int i, int j, double* Ai, int rank, double* Op, int info) {
    info = dvr->GetElmDiagOperator(i,j,Ai,rank,Op);
  }
  
  void DVRpp_GetElmQOpt(DVR* dvr, int PowMax, int i, double* AvgOpt, int info) {
    info = dvr->GetElmQOpt(PowMax, i, AvgOpt);
  }
  
  void DVRpp_GetPop(DVR* dvr, int i, int* ni, int info) {
    info = dvr->GetPop(i,ni);
  }
  
  void DVRpp_GetVector(DVR* dvr, int i, double* Psi, int info) {
    info = dvr->GetVector(i,Psi);
  }

  int DVRpp_GetSize(DVR* dvr) {
    return dvr->GetSize();
  }

}



