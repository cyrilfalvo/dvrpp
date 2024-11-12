
#include "../include/DVR++.h"

extern "C" {
  
  // C wrapper for C++ DVR++ contructor
  DVR* dvrpp_create() {
    return new DVR();
  }

  // C wrapper for C++ DVR+++ destructor
  void dvrpp_destroy(DVR* dvr) {
    delete dvr;
  }

  void dvrpp_SetDimension(DVR* dvr, int DimV) { 
    dvr->SetDimension(DimV);
  }
  
  int dvrpp_AddHermBasis(DVR* dvr, int mode, int N, double Fact, double Delta, double Xeq) {
    return dvr->AddHermBasis(mode,N,Fact,Delta,Xeq);
  }
  
  int dvrpp_AddSineBasis(DVR* dvr, int mode, int N, double Fact, double Delta, double Xeq) {
    return dvr->AddSineBasis(mode,N,Fact,Delta,Xeq);
  }
  
  int dvrpp_AddExpBasis(DVR* dvr, int mode, int N, double Fact, double Delta, double Xeq) {
    return dvr->AddExpBasis(mode,N,Fact,Delta,Xeq);
  }
  
  int dvrpp_PrepareDvr(DVR* dvr) {
    return dvr->PrepareDvr();
  }

  int dvrpp_GetNodeCoord(DVR* dvr, int index, double* R) {
    return dvr->GetNodeCoord(index, R);
  }
  

  int dvrpp_LoadPotential_array(DVR* dvr, double* V) {
    return dvr->LoadPotential(V);
  }
  
  int dvrpp_LoadPotential(DVR* dvr, int index, double V) {
    return dvr->LoadPotential(index,V);
  }
  


}



