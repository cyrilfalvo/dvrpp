/*
 *  BasisHerm.cpp
 *  DVR
 *
 *  Created by Cyril Falvo on 02/12/10.
 *  Copyright 2010 Université Paris-Sud. All rights reserved.
 *
 */

#include "../include/BasisHerm.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

BasisHerm::BasisHerm() {
}

BasisHerm::~BasisHerm() {
}

int BasisHerm::SetBasisSet(int Size, double DeltaV, double X0V){
  BasisSet::SetBasisSet(Size,DeltaV,X0V);
  for(int i=0;i<N-1;i++) {
    Basis[i+N*i] = X0;
    Basis[i+N*(i+1)] = sqrt(double(i+1)/2.)*Delta;
    Basis[(i+1)+N*i] = sqrt(double(i+1)/2.)*Delta;
  }
  Basis[(N-1)*(N+1)] = X0;
  ComputeNodes();
  ComputeKineticOperator();
  return SUCCESS;
}

int BasisHerm::ComputeKineticOperator() {
  if (A     == 0) { return FAILURE; }
  if (Tq    == 0) { return FAILURE; }	
  if (Delta == 0) { return FAILURE; }
  if (N     == 0) { return FAILURE; }
  if (Xi    == 0) { return FAILURE; }
  
  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      Tq[i+j*N] = 0;
      for(int k=0;k<N;k++) {
	Tq[i+j*N] += Basis[k+i*N]*Basis[k+j*N]*double(2*k+1); 
      }
      Tq[i+j*N] *= 0.5/pow(Delta,2)*A;
    }
    Tq[i+i*N] += -0.5*pow(Xi[i]-X0,2)/pow(Delta,4)*A;
  }
  return SUCCESS;
}



