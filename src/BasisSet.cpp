/*
 *  BasisSet.cpp
 *  DVR
 *
 *  Created by Cyril Falvo on 02/12/10.
 *  Copyright 2010 Université Paris-Sud. All rights reserved.
 *
 */

#include "../include/BasisSet.h"
#include <cstdio>
#include "../include/lapack.h"

BasisSet::BasisSet() {
  N     = 0;
  Basis = 0;
  Xi    = 0;
  Delta = 0;
  X0    = 0;
  Tq    = 0;
}

BasisSet::~BasisSet() {
  if(Basis!=0) { delete [] Basis; } 
  if(Xi!=0) { delete [] Xi; }
  if(Tq!=0) { delete [] Tq; }
  Basis = 0;
  Xi    = 0;
  N     = 0;
  Delta = 0;
  X0    = 0;
}

int BasisSet::SetBasisSet(int Size, double DeltaV,double X0V) {
  N=Size;
  Delta = DeltaV;
  X0 = X0V;
  Delta = DeltaV;
  X0 = X0V;
  if(N<=0) { return FAILURE; }
  Basis = new double[N*N];
  Tq    = new double[N*N]; 
  Xi    = new double[N];
  for(int i=0;i<N*N;i++) { Basis[i] = 0.; }
  for(int i=0;i<N*N;i++) { Tq[i] = 0.; }
  for(int i=0;i<N;i++) { Xi[i] = 0.; }
  return SUCCESS;
}

int BasisSet::ComputeNodes() {
  
  int lwork;
  int info;
  double* work;
  
  lwork=-1;
  work = new double[2];  
  dsyev_('V','U',N,Basis,N,Xi,work,lwork,info);
  lwork = (int) work[0];
  delete [] work;
  work = new double[lwork];
  dsyev_('V','U',N,Basis,N,Xi,work,lwork,info);
  delete [] work;
	
  for(int i=0;i<N;i++) {
    if(Basis[0+i*N] < 0) {
      for(int j=0;j<N;j++) {
	Basis[j+i*N] *= -1;
      }
    }
  }
  return SUCCESS;
}

void BasisSet::SetFact(double FactV) {
  Fact = FactV;
}


int BasisSet::ComputeKineticOperator() {
  return SUCCESS;
}

double BasisSet::BasisReal(double x) {
  return SUCCESS;
}


