/*
 *  BasisSine.cpp
 *  DVR
 *
 *  Created by Cyril Falvo on 13/19/12.
 *  Copyright 2012 Université Paris-Sud. All rights reserved.
 *
 */

#include "../include/BasisExp.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

BasisExp::BasisExp() {
  L=0;
}

BasisExp::~BasisExp() {
  L=0;
}

int BasisExp::SetBasisSet(int Size, double DeltaV, double X0V){
  BasisSet::SetBasisSet(Size,DeltaV,X0V);
  L = Delta*(double)N;
  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      Basis[i+j*N] = sqrt(1./L)*cos((double)i*(double)j*PI/(double)N);
    }
  }
  for(int i=0;i<N;i++) {
    Xi[i] = X0 + (double)i*Delta;
  }
  ComputeKineticOperator();

  return SUCCESS;
}

int BasisExp::ComputeKineticOperator() {
  if (A     == 0) { return FAILURE; }
  if (Tq    == 0) { return FAILURE; }	
  if (Delta == 0) { return FAILURE; }
  if (N     == 0) { return FAILURE; }
  if (Xi    == 0) { return FAILURE; }
  
  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      if (i==j) {
        Tq[i+i*N] = PI*PI*A/6./L/L*((double)N*(double)N-1);
      }
      else {
        Tq[i+j*N] = PI*PI*A/L/L*pow(-1.,i-j)*cos(PI*(double)(i-j)/(double)N)/pow( sin(PI*(double)(i-j)/(double)N),2);
      }
    }
  }
  return SUCCESS;
}



