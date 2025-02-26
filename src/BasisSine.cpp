/*
 *  BasisSine.cpp
 *  DVR
 *
 *  Created by Cyril Falvo on 11/10/11.
 *  Copyright 2010 Université Paris-Sud. All rights reserved.
 *
 */

#include "../include/BasisSine.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

BasisSine::BasisSine() {
  L=0;
  M=0;
}

BasisSine::~BasisSine() {
  L=0;
  M=0;
}

int BasisSine::SetBasisSet(int Size, double DeltaV, double X0V){
  BasisSet::SetBasisSet(Size,DeltaV,X0V);
  M = N+1;
  L = Delta*(double)M;
  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      Basis[i+j*N] = sqrt(2./L)*sin(double((i+1)*(j+1))*PI/double(M));
    }
  }
  for(int i=0;i<N;i++) {
    Xi[i] = X0 + (double)i*Delta;
  }
  ComputeKineticOperator();

  return SUCCESS;
}

int BasisSine::ComputeKineticOperator() {
  if (A     == 0) { return FAILURE; }
  if (Tq    == 0) { return FAILURE; }	
  if (Delta == 0) { return FAILURE; }
  if (N     == 0) { return FAILURE; }
  if (Xi    == 0) { return FAILURE; }
  
  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      if (i==j) {
        Tq[i+i*N] = 0.25*PI*PI*A/L/L*
	  (double(2*M*M+1)/3. - 1./pow(sin(PI*double(i+1)/double(M)),2));
      }
      else {
        Tq[i+j*N] = 0.25*PI*PI*A/L/L*double(pow(-1.,i-j))*
	  (1./pow(sin(PI*double(i-j)/double(2*M)),2) - 1./pow(sin(PI*double(i+j+2)/double(2*M)),2) );
      }
    }
  }
  return SUCCESS;
}



