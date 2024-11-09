/*
 *  BasisHerm.h
 *  DVR
 *
 *  Created by Cyril Falvo on 02/12/10.
 *  Copyright 2010 Université Paris-Sud. All rights reserved.
 *
 */

#ifndef _CLASS_BASISEXP_
#define _CLASS_BASISEXP_

#define SUCCESS  0
#define FAILURE  1

#include "BasisSet.h"

class BasisExp : public BasisSet {
	
public:
  BasisExp();
  ~BasisExp();

  double L;

  int SetBasisSet(int Size, double DeltaV, double X0V);
  int ComputeKineticOperator();
  
};

#endif


