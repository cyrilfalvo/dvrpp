/*
 *  BasisSet.h
 *  DVR
 *
 *  Created by Cyril Falvo on 02/12/10.
 *  Copyright 2010 Université Paris-Sud. All rights reserved.
 *
 */

#ifndef _CLASS_BASISSET_
#define _CLASS_BASISSET_

#define SUCCESS  0
#define FAILURE  1

#ifndef PI
#define PI              3.141592653589793
#endif

class BasisSet {

public:
  BasisSet();
  ~BasisSet();
  
  void SetFact(double FactV);
  int ComputeNodes();
  
  virtual int ComputeKineticOperator();
  virtual double BasisReal(double x);
  virtual int SetBasisSet(int Size, double DeltaV, double X0V);
  
  double* Xi;
  double* Tq;  

 protected:

  int N;
  double* Basis;
  double Delta;
  double X0;
  double Fact;
  
};

#endif


