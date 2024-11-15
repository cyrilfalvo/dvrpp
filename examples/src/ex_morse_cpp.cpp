
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "../../include/DVRpp.h"

using namespace std;

#define alpha 0.3

double morse(double q) {
  double PP;
  double fact;
  fact = 1.-exp(-alpha*q);
  PP = 0.5/alpha/alpha*fact*fact;
  return PP;
}

double analytic(int n) {
  double D;
  double xe;
  double energy;
  double nn;

  D = 0.5/alpha/alpha;
  xe = 0.25/D;
 
  nn = (double)n+0.5;
  energy = nn*(1. - xe*nn);
  
  return energy;
}


int main (int argc, char * argv[]) {

  DVR*  Dvr = new DVR();

  int    np = 20;
  int    size;
  
  Dvr->SetDimension(1);
  Dvr->AddHermBasis(0,np,1.,1.2,2.50);
  Dvr->PrepareDvr();
 
  size = Dvr->GetSize();
 
  for(int i=0;i<size;i++) {
    double QQ;
    Dvr->GetNodeCoord(i,&QQ);
    double Potential = morse(QQ);
    Dvr->LoadPotential(i,Potential);
  }
  
  Dvr->SolveHamiltonianLapack();

  for(int i=0;i<10;i++) {
    printf("%3d %16.4f %16.4f %16.4f\n",i,analytic(i),Dvr->GetEnergy(i),(Dvr->GetEnergy(i)-analytic(i)));
  }

  return EXIT_SUCCESS;
  
}
